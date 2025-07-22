#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
overlap_comparison_*.csv を解析し，数値解と理論解の一致度をチェックするスクリプト
使い方:
    python check_overlap_validity.py -f outputs/overlap_comparison_0_1_201606.csv \
                                    --tol 1e-3 --plot
"""

import argparse
import csv
import math
import os
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np


def read_overlap_csv(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """CSV を読み込み，time, sim_overlap, theo_overlap, sim_velocity, theo_velocity の配列を返す"""
    times: List[float] = []
    sim_overlap: List[float] = []
    theo_overlap: List[float] = []
    sim_velocity: List[float] = []
    theo_velocity: List[float] = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)  # 新しいフォーマット：['time', 'sim_overlap', 'theoretical_overlap', 'overlap_abs_error', 'overlap_rel_error(%)', 'sim_velocity', 'theoretical_velocity', 'velocity_abs_error', 'velocity_rel_error(%)']
        for row in reader:
            if len(row) < 7:  # 最低限必要な列数をチェック
                continue
            times.append(float(row[0]))
            sim_overlap.append(float(row[1]))
            theo_overlap.append(float(row[2]))
            sim_velocity.append(float(row[5]))
            theo_velocity.append(float(row[6]))
    return np.array(times), np.array(sim_overlap), np.array(theo_overlap), np.array(sim_velocity), np.array(theo_velocity)


def compute_metrics(sim: np.ndarray, theo: np.ndarray) -> dict:
    """MAE, RMSE, 最大誤差, 相対誤差(%) を返す"""
    abs_err = np.abs(sim - theo)
    mae = abs_err.mean() # 平均絶対誤差
    rmse = math.sqrt(((sim - theo) ** 2).mean()) # 二乗平均平方根誤差
    max_err = abs_err.max() # 最大誤差
    # 理論解が 0 に近い箇所があるのでゼロ除算を回避しつつ相対誤差を計算
    rel_err = (
        (abs_err / np.where(np.abs(theo) > 1e-12, np.abs(theo), 1.0)).mean() * 100.0 # 相対誤差
    )
    return dict(MAE=mae, RMSE=rmse, MAX=max_err, REL_PERCENT=rel_err) # 誤差指標を返す


def judge(metrics: dict, tol: float) -> bool:
    """MAE と RMSE がともに tol 未満なら OK"""
    return metrics["MAE"] < tol and metrics["RMSE"] < tol # MAE と RMSE がともに tol 未満なら OK


def main():
    parser = argparse.ArgumentParser(description="オーバーラップ妥当性チェック")
    parser.add_argument(
        "-f",
        "--file",
        required=True,
        help="解析対象 CSV (outputs/overlap_comparison_*.csv)",
    )
    parser.add_argument(
        "--tol", type=float, default=1e-3, help="許容誤差しきい値 (デフォルト 1e-3)"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        help="グラフを表示 （matplotlib が必要）",
    )
    parser.add_argument(
        "--png",
        type=str,
        default=None,
        help="グラフを PNG 形式で保存するファイル名 (拡張子 .png を含む)"
    )
    args = parser.parse_args()

    if not os.path.exists(args.file):
        print(f"ファイルが見つかりません: {args.file}")
        return

    t, sim_overlap, theo_overlap, sim_velocity, theo_velocity = read_overlap_csv(args.file)
    metrics_overlap = compute_metrics(sim_overlap, theo_overlap)
    metrics_velocity = compute_metrics(sim_velocity, theo_velocity)

    print("=== オーバーラップの誤差指標 ===")
    for k, v in metrics_overlap.items(): # 誤差指標を表示
        if k == "REL_PERCENT": # 相対誤差は小数点以下6桁まで表示
            print(f"{k:12s}: {v: .6f} %") # 小数点以下6桁まで表示
        else:
            print(f"{k:12s}: {v: .6e}") # 指数表示
    
    print("\n=== 速度の誤差指標 ===")
    for k, v in metrics_velocity.items(): # 誤差指標を表示
        if k == "REL_PERCENT": # 相対誤差は小数点以下6桁まで表示
            print(f"{k:12s}: {v: .6f} %") # 小数点以下6桁まで表示
        else:
            print(f"{k:12s}: {v: .6e}") # 指数表示

    # オーバーラップの判定をメインとする
    if judge(metrics_overlap, args.tol):
        print(f"\n→ オーバーラップ判定: OK (許容誤差 {args.tol:g} 以下)")
    else:
        print(f"\n→ オーバーラップ判定: NG (許容誤差 {args.tol:g} を超過)")
    
    if judge(metrics_velocity, args.tol):
        print(f"→ 速度判定: OK (許容誤差 {args.tol:g} 以下)")
    else:
        print(f"→ 速度判定: NG (許容誤差 {args.tol:g} を超過)")

    # プロットを作成する条件: --plot または --png が指定された場合
    if args.plot or args.png:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
        
        # オーバーラップのプロット
        ax1.plot(t, sim_overlap, label="Simulation", lw=1)
        ax1.plot(t, theo_overlap, label="Theory", lw=1, ls="--")
        ax1.set_xlabel("Time [s]")
        ax1.set_ylabel("Overlap [m]")
        ax1.set_title("Overlap: Simulation vs Theory")
        ax1.legend()
        ax1.grid(True)
        
        # 速度のプロット
        ax2.plot(t, sim_velocity, label="Simulation", lw=1)
        ax2.plot(t, theo_velocity, label="Theory", lw=1, ls="--")
        ax2.set_xlabel("Time [s]")
        ax2.set_ylabel("Velocity [m/s]")
        ax2.set_title("Velocity: Simulation vs Theory")
        ax2.legend()
        ax2.grid(True)
        
        plt.tight_layout()

        # PNG 保存
        if args.png:
            plt.savefig(args.png, dpi=150)
            print(f"グラフを保存しました: {args.png}")

        # 表示
        if args.plot:
            plt.show()
        else:
            plt.close()


if __name__ == "__main__":
    main()
