# DEM (Discrete Element Method) Python Simulation

離散要素法（DEM）による粒子シミュレーションのPythonプロジェクトです。

## 概要

このプロジェクトは、粒子間の衝突や接触力を数値計算で解析するDEMシミュレーションツールです。物理的に正確な粒子の運動をシミュレートし、衝突や壁との相互作用を可視化できます。

## 主要機能

- **粒子間衝突シミュレーション**: 多粒子系の衝突を物理的に正確にシミュレート
- **壁との相互作用**: 粒子と壁面との衝突を計算
- **可視化機能**: アニメーション形式でシミュレーション結果を表示
- **データ出力**: CSV形式でシミュレーション結果を保存
- **跳ね返り係数解析**: 理論値と数値計算結果の比較分析

## ファイル構成

- `dem_calc.py`: メインのDEMシミュレーションクラス
- `config.json`: シミュレーションパラメータ設定ファイル
- `main.py`: シミュレーション実行用メインスクリプト
- `visualizer.py`: アニメーション作成機能
- `elements.py`: 粒子と線要素のクラス定義
- `interface.py`: 接触界面の物理パラメータ定義

## 必要な環境

- Python 3.6以上
- numpy
- matplotlib
- tkinter

## インストール

```bash
git clone https://github.com/244x050xShotaTakahashi/DEM-Python.git
cd DEM-Python
pip install numpy matplotlib
```

## 使用方法

### 基本的な実行方法

```bash
python main.py
```

### 設定の変更

`config.json`ファイルを編集してシミュレーションパラメータを調整できます：

```json
{
  "simulation": {
    "dt": 1e-05,
    "max_steps": 10000,
    "output_interval": 1000
  },
  "particles": {
    "radius": 5.0,
    "count": 2,
    "random": false
  }
}
```

## 出力ファイル

シミュレーション実行後、以下のファイルが出力されます：

- `outputs/dem_results.csv.*`: 各ステップでの粒子状態
- `outputs/collision_statistics.csv`: 衝突統計情報
- `outputs/overlap/`: オーバーラップ解析結果
- `dem_simulation.gif`: アニメーション（可視化機能使用時）

## 開発者

Shota Takahashi (244x050xShotaTakahashi)

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。
