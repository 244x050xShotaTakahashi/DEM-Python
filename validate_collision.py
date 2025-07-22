import numpy as np
import pandas as pd

def calculate_theoretical_velocities(v1_initial, v2_initial, e):
    """
    理論的な衝突後の速度を計算する関数
    
    Parameters:
    -----------
    v1_initial : float
        粒子1の初期速度
    v2_initial : float
        粒子2の初期速度
    e : float
        反発係数
    
    Returns:
    --------
    tuple
        (v1_final, v2_final) 衝突後の速度
    """
    v1_final = ((1 - e) * v1_initial + (1 + e) * v2_initial) / 2
    v2_final = ((1 + e) * v1_initial + (1 - e) * v2_initial) / 2
    return v1_final, v2_final

def calculate_relative_error(actual, theoretical):
    """
    相対誤差を計算する関数
    
    Parameters:
    -----------
    actual : float
        実際の値
    theoretical : float
        理論値
    
    Returns:
    --------
    float
        相対誤差（%）
    """
    return abs(actual - theoretical) / abs(theoretical) * 100

def validate_collision(file_path, v1_initial, v2_initial, e):
    """
    衝突計算の妥当性を検証する関数
    
    Parameters:
    -----------
    file_path : str
        結果ファイルのパス
    v1_initial : float
        粒子1の初期速度
    v2_initial : float
        粒子2の初期速度
    e : float
        反発係数
    """
    # ファイルの読み込み
    df = pd.read_csv(file_path)
    
    # 最後のステップのデータを取得
    last_step = df['Step'].max()
    final_data = df[df['Step'] == last_step]
    
    # 粒子の速度を取得
    v1_final = final_data[final_data['Particle'] == 0]['Vx'].values[0]
    v2_final = final_data[final_data['Particle'] == 1]['Vx'].values[0]
    
    # 理論値を計算
    v1_theoretical, v2_theoretical = calculate_theoretical_velocities(v1_initial, v2_initial, e)
    
    # 相対誤差を計算
    error_v1 = calculate_relative_error(v1_final, v1_theoretical)
    error_v2 = calculate_relative_error(v2_final, v2_theoretical)
    
    # 運動量保存の確認
    initial_momentum = v1_initial + v2_initial
    final_momentum = v1_final + v2_final
    momentum_error = calculate_relative_error(final_momentum, initial_momentum)
    
    # 実効的な反発係数の計算
    e_effective = (v2_final - v1_final) / (v1_initial - v2_initial)
    e_error = calculate_relative_error(e_effective, e)
    
    # 結果の表示
    print(f"衝突計算の妥当性検証結果:")
    print(f"初期条件:")
    print(f"  粒子1の初期速度: {v1_initial}")
    print(f"  粒子2の初期速度: {v2_initial}")
    print(f"  反発係数: {e}")
    print(f"\nシミュレーション結果:")
    print(f"  粒子1の最終速度: {v1_final:.6f}")
    print(f"  粒子2の最終速度: {v2_final:.6f}")
    print(f"\n理論値:")
    print(f"  粒子1の理論速度: {v1_theoretical:.6f}")
    print(f"  粒子2の理論速度: {v2_theoretical:.6f}")
    print(f"\n相対誤差:")
    print(f"  粒子1: {error_v1:.6f}%")
    print(f"  粒子2: {error_v2:.6f}%")
    print(f"  運動量: {momentum_error:.4f}%")
    print(f"  反発係数: {e_error:.2f}%")
    print(f"\n実効的な反発係数: {e_effective:.4f}")

if __name__ == "__main__":
    # 検証パラメータの設定
    file_path = "/Users/shota.takahashi/Desktop/cursor_folder/DEM/Python/outputs/dem_results.csv.230"
    v1_initial = 20.0
    v2_initial = 0.0
    e = 0.9999313555136594
    
    # 妥当性検証の実行
    validate_collision(file_path, v1_initial, v2_initial, e) 