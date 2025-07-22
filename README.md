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

## オーバーラップと速度の妥当性検証

このプロジェクトには、数値計算結果と理論解を比較してシミュレーションの精度を検証する機能が含まれています。

### 概要

粒子衝突時のオーバーラップ量と速度の時系列データを理論解（減衰自由振動の解析解）と比較し、以下の誤差指標を計算します：

- **MAE (Mean Absolute Error)**: 平均絶対誤差
- **RMSE (Root Mean Square Error)**: 二乗平均平方根誤差
- **MAX**: 最大誤差
- **REL_PERCENT**: 平均相対誤差（%）

### 前提条件

妥当性検証を実行する前に、DEMシミュレーションを実行して`outputs/overlap/`ディレクトリに`overlap_velocity_comparison_*.csv`ファイルが生成されている必要があります。

### 基本的な使用方法

```bash
# シミュレーション実行（粒子衝突が発生する設定で）
python dem_calc.py

# 妥当性検証の実行（基本）
python check_overlap_validity.py -f outputs/overlap/overlap_velocity_comparison_0_1_12345.csv

# 許容誤差を指定して実行
python check_overlap_validity.py -f outputs/overlap/overlap_velocity_comparison_0_1_12345.csv --tol 1e-4

# グラフを表示しながら実行
python check_overlap_validity.py -f outputs/overlap/overlap_velocity_comparison_0_1_12345.csv --plot

# グラフをPNGファイルに保存
python check_overlap_validity.py -f outputs/overlap/overlap_velocity_comparison_0_1_12345.csv --png result_validation.png
```

### オプション詳細

| オプション | 説明 | デフォルト |
|-----------|------|-----------|
| `-f`, `--file` | 解析対象のCSVファイルパス（必須） | - |
| `--tol` | 許容誤差しきい値 | 1e-3 |
| `--plot` | グラフをウィンドウで表示 | False |
| `--png` | グラフをPNGファイルに保存（ファイル名指定） | None |

### 出力例

```
=== オーバーラップの誤差指標 ===
MAE         :  1.234567e-05
RMSE        :  2.345678e-05
MAX         :  5.678901e-05
REL_PERCENT :  0.123456 %

=== 速度の誤差指標 ===
MAE         :  3.456789e-04
RMSE        :  4.567890e-04
MAX         :  8.901234e-04
REL_PERCENT :  0.234567 %

→ オーバーラップ判定: OK (許容誤差 1e-3 以下)
→ 速度判定: OK (許容誤差 1e-3 以下)
```

### パラメータ掃引での自動実行

`parameter_sweep.py`を実行すると、各パラメータ組み合わせでシミュレーションを実行後、自動的に妥当性検証が実行され、結果がPNGファイルに保存されます：

```bash
python parameter_sweep.py
```

この場合、`outputs/overlap/`ディレクトリ内に`e=0.5_v=20_k=100.png`のような名前でグラフが保存されます。

### 判定基準

- **OK**: MAEとRMSEの両方が指定した許容誤差（`--tol`）未満
- **NG**: MAEまたはRMSEのいずれかが許容誤差を超過

この機能により、シミュレーションパラメータの妥当性や数値積分手法の精度を定量的に評価できます。

## 開発者

Shota Takahashi (244x050xShotaTakahashi)
