import sys
import time
import json
from config import Config
from dem_simulation import DEMSimulation
from visualization import Visualizer

def main():
    # 設定ファイルの読み込み
    try:
        config = Config('config.json')
    except FileNotFoundError:
        print("エラー: config.jsonファイルが見つかりません。")
        sys.exit(1)
    except json.JSONDecodeError:
        print("エラー: config.jsonの形式が正しくありません。")
        sys.exit(1)

    # シミュレーションの初期化
    simulation = DEMSimulation(config)

    # 可視化の初期化
    visualizer = Visualizer(simulation)

    # シミュレーションのメインループ
    try:
        while True:
            start_time = time.time()

            # シミュレーションステップの実行
            simulation.step()

            # 可視化の更新
            visualizer.redraw()

            # フレームレート制御
            elapsed_time = time.time() - start_time
            target_frame_time = 1.0 / config.get('visualization', {}).get('fps', 60)
            if elapsed_time < target_frame_time:
                time.sleep(target_frame_time - elapsed_time)

            # 経過時間の表示（オプション）
            if config.get('debug', {}).get('show_elapsed_time', False):
                print(f"Step time: {elapsed_time:.4f} seconds")

    except KeyboardInterrupt:
        print("\nシミュレーションが中断されました。")
    finally:
        # クリーンアップ処理（必要に応じて）
        visualizer.close()

if __name__ == '__main__':
    main()