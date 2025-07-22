import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
import glob
import os
import json
import math

def numerical_sort(value: str) -> int:
    """
    'outputs/dem_results.csv.12345' -> 12345
    末尾に余計な空白やコピー番号が付いていても無視する。
    """
    tail = value.split('.')[-1].strip()      # 末尾を取り出し空白除去
    tail = tail.split()[0]                   # 空白以降を捨てる
    return int(tail) if tail.isdigit() else 0

def create_dem_animation():
    # 設定ファイルの読み込み
    config_path = './config.json'
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # 粒子の半径をconfig.jsonから取得
    PARTICLE_RADIUS = config['particles']['radius']
    
    # ファイル名を数値順にソート
    files = sorted(glob.glob('./outputs/dem_results.csv.*'), key=numerical_sort)
    
    if not files:
        print("エラー: outputsディレクトリにdem_results.csv.*ファイルが見つかりません")
        return
    
    # フレーム数を制限する処理を追加
    frame_interval = max(1, len(files) // 100)  # 約1000フレームになるように間引く
    files = files[::frame_interval]
    
    total_frames = len(files)
    print(f"アニメーション作成開始: 合計{total_frames}フレーム（間引き率: {frame_interval}）")
    
    df = pd.read_csv(files[0])
    
    # デバッグ用：利用可能な列名を表示
    # print("利用可能な列名:", df.columns.tolist())
    
    fig, ax = plt.subplots(figsize=(10, 10))
    # config.jsonの値に基づいて表示範囲を設定
    x_max = max(max(line['x1'], line['x2']) for line in config['lines'])
    x_min = min(min(line['x1'], line['x2']) for line in config['lines'])
    y_max = max(max(line['y1'], line['y2']) for line in config['lines'])
    y_min = min(min(line['y1'], line['y2']) for line in config['lines'])
    ax.set_xlabel('X', fontsize=20)
    ax.set_ylabel('Y', fontsize=20)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_aspect('equal')

    # 粒子を円として描画するための初期化
    particles = []
    particle_numbers = []  # 粒子番号表示用のテキストオブジェクトを格納
    for i in range(len(df)):
        circle = patches.Circle((0, 0), radius=PARTICLE_RADIUS,
                            fc='blue', alpha=0.6)
        line = plt.Line2D([0, 0], [0, 0], color='black')
        number = ax.text(0, 0, str(i+1), ha='center', va='center', color='yellow', fontsize=10)
        particles.append((circle, line))
        particle_numbers.append(number)
        ax.add_patch(circle)
        ax.add_line(line)
    
    # 壁（線）の描画
    for line in config['lines']:
        x = [line['x1'], line['x2']]
        y = [line['y1'], line['y2']]
        ax.plot(x, y, 'k-', linewidth=1)
    
    # 目盛りのフォントサイズを設定
    ax.tick_params(axis='both', which='major', labelsize=20)
    title = ax.set_title('DEM Simulation (Step: 0)', fontsize=24)

    def init():
        for circle, line in particles:
            circle.center = (0, 0)
            line.set_data([0, 0], [0, 0])
        for number in particle_numbers:
            number.set_position((0, 0))
        return [p[0] for p in particles] + [p[1] for p in particles] + particle_numbers + [title]

    
    def update(frame):
        # 進捗状況を表示（パーセンテージ付き）
        progress = (frame + 1) / total_frames * 100
        print(f"\r処理中... {frame + 1}/{total_frames} フレーム完了 ({progress:.1f}%)", end="")
        
        df = pd.read_csv(files[frame])
        
        # 各粒子の位置と回転を更新
        for i, ((circle, line), number) in enumerate(zip(particles, particle_numbers)):
            # 文字列を浮動小数点数に変換
            x = float(df['X'].iloc[i])   
            y = float(df['Y'].iloc[i])  
            angle = float(df['Angle'].iloc[i]) if 'Angle' in df.columns else 0
            
            circle.center = (x, y)
            number.set_position((x, y))  # 粒子番号の位置を更新
            
            # 回転を示す線の更新
            line_x = [x, x + PARTICLE_RADIUS * math.cos(angle)]
            line_y = [y, y + PARTICLE_RADIUS * math.sin(angle)]
            line.set_data(line_x, line_y)
        
        step = int(df['Step'].iloc[0]) if 'Step' in df.columns else frame
        title.set_text(f'DEM Simulation (Step: {step})')
        title.set_fontsize(24)
        
        return [p[0] for p in particles] + [p[1] for p in particles] + particle_numbers + [title]
    
    anim = animation.FuncAnimation(
        fig, update, init_func=init,
        frames=len(files), interval=100, blit=True
    )
    
    # 保存設定を変更
    print("\nGIFファイルを保存中...")
    writer = animation.PillowWriter(fps=30, bitrate=1800)
    anim.save('dem_simulation.gif', writer=writer, dpi=100)  # dpiを下げてファイルサイズを削減
    print("\nアニメーション作成完了: dem_simulation.gif")
    plt.close()

if __name__ == "__main__":
    create_dem_animation()