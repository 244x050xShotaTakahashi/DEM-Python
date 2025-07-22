import os
import json
import subprocess
import csv
import shutil


def run_parameter_sweep():
    """dem_calc.py を複数パラメータで実行し、エネルギー保存則と跳ね返り係数の結果を集計する"""

    # スクリプトが置かれているディレクトリ (DEM/Python) を基準とする
    base_dir = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(base_dir, 'config.json')
    outputs_dir = os.path.join(base_dir, 'outputs')
    dem_script = os.path.join(base_dir, 'dem_calc.py')

    # 探索するパラメータ
    e_list = [0.999999, 0.5, 0.25]
    vx_list = [20, 50, 100]  # 初速度 
    kn_list = [1e2, 1e5, 1e8]  # バネ定数 

    # summary ファイルは outputs 配下に置くと dem_calc が削除してしまうため、
    # base_dir 直下に配置する
    energy_summary_path = os.path.join(base_dir, 'energy_conservation_summary.csv')
    restitution_summary_path = os.path.join(base_dir, 'restitution_analysis_summary.csv')

    # config.json のオリジナルを退避
    with open(config_path, 'r', encoding='utf-8') as f:
        original_config = json.load(f)

    # 出力 CSV の準備
    os.makedirs(outputs_dir, exist_ok=True)
    
    # エネルギー保存則の結果用CSV
    with open(energy_summary_path, 'w', newline='', encoding='utf-8') as energy_file:
        energy_header_writer = csv.writer(energy_file)
        energy_header_writer.writerow([
            'restitution_e',
            'initial_vx',
            'kn',
            'pre_energy',
            'post_energy',
            'energy_difference',
            'theoretical_energy',
            'relative_error_percent'
        ])
    
    # 跳ね返り係数解析結果用CSV
    with open(restitution_summary_path, 'w', newline='', encoding='utf-8') as restitution_file:
        rest_writer = csv.writer(restitution_file)
        rest_writer.writerow([
            'restitution_e',
            'initial_vx',
            'kn',
            'set_restitution',
            'calculated_restitution',
            'relative_error_percent',
            'pre_relative_velocity',
            'post_relative_velocity'
        ])

        # 各パラメータの組み合わせでループ
        for e_val in e_list:
            for vx_val in vx_list:
                for kn_val in kn_list:
                    print(f"==== 計算開始: e={e_val}, vx={vx_val}, kn={kn_val} ====")

                    # config を読み込み直して毎回クリーンな状態から編集
                    with open(config_path, 'r', encoding='utf-8') as cf:
                        config_data = json.load(cf)

                    # パラメータを上書き
                    config_data['interface']['e'] = e_val
                    # 初速度は固定粒子 0 番の vx を対象とする
                    if config_data['particles']['fixed']:
                        config_data['particles']['fixed'][0]['vx'] = vx_val
                    # kn を更新
                    config_data['interface']['particle_particle']['kn'] = kn_val

                    # config.json を保存
                    with open(config_path, 'w', encoding='utf-8') as cf:
                        json.dump(config_data, cf, indent=2, ensure_ascii=False)

                    # outputs ディレクトリをクリーンアップ
                    if os.path.exists(outputs_dir):
                        for fname in os.listdir(outputs_dir):
                            fpath = os.path.join(outputs_dir, fname)
                            if os.path.isfile(fpath):
                                os.remove(fpath)
                    else:
                        os.makedirs(outputs_dir, exist_ok=True)

                    # シミュレーションを実行
                    try:
                        subprocess.run(
                            ['python', dem_script],
                            cwd=base_dir,
                            check=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True
                        )
                    except subprocess.CalledProcessError as e:
                        print("シミュレーション実行中にエラーが発生しました:", e)
                        print("stderr:\n", e.stderr)
                        # エラー時は次へ
                        continue

                    # 衝突解析結果が書かれたファイルを解析
                    collision_csv = os.path.join(outputs_dir, 'collision_statistics.csv')
                    if not os.path.exists(collision_csv):
                        print("collision_statistics.csv が見つかりません。スキップします。")
                        continue

                    # ファイルを読み込んで各セクションを検索
                    pre_energy = post_energy = theo_energy = rel_error = None
                    set_restitution = calculated_restitution = restitution_error = None
                    pre_relative_velocity = post_relative_velocity = None
                    
                    with open(collision_csv, 'r', encoding='utf-8') as f:
                        lines = [line.strip() for line in f]

                    # エネルギー保存則セクションの解析
                    for idx, line in enumerate(lines):
                        if line.startswith('エネルギー保存則の検証'):
                            # 2行後に数値が並ぶ想定
                            if idx + 2 < len(lines):
                                vals = lines[idx + 2].split(',')
                                try:
                                    # ファイルには [pre_energy, post_energy, error_percent] の3列構成
                                    # または 4列構成で理論値が含まれている場合がある
                                    pre_energy = float(vals[0])
                                    post_energy = float(vals[1])
                                    if len(vals) >= 4:
                                        # 理論値・誤差のいずれかが追加で存在
                                        # len==4 の場合は [pre, post, theo, error] と仮定
                                        theo_energy = float(vals[2])
                                        rel_error = float(vals[3])
                                    else:
                                        # 3列の場合 val[2] は誤差
                                        rel_error = float(vals[2])
                                        # 理論値は入力 e の二乗を用いて近似計算
                                        theo_energy = pre_energy * (e_val ** 2)
                                except ValueError:
                                    pass
                            break

                    # 跳ね返り係数セクションの解析
                    for idx, line in enumerate(lines):
                        if line.startswith('跳ね返り係数の検証'):
                            if idx + 2 < len(lines):
                                vals = lines[idx + 2].split(',')
                                try:
                                    set_restitution = float(vals[0])
                                    calculated_restitution = float(vals[1])
                                    restitution_error = float(vals[2])
                                except ValueError:
                                    pass
                            break

                    # 相対速度の計算（衝突前後の速度セクションから）
                    for idx, line in enumerate(lines):
                        if line.startswith('衝突前後の速度:'):
                            if idx + 3 < len(lines):  # ヘッダー行と2つの粒子のデータ行を確認
                                try:
                                    p1_vals = lines[idx + 2].split(',')
                                    p2_vals = lines[idx + 3].split(',')
                                    pre_relative_velocity = abs(float(p1_vals[1]) - float(p2_vals[1]))
                                    post_relative_velocity = abs(float(p1_vals[2]) - float(p2_vals[2]))
                                except (ValueError, IndexError):
                                    pass
                            break

                    # エネルギー保存則の結果を書き込み
                    if None not in (pre_energy, post_energy, theo_energy, rel_error):
                        energy_diff = abs(post_energy - pre_energy)
                        with open(energy_summary_path, 'a', newline='', encoding='utf-8') as energy_file:
                            energy_writer = csv.writer(energy_file)
                            energy_writer.writerow([e_val, vx_val, kn_val, pre_energy, post_energy, 
                                                   energy_diff, theo_energy, rel_error])
                        print(f"  -> エネルギー保存則: エラー率 {rel_error:.2f}% , ΔE={energy_diff:.6f} を記録しました。")

                    # 跳ね返り係数の結果を書き込み
                    if None not in (set_restitution, calculated_restitution, restitution_error):
                        # rest_writer は restitution_summary_path のヘッダー書き込み時に生成したものを使用
                        rest_writer.writerow([e_val, vx_val, kn_val, set_restitution, 
                                               calculated_restitution, restitution_error,
                                               pre_relative_velocity, post_relative_velocity])
                        restitution_file.flush()
                        print(f"  -> 跳ね返り係数: 設定値={set_restitution:.6f}, 計算値={calculated_restitution:.6f}, "
                              f"誤差={restitution_error:.2f}% を記録しました。")

                    # ------------------- 追加: オーバーラップ妥当性グラフ作成 -------------------
                    # overlap ディレクトリに生成された比較 CSV を取得
                    overlap_dir = os.path.join(outputs_dir, 'overlap')
                    if os.path.isdir(overlap_dir):
                        csv_files = [os.path.join(overlap_dir, f) for f in os.listdir(overlap_dir) if f.startswith('overlap_velocity_comparison_') and f.endswith('.csv')]
                        if csv_files:
                            # 最新のファイルを選択（mtime が最大）
                            latest_csv = max(csv_files, key=os.path.getmtime)

                            # PNG の保存先ファイル名を作成
                            png_name = f"e={e_val}_v={vx_val}_k={kn_val}.png"
                            png_path = os.path.join(overlap_dir, png_name)

                            # check_overlap_validity.py を呼び出し (PNG 保存のみ、画面表示なし)
                            try:
                                subprocess.run(
                                    ['python', 'check_overlap_validity.py', '-f', latest_csv, '--png', png_path],
                                    cwd=base_dir,
                                    check=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    text=True
                                )
                                print(f"  -> オーバーラップ妥当性グラフを保存: {png_path}")
                            except subprocess.CalledProcessError as e:
                                print("    オーバーラップ妥当性チェック中にエラーが発生しました:", e)
                                print("    stderr:\n", e.stderr)
                    # --------------------------------------------------------------------------

    # config.json を元に戻す
    with open(config_path, 'w', encoding='utf-8') as f:
        json.dump(original_config, f, indent=2, ensure_ascii=False)

    print("====== パラメータ掃引が完了しました ======")
    print(f"エネルギー保存則の結果は {energy_summary_path} に保存されています。")
    print(f"跳ね返り係数解析の結果は {restitution_summary_path} に保存されています。")


if __name__ == '__main__':
    run_parameter_sweep() 