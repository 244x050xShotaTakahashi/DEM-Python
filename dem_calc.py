# -*- coding: utf-8 -*-

import sys
import math
import random
import time
import tkinter
import csv
import json
import os
import glob  # ファイル検索のために追加
import numpy as np  # 配列計算のために追加


class Element(object):
    def __init__(self):
        global elem_num
        self.n = elem_num #要素No.
        elem_num += 1
        self.r = 0 #半径
        self.x = 0 #X座標
        self.y = 0 #Y座標
        self.a = 0  #角度
        self.dx = 0 #X方向増加量
        self.dy = 0 #Y方向増加量
        self.da = 0 #角度増加量
        self.vx = 0 #X方向速度
        self.vy = 0 #Y方向速度
        self.va = 0 #角速度
        self.fy = 0
        self.fx = 0
        self.fm = 0
        
        self.en = [] #弾性力（直方向）
        self.es = [] #弾性力（せん断方向）

class Particle(Element):
    def __init__(self,x,y,vx=0,vy=0,r=0):
        super(Particle,self).__init__()
        self.type = 1
        
        self.x = x #X座標
        self.y = y #Y座標
        self.vx = vx #X方向速度
        self.vy = vy #Y方向速度
        
        self.r = r #半径
        self.m = 4.0/3.0*math.pi*rho*self.r**3 # 質量
        self.Ir = math.pi*rho*self.r**4/2.0 #慣性モーメント

    def config(self):
        self.en = [0 for i in range(elem_num)]
        self.es = [0 for i in range(elem_num)]
        
    def nextStep(self,dt):
        #位置更新（オイラー差分）
        ax = self.fx/self.m
        ay = self.fy/self.m
        aa = self.fm/self.Ir
        self.vx += ax*dt
        self.vy += ay*dt
        self.va += aa*dt
        self.dx = self.vx*dt
        self.dy = self.vy*dt
        self.da = self.va*dt
        self.x += self.dx
        self.y += self.dy
        self.a += self.da

class Line(Element):
    def __init__(self,x1,y1,x2,y2):
        super(Line,self).__init__()
        self.type = 2
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2

class Interface():
    def __init__(self, config):
        self.config = config
        self.e = config['e']
        # self.e = 0

    def set(self, elem1, elem2):
        if elem1.type == 1 and elem2.type == 1:
            # 粒子間
            self.kn = self.config['particle_particle']['kn']
            self.ks = self.config['particle_particle']['ks']
            self.etas = self.config['particle_particle']['etas']
            self.frc = self.config['particle_particle']['frc']
            self.etan = -2*math.log(self.e)*math.sqrt((elem1.m*self.kn)/((math.pi**2)+((math.log(self.e))**2)))
            # self.etan = self.config['particle_particle']['etan']
        elif elem1.type == 1 and elem2.type == 2:
            # 粒子-Line間
            self.kn = self.config['particle_line']['kn']
            self.ks = self.config['particle_line']['ks']
            self.etas = self.config['particle_line']['etas']
            self.frc = self.config['particle_line']['frc']
            self.etan = -2*math.log(self.e)*math.sqrt((elem1.m*self.kn)/((math.pi**2)+((math.log(self.e))**2)))
            # self.etan = self.config['particle_line']['etan']

    def get_kn_etan(self, elem1, elem2):
        self.set(elem1, elem2)
        return self.kn, self.etan

class DEM():
    def __init__(self):
        # 設定ファイルの読み込み
        with open('config.json', 'r') as f:
            self.config = json.load(f)
        
        # グローバル変数の設定
        global G, rho, elem_num
        G = self.config['simulation']['G']
        rho = self.config['simulation']['rho']
        elem_num = 0
        
        # シミュレーションパラメータの設定
        self.dt = self.config['simulation']['dt']
        # self.e = self.config['simulation']['e']
        self.e = 0
        self.max_steps = self.config['simulation']['max_steps']  
        
        self.step_times = [] # ステップごとの計算時間を記録
        self.calc_start_time = None   # 計算開始時間
        self.calc_time = 0 # 総計算時間
        self.last_time = time.time()
        self.step = 0
        self.h0 = 0
        
        # 接触追跡用の変数を追加
        self.contact_pairs = {}  # 接触中の粒子ペアを記録 {粒子ペアID: 接触開始ステップ}
        self.contact_durations = []  # 接触期間を記録 [接触開始ステップ, 接触終了ステップ, 粒子1ID, 粒子2ID]
        
        # 粒子と壁の接触追跡用の変数を追加
        self.wall_contact_pairs = {}  # 接触中の粒子-壁ペアを記録 {粒子ID_壁ID: 接触開始ステップ}
        self.wall_contact_durations = []  # 粒子-壁接触期間を記録 [接触開始ステップ, 接触終了ステップ, 粒子ID, 壁ID]
        
        # オーバーラップ追跡用の変数を追加
        self.overlap_data = {}  # 接触中の粒子ペアの最大オーバーラップを記録 {粒子ペアID: 最大オーバーラップ}
        self.max_overlap_ratios = []  # 接触終了時の最大オーバーラップ比率 [粒子1ID, 粒子2ID, 最大オーバーラップ, 比率(%)]
        
        # 粒子-壁のオーバーラップ追跡用の変数を追加
        self.wall_overlap_data = {}  # 接触中の粒子-壁ペアの最大オーバーラップを記録 {粒子ID_壁ID: 最大オーバーラップ}
        self.wall_max_overlap_ratios = []  # 接触終了時の最大オーバーラップ比率 [粒子ID, 壁ID, 最大オーバーラップ, 比率(%)]
        
        # 粒子ペアごとのオーバーラップ時系列データを保持
        self.overlap_series = {}  # {pair_id: {'time':[], 'overlap':[], 'initial_vn': float, 'p1_id': int, 'p2_id': int}}
        
        # 理論解との比較結果 (平均絶対誤差など) を保持
        self.overlap_analysis_results = []  # [[pair_id, mae], ...]
        
        # ---- 追加: 前ステップ速度保持用 ----
        self.prev_velocities = {}  # {particle_id: (vx, vy)}
        
        self.pars = []
        self.lines = []
        self.beams = []
        self.interface = Interface(self.config['interface'])
        
        # CSVファイルの設定
        self.setup_csv()
        
        # パーティクルの生成
        self.generate_particles()
        
        # 線の生成
        self.generate_lines()
        
        for p in self.pars:
            p.config()
        
        # kn と etan を取得(粒子-線間)
        kn, etan = self.interface.get_kn_etan(self.pars[0], self.lines[0])
        print("--------------------------------")
        print(f"粒子-線間のkn: {kn}, etan: {etan}")
        # はねかえり係数の計算と表示
        self.wall_restitution = self.calculate_restitution_coefficient(kn, etan)
        print(f"はねかえり係数: {self.wall_restitution}")
        
        print("--------------------------------")
        
        # 粒子数が2個以上の場合のみ粒子-粒子間の係数を計算
        if len(self.pars) >= 2:
            # kn と etan を取得(粒子-粒子間)
            kn, etan = self.interface.get_kn_etan(self.pars[0], self.pars[1])
            print(f"粒子-粒子間のkn: {kn}, etan: {etan}")
            # はねかえり係数の計算と表示
            self.particle_restitution = self.calculate_restitution_coefficient(kn, etan)
            print(f"はねかえり係数: {self.particle_restitution}")
            print("--------------------------------")
            
            # 設定された跳ね返り係数を計算値で上書き
            self.e = self.particle_restitution
        else:
            print("粒子数が1個のため、粒子-粒子間の係数計算をスキップします")
            # 粒子が1個の場合は壁との係数を使用
            self.particle_restitution = self.wall_restitution
            self.e = self.wall_restitution
            print("--------------------------------")
        
        # 衝突解析用の変数を追加
        self.collision_analysis = {
            'has_collision': False,
            'pre_collision_velocities': {},  # 衝突前の速度
            'post_collision_velocities': {},  # 衝突後の速度
            'theoretical_velocities': {},  # 理論的な速度
            'relative_errors': {},  # 相対誤差
            'energy_conservation': {}  # エネルギー保存の検証
        }
    def setup_csv(self):
        # 出力ディレクトリの作成
        os.makedirs('outputs', exist_ok=True)
        
        # outputsディレクトリ内の既存ファイルをすべて削除
        for filename in os.listdir('outputs'):
            file_path = os.path.join('outputs', filename)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                    # print(f"削除: {file_path}")
            except Exception as e:
                print(f"ファイル削除エラー {file_path}: {e}")

    def save_state(self):
        """現在の粒子の状態をCSVファイルに保存"""
        filename = f'outputs/dem_results.csv.{self.step}'
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            # ヘッダーの書き込み
            writer.writerow(['Step', 'Particle', 'X', 'Y', 'Angle', 'Vx', 'Vy', 'Va', 'Fx', 'Fy', 'Fm'])
            
            # 各粒子のデータを書き込み
            for p in self.pars:
                writer.writerow([
                    self.step,      # ステップ数
                    p.n,           # 粒子番号
                    p.x,           # X座標
                    p.y,           # Y座標
                    p.a,           # 角度
                    p.vx,          # X方向速度
                    p.vy,          # Y方向速度
                    p.va,          # 角速度
                    p.fx,          # X方向力
                    p.fy,          # Y方向力
                    p.fm           # モーメント
                ])
            
    def generate_particles(self):
        particle_config = self.config['particles']
        if particle_config['random']:
            self.generate_random_particles(particle_config['count'], particle_config['x_range'], particle_config['y_range'])
        else:
            for p in particle_config['fixed']:
                self.pars.append(Particle(p['x'], p['y'], p['vx'], p['vy'], p['r']))

    def generate_random_particles(self, count, x_range, y_range):
        max_attempts = 1000
        particles_placed = 0
        while particles_placed < count:
            for attempt in range(max_attempts):
                x = random.uniform(x_range[0], x_range[1])
                y = random.uniform(y_range[0], y_range[1])
                r = self.config['particles']['radius']
                p = Particle(x, y, 0, 0, r)
                
                if not self._hitParticle(x, y, p.r, self.pars) and not self._hitLine(x, y, p.r, self.lines):
                    self.pars.append(p)
                    particles_placed += 1
                    break
                else:
                    print(f"警告: 粒子 {particles_placed} の配置に失敗しました。再試行します。")
                    if attempt > max_attempts:
                        break
        
        print(f"配置された粒子の総数: {particles_placed}")

    def generate_lines(self):
        for l in self.config['lines']:
            self.lines.append(Line(l['x1'], l['y1'], l['x2'], l['y2']))

    def calculate_restitution_coefficient(self, kn, etan):
        r = self.config['particles']['radius']
        restitution_coefficient = math.exp(-math.pi / (2 * math.sqrt(((4.0 / 3.0 * math.pi * rho * r ** 3) * kn / (etan) ** 2) - (1 / 4))))
        return restitution_coefficient  # 計算した跳ね返り係数を返すように変更
            
    def _hitParticle(self,x,y,r,pars):
        hit = False
        for p in pars:
            lx = p.x - x
            ly = p.y - y
            ld = (lx**2+ly**2)**0.5
            if (p.r+r)>=ld:
                hit = True
                break
        return hit
    
    def _hitLine(self,px,py,pr,lines):
        hit = False
        for l in lines:
            th0 = math.atan2(l.y2-l.y1,l.x2-l.x1)
            th1 = math.atan2(py-l.y1,px-l.x1)
            a = math.sqrt((px-l.x1)**2+(py-l.y1)**2)
            d = abs(a*math.sin(th1-th0))
            if d < pr:
                b = math.sqrt((px-l.x2)**2+(py-l.y2)**2)
                s = math.sqrt((l.x2-l.x1)**2+(l.y2-l.y1)**2)
                if a < s and b < s:
                    hit = True
                elif a < b and a < pr:
                    hit = True
                elif b < pr:
                    hit = True
                if hit:
                    break
        return hit
    
    def calcForce(self):
        # === 追加: 衝突判定に入る前に前ステップの速度を保存 ===
        self.prev_velocities = {p.n: (p.vx, p.vy) for p in self.pars}
        
        # for idx, l in enumerate(self.lines):
        #     print(f"Python: 線要素インデックス: {idx}, 内部ID: {l.n}")  
    
        #2粒子間の接触判定
        
        #近傍粒子探索をする場合はここに追加実装
        
        # 現在のステップで接触している粒子ペアを記録
        current_contacts = set()
        
        combs = [(p1, p2) for p1 in self.pars for p2 in self.pars]
        for p1,p2 in combs:
            if p1.n == p2.n:
                continue
            
            lx = p1.x - p2.x
            ly = p1.y - p2.y
            ld = math.sqrt(lx**2+ly**2)
            
            if ld == 0:
                continue  # 距離がゼロの場合はスキップ
            
            # 粒子ペアのIDを作成（小さい番号を先に）
            pair_id = f"{min(p1.n, p2.n)}_{max(p1.n, p2.n)}"
            
            if (p1.r+p2.r)>ld:
                # オーバーラップを計算
                overlap = (p1.r + p2.r) - ld
                
                # --------- オーバーラップ時系列の毎ステップ記録 ---------
                cos_tmp = lx / ld
                sin_tmp = ly / ld
                # 法線方向の相対速度を計算
                vn_current = (p1.vx - p2.vx) * cos_tmp + (p1.vy - p2.vy) * sin_tmp
                
                if pair_id not in self.overlap_series:
                    # ---- 追加: 衝突直前(前ステップ)の情報を初期値として格納 ----
                    v1_prev = self.prev_velocities.get(p1.n, (p1.vx, p1.vy))
                    v2_prev = self.prev_velocities.get(p2.n, (p2.vx, p2.vy))
                    vn_prev = (v1_prev[0] - v2_prev[0]) * cos_tmp + (v1_prev[1] - v2_prev[1]) * sin_tmp
                    t_prev = (self.step - 1) * self.dt  # 衝突開始直前の時刻
                    self.overlap_series[pair_id] = {
                        'time': [t_prev],
                        'overlap': [0.0],
                        'velocity': [vn_prev],
                        'initial_vn': vn_prev,
                        'p1_id': p1.n,
                        'p2_id': p2.n
                    }
                # ここからは従来どおり現在ステップのデータを追加
                self.overlap_series[pair_id]['time'].append(self.step * self.dt)
                self.overlap_series[pair_id]['overlap'].append(overlap)
                self.overlap_series[pair_id]['velocity'].append(vn_current)
                # --------- ここまで追加 ---------
                
                cos_a = lx/ld
                sin_a = ly/ld
                self.force2par(p1,p2,cos_a,sin_a)
                
                # 現在接触しているペアを記録
                current_contacts.add(pair_id)
                
                # 新しく接触が始まった場合
                if pair_id not in self.contact_pairs:
                    self.contact_pairs[pair_id] = self.step
                    # 衝突開始時にCSVファイルを保存
                    self.save_state()
                    print(f"\n粒子{p1.n}と粒子{p2.n}の衝突開始: ステップ{self.step}")
            else:
                p1.en[p2.n] = 0.0
                p1.es[p2.n] = 0.0
                
                # 接触が終了した場合
                if pair_id in self.contact_pairs:
                    start_step = self.contact_pairs[pair_id]
                    self.contact_durations.append([start_step, self.step, p1.n, p2.n])
                    
                    # 最大オーバーラップと比率を記録
                    if pair_id in self.overlap_data:
                        max_overlap = self.overlap_data[pair_id]
                        # 粒子の半径の小さい方を基準にする
                        reference_radius = min(p1.r, p2.r)
                        overlap_ratio = (max_overlap / reference_radius) * 100.0  # パーセント表示
                        self.max_overlap_ratios.append([p1.n, p2.n, max_overlap, overlap_ratio])
                        print(f"粒子{p1.n}と粒子{p2.n}の最大オーバーラップ: {max_overlap:.6f}, 半径比: {overlap_ratio:.2f}%")
                        del self.overlap_data[pair_id]
                    
                    del self.contact_pairs[pair_id]
                    # オーバーラップ時系列をCSVへ保存
                    if pair_id in self.overlap_series:
                        self.save_overlap_series(pair_id)
                        del self.overlap_series[pair_id]
                    # 衝突終了時にCSVファイルを保存
                    self.save_state()
                    print(f"\n粒子{p1.n}と粒子{p2.n}の衝突終了: ステップ{self.step}")
        
        #粒子と線の接触判定
        # 現在のステップで接触している粒子-壁ペアを記録
        current_wall_contacts = set()
                
        combs = [(l, p) for l in self.lines for p in self.pars]
        for l,p in combs:
            hit = False
            th0 = math.atan2(l.y2-l.y1,l.x2-l.x1)
            th1 = math.atan2(p.y-l.y1,p.x-l.x1)
            a = math.sqrt((p.x-l.x1)**2+(p.y-l.y1)**2)
            d = abs(a*math.sin(th1-th0))
            if d < p.r:
                b = math.sqrt((p.x-l.x2)**2+(p.y-l.y2)**2)
                s = math.sqrt((l.x2-l.x1)**2+(l.y2-l.y1)**2)
                if a < s and b < s:
                    s1 = math.sqrt(a**2-d**2)
                    x = l.x1 + s1*math.cos(th0)
                    y = l.y1 + s1*math.sin(th0)
                    hit = True
                elif a < b and a < p.r:
                    x = l.x1
                    y = l.y1
                    hit = True
                elif b < p.r:
                    x = l.x2
                    y = l.y2
                    hit = True       
            if hit:
                # 壁とのオーバーラップを計算
                lx = x - p.x
                ly = y - p.y
                ld = math.sqrt(lx**2+ly**2)
                overlap = p.r - ld
                
                # 壁とのオーバーラップを記録
                wall_pair_id = f"{p.n}_{l.n}"
                if wall_pair_id in self.wall_overlap_data:
                    self.wall_overlap_data[wall_pair_id] = max(self.wall_overlap_data[wall_pair_id], overlap)
                else:
                    self.wall_overlap_data[wall_pair_id] = overlap
                
                cos_a = lx/ld
                sin_a = ly/ld
                self.force2line(p,l,cos_a,sin_a)
                # print(f"Python: 接触あり - 粒子: {p.n}, 線: {self.lines.index(l)}, 線内部ID: {l.n}")
                # 粒子-壁ペアのIDを作成
                wall_pair_id = f"{p.n}_{l.n}"
                
                current_wall_contacts.add(wall_pair_id)
                
                # 新しく接触が始まった場合
                if wall_pair_id not in self.wall_contact_pairs:
                    self.wall_contact_pairs[wall_pair_id] = self.step
                    
            else:
                p.en[l.n] = 0.0
                p.es[l.n] = 0.0
                
                # 接触が終了した場合
                wall_pair_id = f"{p.n}_{l.n}"
                if wall_pair_id in self.wall_contact_pairs:
                    start_step = self.wall_contact_pairs[wall_pair_id]
                    self.wall_contact_durations.append([start_step, self.step, p.n, l.n])
                    
                    # 最大オーバーラップと比率を記録
                    if wall_pair_id in self.wall_overlap_data:
                        max_overlap = self.wall_overlap_data[wall_pair_id]
                        overlap_ratio = (max_overlap / p.r) * 100.0  # パーセント表示
                        self.wall_max_overlap_ratios.append([p.n, l.n, max_overlap, overlap_ratio])
                        print(f"粒子{p.n}と壁{l.n}の最大オーバーラップ: {max_overlap:.6f}, 半径比: {overlap_ratio:.2f}%")
                        del self.wall_overlap_data[wall_pair_id]
                    
                    del self.wall_contact_pairs[wall_pair_id]
        
        #外力
        for p in self.pars:
            p.fy += -G*p.m #重力
        for b in self.beams:
            b.fy += -G*b.m
            
    def force2par(self,p1,p2,cos_a,sin_a):
        
        #相対的変位増分
        un = +(p1.dx-p2.dx)*cos_a+(p1.dy-p2.dy)*sin_a
        us = -(p1.dx-p2.dx)*sin_a+(p1.dy-p2.dy)*cos_a+(p1.r*p1.da+p2.r*p2.da)
        #相対的速度増分
        vn = +(p1.vx-p2.vx)*cos_a+(p1.vy-p2.vy)*sin_a
        vs = -(p1.vx-p2.vx)*sin_a+(p1.vy-p2.vy)*cos_a+(p1.r*p1.va+p2.r*p2.va)
        
        self.interface.set(p1,p2)
        
        #合力（局所座標系）
        p1.en[p2.n] += self.interface.kn*un
        p1.es[p2.n] += self.interface.ks*us
        hn = p1.en[p2.n] + self.interface.etan*vn
        hs = p1.es[p2.n] + self.interface.etas*vs
        
        if hn <= 0.0: 
            #法線力がなければ、せん断力は０             
            hs = 0.0         
        elif abs(hs) >= self.interface.frc*hn:
            #摩擦力以上のせん断力は働かない
            if abs(hs) > 1e-12:  # 0 で割ることを回避
                hs = self.interface.frc*abs(hn)*hs/abs(hs)
            else:
                hs = 0.0        
        
        #全体合力（全体座標系）
        p1.fx += -hn*cos_a + hs*sin_a
        p1.fy += -hn*sin_a - hs*cos_a
        p1.fm -=  p1.r*hs
        
    def force2line(self,p1,line,cos_a,sin_a):
        
        #相対的変位増分
        un = +(p1.dx)*cos_a+(p1.dy)*sin_a
        us = -(p1.dx)*sin_a+(p1.dy)*cos_a+(p1.r*p1.da)
        #相対的速度増分
        vn = +(p1.vx)*cos_a+(p1.vy)*sin_a
        vs = -(p1.vx)*sin_a+(p1.vy)*cos_a+(p1.r*p1.va)
        
        # interfaceパラメータの設定を追加
        self.interface.set(p1,line)
        
        #合力
        p1.en[line.n] += self.interface.kn*un
        p1.es[line.n] += self.interface.ks*us
        
        # 摩擦力の計算
        hn = p1.en[line.n] + self.interface.etan*vn
        # せん断力の計算
        hs = p1.es[line.n] + self.interface.etas*vs
        
        if hn <= 0.0:
            #法線力がなければ、せん断力は０
            hs = 0.0
        elif abs(hs)-self.interface.frc*hn >= 0.0:
            #摩擦力以上のせん断力は働かない
            if abs(hs) > 1e-12:
                hs = self.interface.frc*abs(hn)*hs/abs(hs)
            else:
                hs = 0.0        
        
        #合力
        p1.fx += -hn*cos_a + hs*sin_a
        p1.fy += -hn*sin_a - hs*cos_a
        p1.fm -=  p1.r*hs
        
    def debug_output(self, tag):
        """現在のstepとタグ、各粒子の状態（位置、速度、力など）を出力する"""
        filename = "./debug/debug_log.txt"
        with open(filename, mode="a", encoding="utf-8") as f:
            f.write(f"Step {self.step} {tag}\n")
            for p in self.pars:
                # 各粒子の状態を出力（必要な項目を追加）
                f.write(f"Particle {p.n}: x={p.x:.7f}, y={p.y:.7f}, "
                        f"vx={p.vx:.7f}, vy={p.vy:.7f}, "
                        f"Fx={p.fx:.5f}, Fy={p.fy:.5f}, Fm={p.fm:.5f}\n")
            f.write("\n")
            
    def calculate_theoretical_velocity(self, particle, time):
        # 粒子の情報を使用してfrc_thresholdを計算
        g = G  # 重力加速度
        flag = 0 # if文の判別用
        theta = math.atan2(self.lines[0].y2 - self.lines[0].y1, self.lines[0].x2 - self.lines[0].x1)
        bottom = (particle.m * particle.r**2)+ particle.Ir
        top = particle.m * g * particle.r**2*math.sin(theta)
        # 慣性モーメントを使用してfrc_thresholdを計算
        frc_threshold = (particle.Ir / bottom)*math.tan(theta)
    
        if self.interface.frc > frc_threshold:
            # 摩擦が大きい場合(転がる)場合の理論解(flag==0)
            theoretical_velocity = (top/bottom)*time
            return theoretical_velocity ,flag
        else:
            # 摩擦が小さい場合(滑る)場合の理論解(flag==1)
            theoretical_velocity = g*time*(math.sin(theta) - self.interface.frc*math.cos(theta))
            flag += 1
            return theoretical_velocity ,flag
        
    def calcStep(self):
        # 出力間隔を1000ステップに設定
        self.output_interval = self.config['simulation']['output_interval']
    
        if self.step == 0:
            print("--------------------------------")
            print("計算開始")
            # 初期状態を保存
            self.save_state()
            # 計算開始時間を記録
            self.calc_start_time = time.time()

        # 計算時間の計測開始
        step_start_time = time.time()
            
        # 既存の計算処理
        for p in self.pars:
            p.fx = 0
            p.fy = 0
            p.fm = 0
        self.calcForce()
        
        # 更新処理（位置・速度更新等）
        for p in self.pars:
            p.nextStep(self.dt)
        self.step += 1
        
        # ステップ計算時間を計測
        step_end_time = time.time()
        step_time = step_end_time - step_start_time
        self.step_times.append(step_time)
        self.calc_time += step_time
            
        # 指定間隔でCSVファイルに出力
        if self.step % self.output_interval == 0:
            self.save_state()
            print(f"\rStep: {self.step}/{self.max_steps} ({(self.step/self.max_steps)*100:.1f}%)", end="")

        # 計算終了時の処理
        if self.step >= self.max_steps:
            # 総計算時間を最終計測
            total_elapsed = time.time() - self.calc_start_time
            print(f"\n計算完了")
            print("--------------------------------")
            print(f"総ステップ数: {self.max_steps}")
            print(f"総計算時間: {total_elapsed:.2f} 秒")
            print(f"1ステップあたりの平均時間: {(total_elapsed/self.max_steps)*1000:.3f} ミリ秒")
            print("--------------------------------")
            
            # 接触期間の統計情報を出力
            print("--------------------------------")
            print("接触期間の統計情報:")
            
            # 粒子-粒子間の接触統計
            if self.contact_durations:
                total_contacts = len(self.contact_durations)
                total_steps = sum(duration[1] - duration[0] for duration in self.contact_durations)
                avg_duration = total_steps / total_contacts if total_contacts > 0 else 0
                
                print("粒子-粒子間の接触:")
                print(f"  総接触回数: {total_contacts}")
                print(f"  平均接触ステップ数: {avg_duration:.2f}")
                
                # オーバーラップの統計情報を追加
                if self.max_overlap_ratios:
                    overlap_ratios = [data[3] for data in self.max_overlap_ratios]
                    avg_overlap = sum(overlap_ratios) / len(overlap_ratios)
                    max_overlap = max(overlap_ratios)
                    print(f"  平均最大オーバーラップ比: {avg_overlap:.2f}%")
                    print(f"  最大オーバーラップ比: {max_overlap:.2f}%")
            else:
                print("粒子-粒子間の接触は発生しませんでした")
            
            # 粒子-壁間の接触統計
            if self.wall_contact_durations:
                total_wall_contacts = len(self.wall_contact_durations)
                total_wall_steps = sum(duration[1] - duration[0] for duration in self.wall_contact_durations)
                avg_wall_duration = total_wall_steps / total_wall_contacts if total_wall_contacts > 0 else 0
                
                print("粒子-壁間の接触:")
                print(f"  総接触回数: {total_wall_contacts}")
                print(f"  平均接触ステップ数: {avg_wall_duration:.2f}")
                
                # オーバーラップの統計情報を追加
                if self.wall_max_overlap_ratios:
                    wall_overlap_ratios = [data[3] for data in self.wall_max_overlap_ratios]
                    avg_wall_overlap = sum(wall_overlap_ratios) / len(wall_overlap_ratios)
                    max_wall_overlap = max(wall_overlap_ratios)
                    print(f"  平均最大オーバーラップ比: {avg_wall_overlap:.2f}%")
                    print(f"  最大オーバーラップ比: {max_wall_overlap:.2f}%")
            else:
                print("粒子-壁間の接触は発生しませんでした")
            
            # 接触期間の詳細をCSVファイルに出力
            self.save_contact_durations()
            print("--------------------------------")
            
            # アニメーション作成
            try:
                from visualizer import create_dem_animation
                print("アニメーション作成中...")
                create_dem_animation()
            except ImportError:
                print("visualizer.pyが見つかりません。アニメーション作成をスキップします。")
            
            sys.exit()
    
    def save_contact_durations(self):
        """接触期間の詳細をCSVファイルに保存"""
        # 粒子-粒子間の接触データを保存
        filename = 'outputs/collision_statistics.csv'
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            # ヘッダーの書き込み
            writer.writerow(['粒子-粒子間の接触'])
            writer.writerow(['接触開始ステップ', '接触終了ステップ', '接触ステップ数', '粒子1ID', '粒子2ID', '最大オーバーラップ', 'オーバーラップ/半径比(%)'])
            
            # 各接触期間のデータを書き込み
            for i, duration in enumerate(self.contact_durations):
                start_step, end_step, p1_id, p2_id = duration
                contact_steps = end_step - start_step
                
                # オーバーラップデータの取得（ある場合）
                max_overlap = 0.0
                overlap_ratio = 0.0
                for overlap_data in self.max_overlap_ratios:
                    if (overlap_data[0] == p1_id and overlap_data[1] == p2_id) or (overlap_data[0] == p2_id and overlap_data[1] == p1_id):
                        max_overlap = overlap_data[2]
                        overlap_ratio = overlap_data[3]
                        break
                
                writer.writerow([start_step, end_step, contact_steps, p1_id, p2_id, max_overlap, overlap_ratio])
            
            writer.writerow(['--------------------------------'])
            
            writer.writerow(['粒子-壁間の接触'])
            writer.writerow(['接触開始ステップ', '接触終了ステップ', '接触ステップ数', '粒子ID', '壁ID', '最大オーバーラップ', 'オーバーラップ/半径比(%)'])
            for i, duration in enumerate(self.wall_contact_durations):
                start_step, end_step, p_id, l_id = duration
                contact_steps = end_step - start_step
                
                # オーバーラップデータの取得（ある場合）
                max_overlap = 0.0
                overlap_ratio = 0.0
                for overlap_data in self.wall_max_overlap_ratios:
                    if overlap_data[0] == p_id and overlap_data[1] == l_id:
                        max_overlap = overlap_data[2]
                        overlap_ratio = overlap_data[3]
                        break
                        
                writer.writerow([start_step, end_step, contact_steps, p_id, l_id, max_overlap, overlap_ratio])
        
        # オーバーラップの統計情報を出力
        self.save_overlap_statistics()
        
        # 衝突解析を実行
        self.analyze_collision()
    
    def save_overlap_statistics(self):
        """オーバーラップの統計情報をCSVファイルに保存"""
        overlap_dir = os.path.join('outputs', 'overlap')
        os.makedirs(overlap_dir, exist_ok=True)
        filename = os.path.join(overlap_dir, 'overlap_statistics.csv')
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # 粒子-粒子間のオーバーラップ統計
            writer.writerow(['粒子-粒子間の最大オーバーラップ統計'])
            writer.writerow(['粒子1ID', '粒子2ID', '最大オーバーラップ', 'オーバーラップ/半径比(%)'])
            
            for data in self.max_overlap_ratios:
                writer.writerow(data)
                
            writer.writerow(['--------------------------------'])
            
            # 粒子-壁間のオーバーラップ統計
            writer.writerow(['粒子-壁間の最大オーバーラップ統計'])
            writer.writerow(['粒子ID', '壁ID', '最大オーバーラップ', 'オーバーラップ/半径比(%)'])
            
            for data in self.wall_max_overlap_ratios:
                writer.writerow(data)
                
            # オーバーラップの要約統計量
            writer.writerow(['--------------------------------'])
            writer.writerow(['オーバーラップ要約統計'])
            
            # 粒子-粒子間
            particle_overlaps = [data[3] for data in self.max_overlap_ratios]
            if particle_overlaps:
                avg_particle_overlap = sum(particle_overlaps) / len(particle_overlaps)
                max_particle_overlap = max(particle_overlaps) if particle_overlaps else 0
                min_particle_overlap = min(particle_overlaps) if particle_overlaps else 0
                writer.writerow(['粒子-粒子間の平均オーバーラップ比(%)', avg_particle_overlap])
                writer.writerow(['粒子-粒子間の最大オーバーラップ比(%)', max_particle_overlap])
                writer.writerow(['粒子-粒子間の最小オーバーラップ比(%)', min_particle_overlap])
            else:
                writer.writerow(['粒子-粒子間のオーバーラップデータなし'])
            
            # 粒子-壁間
            wall_overlaps = [data[3] for data in self.wall_max_overlap_ratios]
            if wall_overlaps:
                avg_wall_overlap = sum(wall_overlaps) / len(wall_overlaps)
                max_wall_overlap = max(wall_overlaps) if wall_overlaps else 0
                min_wall_overlap = min(wall_overlaps) if wall_overlaps else 0
                writer.writerow(['粒子-壁間の平均オーバーラップ比(%)', avg_wall_overlap])
                writer.writerow(['粒子-壁間の最大オーバーラップ比(%)', max_wall_overlap])
                writer.writerow(['粒子-壁間の最小オーバーラップ比(%)', min_wall_overlap])
            else:
                writer.writerow(['粒子-壁間のオーバーラップデータなし'])

    def analyze_collision(self):
        """衝突前後の速度を分析し、理論値との比較を行う"""
        if len(self.contact_durations) == 0:
            print("分析する衝突データがありません。")
            return
        
        # 最初の衝突を取得
        first_collision = self.contact_durations[0]
        start_step, end_step, p1_id, p2_id = first_collision
        
        print(f"衝突分析: 粒子{p1_id}と粒子{p2_id}の衝突 (ステップ {start_step}-{end_step})")
        
        # 衝突前のデータは衝突開始ステップのデータを使用
        pre_collision_file = f'outputs/dem_results.csv.0'
        
        # 衝突後のデータは衝突終了ステップのデータを使用
        post_collision_file = f'outputs/dem_results.csv.{end_step}'
        
        # ファイルの存在確認
        if not os.path.exists(pre_collision_file):
            print(f"衝突前のデータファイルが見つかりません: {pre_collision_file}")
            return
        
        if not os.path.exists(post_collision_file):
            print(f"衝突後のデータファイルが見つかりません: {post_collision_file}")
            return
        
        print(f"衝突前: {pre_collision_file}, 衝突後: {post_collision_file}")
        
        # 衝突前後の速度を取得
        pre_velocities = self.get_particle_velocities(pre_collision_file, [p1_id, p2_id])
        post_velocities = self.get_particle_velocities(post_collision_file, [p1_id, p2_id])
        
        # 理論解を計算
        theoretical_velocities = self.calculate_theoretical_post_collision_velocities(
            pre_velocities[p1_id], pre_velocities[p2_id], 
            self.pars[p1_id].m, self.pars[p2_id].m, 
            self.particle_restitution, p1_id, p2_id)
        
        # 相対誤差を計算
        relative_errors = {}
        for pid, theo_vel in theoretical_velocities.items():
            v_theo = theo_vel[0]  # x方向速度のみを考慮
            v_calc = post_velocities[pid][0]  # x方向速度のみを考慮
            if abs(v_theo) > 1e-10:  # ゼロ除算を避ける
                rel_error = abs((v_calc - v_theo) / v_theo) * 100  # パーセント表示
            else:
                rel_error = 0.0 if abs(v_calc) < 1e-10 else float('inf')
            relative_errors[pid] = rel_error
            # 個別に結果を表示
            print(f"粒子{pid}: 衝突前速度={pre_velocities[pid][0]:.2f}, 衝突後速度={v_calc:.2f}, 理論値={v_theo:.2f}, 相対誤差={rel_error:.2f}%")
        
        # --- エネルギー保存則の検証（一次元理論式に基づく） ---
        pre_energy = self.calculate_kinetic_energy(pre_velocities, [p1_id, p2_id])
        post_energy = self.calculate_kinetic_energy(post_velocities, [p1_id, p2_id])
        theoretical_post_energy = self.calculate_kinetic_energy(theoretical_velocities, [p1_id, p2_id])
        
        if theoretical_post_energy > 1e-10:
            energy_error = abs((post_energy - theoretical_post_energy) / theoretical_post_energy) * 100
        else:
            energy_error = 0.0
        
        print(f"エネルギー保存: 衝突前={pre_energy:.2f}, 衝突後={post_energy:.2f}, 理論値={theoretical_post_energy:.2f}, 誤差={energy_error:.2f}%")
        
        # --- 跳ね返り係数の逆算と誤差計算 ---
        # 相対速度の計算
        pre_relative_velocity = pre_velocities[p1_id][0] - pre_velocities[p2_id][0]
        post_relative_velocity = post_velocities[p1_id][0] - post_velocities[p2_id][0]
        
        # 跳ね返り係数の逆算（相対速度の比）
        if abs(pre_relative_velocity) > 1e-10:  # ゼロ除算を避ける
            calculated_restitution = abs(post_relative_velocity / pre_relative_velocity)
        else:
            calculated_restitution = 0.0
        
        # 設定値との誤差計算
        restitution_error = abs((calculated_restitution - self.particle_restitution) / self.particle_restitution) * 100
        
        print(f"跳ね返り係数: 設定値={self.particle_restitution:.6f}, 計算値={calculated_restitution:.6f}, 誤差={restitution_error:.2f}%")
        
        # 結果を保存
        self.collision_analysis['has_collision'] = True
        self.collision_analysis['pre_collision_velocities'] = pre_velocities
        self.collision_analysis['post_collision_velocities'] = post_velocities
        self.collision_analysis['theoretical_velocities'] = theoretical_velocities
        self.collision_analysis['relative_errors'] = relative_errors
        self.collision_analysis['energy_conservation'] = {
            'pre_energy': pre_energy,
            'post_energy': post_energy,
            'theoretical_post_energy': theoretical_post_energy,
            'error_percent': energy_error
        }
        self.collision_analysis['restitution'] = {
            'set_value': self.particle_restitution,
            'calculated_value': calculated_restitution,
            'error_percent': restitution_error
        }
        
        # 解析結果を出力
        self.save_collision_analysis()
    
    def get_particle_velocities(self, file_path, particle_ids):
        """指定されたCSVファイルから粒子の速度を取得"""
        velocities = {}
        with open(file_path, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # ヘッダーをスキップ
            for row in reader:
                pid = int(row[1])
                if pid in particle_ids:
                    vx = float(row[5])
                    vy = float(row[6])
                    velocities[pid] = (vx, vy)
        return velocities
    
    def calculate_theoretical_post_collision_velocities(self, v1, v2, m1, m2, e, p1_id, p2_id):
        """一次元衝突後の理論的速度を計算"""
        # 一次元衝突のみを考慮（x方向）
        v1x, _ = v1
        v2x, _ = v2
        
        # 衝突前の速度を初速度として扱う
        v1x_initial = v1x
        v2x_initial = v2x
        
        # 完全弾性衝突の場合の理論解
        v1x_new = ((m1 - e*m2)*v1x_initial + (1+e)*m2*v2x_initial) / (m1 + m2)
        v2x_new = ((1+e)*m1*v1x_initial + (m2 - e*m1)*v2x_initial) / (m1 + m2)
        
        return {
            p1_id: (v1x_new, 0),
            p2_id: (v2x_new, 0)
        }
    
    def calculate_kinetic_energy(self, velocities, particle_ids):
        """粒子の運動エネルギーを計算"""
        energy = 0.0
        for pid in particle_ids:
            vx, vy = velocities[pid]
            m = self.pars[pid].m
            energy += 0.5 * m * (vx**2 + vy**2)
        return energy
    
    def save_collision_analysis(self):
        """衝突解析結果をCSVファイルに保存"""
        if not self.collision_analysis['has_collision']:
            return
            
        # 既存のcontact_durationsファイルを読み込む
        existing_lines = []
        if os.path.exists('outputs/collision_statistics.csv'):
            with open('outputs/collision_statistics.csv', 'r') as f:
                existing_lines = f.readlines()
        
        # 衝突解析情報を追加
        with open('outputs/collision_statistics.csv', 'w') as f:
            # 既存の内容を書き込む
            for line in existing_lines:
                f.write(line)
            
            # 区切り線
            f.write('\n--------------------------------\n')
            f.write('衝突解析結果\n')
            
            # 衝突前後の速度
            f.write('\n衝突前後の速度:\n')
            f.write('粒子ID,衝突前Vx,衝突後Vx,理論解Vx,相対誤差(%)\n')
            
            for pid in self.collision_analysis['pre_collision_velocities'].keys():
                pre_vx = self.collision_analysis['pre_collision_velocities'][pid][0]
                post_vx = self.collision_analysis['post_collision_velocities'][pid][0]
                theo_vx = self.collision_analysis['theoretical_velocities'][pid][0]
                rel_error = self.collision_analysis['relative_errors'][pid]
                
                f.write(f'{pid},{pre_vx:.6f},{post_vx:.6f},{theo_vx:.6f},{rel_error:.6f}\n')
            
            # エネルギー保存
            f.write('\nエネルギー保存則の検証:\n')
            energy_data = self.collision_analysis['energy_conservation']
            f.write(f'衝突前エネルギー,理論衝突後エネルギー,計算衝突後エネルギー,相対誤差(%)\n')
            f.write(f'{energy_data["pre_energy"]:.6f},{energy_data["theoretical_post_energy"]:.6f},{energy_data["post_energy"]:.6f},{energy_data["error_percent"]:.6f}\n')
            
            # 跳ね返り係数の検証
            f.write('\n跳ね返り係数の検証:\n')
            restitution_data = self.collision_analysis['restitution']
            f.write(f'設定値,計算値,相対誤差(%)\n')
            f.write(f'{restitution_data["set_value"]:.6f},{restitution_data["calculated_value"]:.6f},{restitution_data["error_percent"]:.6f}\n')

    def compute_theoretical_overlap(self, t_list, x0, v0, k, c, m_eff):
        """減衰自由振動 (一次元) の理論式を用いてオーバーラップ量を計算"""
        result = []
        if m_eff <= 0 or k <= 0:
            return [0.0 for _ in t_list]
        omega_n = math.sqrt(k / m_eff)
        zeta = c / (2.0 * math.sqrt(k * m_eff))
        for t in t_list:
            if zeta < 1.0:
                # 減衰自由振動(zeta<1.0)
                omega_d = omega_n * math.sqrt(1 - zeta ** 2) # 減衰固有角振動数
                exp_term = math.exp(-zeta * omega_n * t) # 減衰項
                overlap_t = exp_term * (
                    x0 * math.cos(omega_d * t) +
                    (v0 + zeta * omega_n * x0) / omega_d * math.sin(omega_d * t)
                ) # 減衰自由振動の解
            elif abs(zeta - 1.0) < 1e-8:
                # 臨界減衰(zeta=1.0)
                exp_term = math.exp(-omega_n * t) # 減衰項
                overlap_t = exp_term * (x0 + (v0 + omega_n * x0) * t) # 臨界減衰の解
            else:
                # 過減衰(zeta>1.0)
                r1 = -omega_n * (zeta - math.sqrt(zeta ** 2 - 1)) # 過減衰の解の係数
                r2 = -omega_n * (zeta + math.sqrt(zeta ** 2 - 1)) # 過減衰の解の係数
                A = (v0 - r2 * x0) / (r1 - r2) # 過減衰の解の係数
                B = x0 - A # 過減衰の解の係数
                overlap_t = A * math.exp(r1 * t) + B * math.exp(r2 * t) # 過減衰の解
            result.append(overlap_t)
        return result

    def compute_theoretical_velocity(self, t_list, x0, v0, k, c, m_eff):
        """減衰自由振動 (一次元) の理論式を用いて速度を計算"""
        result = []
        if m_eff <= 0 or k <= 0:
            return [0.0 for _ in t_list]
        omega_n = math.sqrt(k / m_eff)
        zeta = c / (2.0 * math.sqrt(k * m_eff))
        for t in t_list:
            if zeta < 1.0:
                # 減衰自由振動(zeta<1.0)の速度解
                omega_d = omega_n * math.sqrt(1 - zeta ** 2)
                exp_term = math.exp(-zeta * omega_n * t)
                # v(t) = exp(-ζωₙt) × [(-ζωₙx₀ - (v₀+ζωₙx₀))cos(ω_d t) + (ω_d x₀ - ζωₙ(v₀+ζωₙx₀)/ω_d)sin(ω_d t)]
                velocity_t = exp_term * (
                    (-zeta * omega_n * x0 - (v0 + zeta * omega_n * x0)) * math.cos(omega_d * t) +
                    (omega_d * x0 - zeta * omega_n * (v0 + zeta * omega_n * x0) / omega_d) * math.sin(omega_d * t)
                )
            elif abs(zeta - 1.0) < 1e-8:
                # 臨界減衰(zeta=1.0)の速度解
                exp_term = math.exp(-omega_n * t)
                # v(t) = e^{-ωₙt} × [-ωₙ(x₀ + (v₀+ωₙx₀)t) + (v₀+ωₙx₀)]
                velocity_t = exp_term * (
                    -omega_n * (x0 + (v0 + omega_n * x0) * t) + (v0 + omega_n * x0)
                )
            else:
                # 過減衰(zeta>1.0)の速度解
                r1 = -omega_n * (zeta - math.sqrt(zeta ** 2 - 1))
                r2 = -omega_n * (zeta + math.sqrt(zeta ** 2 - 1))
                A = (v0 - r2 * x0) / (r1 - r2)
                B = x0 - A
                # v(t) = A r₁ e^{r₁t} + B r₂ e^{r₂t}
                velocity_t = A * r1 * math.exp(r1 * t) + B * r2 * math.exp(r2 * t)
            result.append(velocity_t)
        return result

    def save_overlap_series(self, pair_id):
        """シミュレーション結果と理論解を比較しCSVへ出力する"""
        data = self.overlap_series.get(pair_id)
        if data is None:
            return

        time_array = [t - data['time'][0] for t in data['time']]  # 開始時刻を 0 とする
        x0 = data['overlap'][0]
        v0 = data['initial_vn']

        # 対象粒子を取得して有効質量を計算
        p1 = self.pars[data['p1_id']]
        p2 = self.pars[data['p2_id']]
        self.interface.set(p1, p2)
        k = self.interface.kn
        c = self.interface.etan
        m_eff = (p1.m * p2.m) / (p1.m + p2.m)

        # 理論通りのオーバーラップと速度を計算
        theo_overlap_raw = self.compute_theoretical_overlap(time_array, x0, v0, k, c, m_eff)
        theo_velocity_raw = self.compute_theoretical_velocity(time_array, x0, v0, k, c, m_eff)
        
        # プロット用に符号を調整（シミュレーションデータと合わせるため）
        theo_overlap = [-x for x in theo_overlap_raw]
        theo_velocity = [-v for v in theo_velocity_raw]

        # 平均絶対誤差を計算
        mae_overlap = sum(abs(s - t) for s, t in zip(data['overlap'], theo_overlap)) / len(theo_overlap)
        mae_velocity = sum(abs(s - t) for s, t in zip(data['velocity'], theo_velocity)) / len(theo_velocity)
        
        # 相対誤差の平均も計算
        relative_errors_overlap = []
        relative_errors_velocity = []
        for sim_o, theo_o, sim_v, theo_v in zip(data['overlap'], theo_overlap, data['velocity'], theo_velocity):
            if abs(theo_o) > 1e-12:
                rel_err_o = abs((sim_o - theo_o) / theo_o) * 100
                relative_errors_overlap.append(rel_err_o)
            if abs(theo_v) > 1e-12:
                rel_err_v = abs((sim_v - theo_v) / theo_v) * 100
                relative_errors_velocity.append(rel_err_v)
        
        avg_rel_error_overlap = sum(relative_errors_overlap) / len(relative_errors_overlap) if relative_errors_overlap else 0
        avg_rel_error_velocity = sum(relative_errors_velocity) / len(relative_errors_velocity) if relative_errors_velocity else 0

        self.overlap_analysis_results.append([pair_id, mae_overlap, mae_velocity, avg_rel_error_overlap, avg_rel_error_velocity])

        # 出力ディレクトリを用意
        overlap_dir = os.path.join('outputs', 'overlap')
        os.makedirs(overlap_dir, exist_ok=True)

        # オーバーラップと速度の比較結果をCSVに保存
        filename = os.path.join(overlap_dir, f"overlap_velocity_comparison_{pair_id}_{self.step}.csv")
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['time', 'sim_overlap', 'theoretical_overlap', 'overlap_abs_error', 'overlap_rel_error(%)', 
                           'sim_velocity', 'theoretical_velocity', 'velocity_abs_error', 'velocity_rel_error(%)'])
            for i, t in enumerate(time_array):
                s_o, t_o = data['overlap'][i], theo_overlap[i]
                s_v, t_v = data['velocity'][i], theo_velocity[i]
                abs_err_o = abs(s_o - t_o)
                abs_err_v = abs(s_v - t_v)
                rel_err_o = abs((s_o - t_o) / t_o) * 100 if abs(t_o) > 1e-12 else 0
                rel_err_v = abs((s_v - t_v) / t_v) * 100 if abs(t_v) > 1e-12 else 0
                writer.writerow([t, s_o, t_o, abs_err_o, rel_err_o, s_v, t_v, abs_err_v, rel_err_v])

        # 詳細な跳ね返り係数分析を実行
        self.analyze_detailed_restitution(pair_id, data, time_array, theo_velocity)

        print(f"オーバーラップ・速度比較結果を {filename} に保存")
        print(f"  オーバーラップMAE={mae_overlap:.3e}, 平均相対誤差={avg_rel_error_overlap:.2f}%")
        print(f"  速度MAE={mae_velocity:.3e}, 平均相対誤差={avg_rel_error_velocity:.2f}%")

    def analyze_detailed_restitution(self, pair_id, data, time_array, theo_velocity):
        """詳細な跳ね返り係数分析を行う"""
        if len(data['velocity']) < 2 or len(theo_velocity) < 2:
            return
            
        # 初速度と最終速度から跳ね返り係数を計算
        initial_velocity = data['velocity'][0]
        final_velocity = data['velocity'][-1]
        initial_velocity_theo = theo_velocity[0]
        final_velocity_theo = theo_velocity[-1]
        
        # 数値解からの跳ね返り係数
        if abs(initial_velocity) > 1e-12:
            restitution_numerical = abs(final_velocity / initial_velocity)
        else:
            restitution_numerical = 0.0
            
        # 理論解からの跳ね返り係数
        if abs(initial_velocity_theo) > 1e-12:
            restitution_theoretical = abs(final_velocity_theo / initial_velocity_theo)
        else:
            restitution_theoretical = 0.0
        
        # 設定値との比較
        set_restitution = self.e
        error_numerical = abs((restitution_numerical - set_restitution) / set_restitution) * 100 if set_restitution > 0 else 0
        error_theoretical = abs((restitution_theoretical - set_restitution) / set_restitution) * 100 if set_restitution > 0 else 0
        
        # 詳細分析結果をCSVに保存
        overlap_dir = os.path.join('outputs', 'overlap')
        restitution_filename = os.path.join(overlap_dir, f"restitution_analysis_{pair_id}_{self.step}.csv")
        
        with open(restitution_filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['分析項目', '値'])
            writer.writerow(['設定跳ね返り係数', set_restitution])
            writer.writerow(['初速度(数値)', initial_velocity])
            writer.writerow(['最終速度(数値)', final_velocity])
            writer.writerow(['跳ね返り係数(数値)', restitution_numerical])
            writer.writerow(['誤差(%)(数値)', error_numerical])
            writer.writerow(['初速度(理論)', initial_velocity_theo])
            writer.writerow(['最終速度(理論)', final_velocity_theo])
            writer.writerow(['跳ね返り係数(理論)', restitution_theoretical])
            writer.writerow(['誤差(%)(理論)', error_theoretical])
            
            # 時系列での跳ね返り係数変化を計算
            writer.writerow([''])
            writer.writerow(['時間', '速度(数値)', '速度(理論)', '瞬時跳ね返り係数(数値)', '瞬時跳ね返り係数(理論)'])
            
            for i, t in enumerate(time_array):
                v_num = data['velocity'][i]
                v_theo = theo_velocity[i]
                
                # 瞬時跳ね返り係数（初速度に対する比率）
                instant_rest_num = abs(v_num / initial_velocity) if abs(initial_velocity) > 1e-12 else 0
                instant_rest_theo = abs(v_theo / initial_velocity_theo) if abs(initial_velocity_theo) > 1e-12 else 0
                
                writer.writerow([t, v_num, v_theo, instant_rest_num, instant_rest_theo])
        
        print(f"  跳ね返り係数詳細分析を {restitution_filename} に保存")
        print(f"  数値解跳ね返り係数={restitution_numerical:.6f}, 誤差={error_numerical:.2f}%")
        print(f"  理論解跳ね返り係数={restitution_theoretical:.6f}, 誤差={error_theoretical:.2f}%")

class MainWindow(tkinter.Tk):
    def __init__(self):
        tkinter.Tk.__init__(self)
        
        self.dem = DEM()
        self.canvas = tkinter.Canvas(self, bg="white")
        self.canvas.pack(fill=tkinter.BOTH,expand=True)
        self.geometry('900x600')
        self.redraw()
        print('init')
        print(f'目標ステップ数: {self.dem.max_steps}')
        print(f'時間刻み幅: {self.dem.dt}''s')
        # CSVファイル出力用のディレクトリを作成
        os.makedirs('outputs', exist_ok=True)

    def on_closing(self):
        if hasattr(self.dem, 'csvfile'):
            self.dem.csvfile.close()  # CSVファイルが存在する場合のみ閉じる
        self.destroy()

    def redraw(self): # 描画関数
        self.canvas.delete('elem')
        w = self.canvas.winfo_width()  # キャンバスの幅を取得
        h = self.canvas.winfo_height()  # キャンバスの高さを取得
        scale = min(w, h) / 200  # スケールを計算（200は元のサイズ）
        
        for p in self.dem.pars:
            x1 = round((p.x - p.r) * scale)
            y1 = round((h - (p.y + p.r) * scale))
            x2 = round((p.x + p.r) * scale)
            y2 = round((h - (p.y - p.r) * scale))
            self.canvas.create_oval(x1, y1, x2, y2, tags='elem')
            x1 = round(p.x * scale)
            y1 = round(h - p.y * scale)
            x2 = round((p.x + p.r * math.cos(p.a)) * scale)
            y2 = round((h - (p.y + p.r * math.sin(p.a)) * scale))
            self.canvas.create_line(x1, y1, x2, y2, tags='elem')
        
        for l in self.dem.lines:
            x1 = round(l.x1 * scale)
            y1 = round(h - l.y1 * scale)
            x2 = round(l.x2 * scale)
            y2 = round(h - l.y2 * scale)
            self.canvas.create_line(x1, y1, x2, y2, width=1, tags='elem')
            
        self.update_idletasks()
        self.after(100,self.redraw)
        
    def calcloop(self): # 計算ループ関数
        self.dem.calcStep()
        if self.dem.step % 1000 == 0:
            self.st = time.time()
        self.after(0,self.calcloop)

def main():
    dem = DEM()
    print('シミュレーション開始')
    print(f'目標ステップ数: {dem.max_steps}')
    print(f'時間刻み幅: {dem.dt}''s')
    
    # メインの計算ループ
    while True:
        dem.calcStep()
        if dem.step >= dem.max_steps:
            break

if __name__ == '__main__':
    main()
