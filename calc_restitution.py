#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import json
import math

def calculate_particle_mass(particle_radius, rho):
    """
    粒子の質量を計算する関数
    質量 = (4/3) * π * ρ * r^3
    """
    return (4.0 / 3.0) * math.pi * rho * (particle_radius ** 3)

def calculate_restitution_coefficient(particle_radius, rho, kn, etan):
    """
    Fortran の calculate_restitution_coefficient サブルーチンと同様の式を用いて、
    跳ね返り係数を計算する。
    restitution_coefficient = exp( -pi / (2 * sqrt( ((4/3*pi*rho*r^3)*kn/(etan^2)) - 1/4 )) )
    """
    pi = 3.14159265358979
    term = ((4.0/3.0) * pi * rho * (particle_radius**3)) * kn / (etan**2)
    inner_sqrt = term - 0.25  # (1.0 / 4.0)
    restitution_coefficient = math.exp(-pi / (2.0 * math.sqrt(inner_sqrt)))
    return restitution_coefficient

def main():
    """
    Usage:
        python calc_restitution.py config.json
    で実行すると、config.json から粒子半径、密度(rho)、kn、etan を読み込み、
    跳ね返り係数を計算して表示する。
    """
    if len(sys.argv) < 2:
        print(f"使い方: {sys.argv[0]} <configファイルパス>")
        sys.exit(1)

    config_path = sys.argv[1]

    # JSON ファイルの読み込み
    with open(config_path, 'r', encoding='utf-8') as f:
        config = json.load(f)

    # config.json 内の関連パラメータ取得
    rho = config["simulation"]["rho"]
    particle_radius = config["particles"]["radius"]
    kn = config["interface"]["particle_particle"]["kn"]
    etan = config["interface"]["particle_particle"]["etan"]

    # 粒子質量の計算
    mass = calculate_particle_mass(
        particle_radius=particle_radius,
        rho=rho
    )
    
    # 跳ね返り係数の計算
    restitution = calculate_restitution_coefficient(
        particle_radius=particle_radius,
        rho=rho,
        kn=kn,
        etan=etan
    )

    # 計算結果を表示
    print(f"粒子間の跳ね返り係数: {restitution:.6f}")
    print(f"粒子質量: {mass:.6f}")
    

if __name__ == "__main__":
    main()