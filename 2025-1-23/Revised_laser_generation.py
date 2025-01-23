# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 12:49:08 2025

@author: iopca
"""

"""
Electric Field Pulse Generator for TDDFT Simulations (Revised)

Generates time-dependent electric field data in TDEFIELD.in format with explicit
unit conversion from Ry/bohr/e to V/Å. The electric field follows a Gaussian-enveloped
cosine form: E(t) = E0 * cos(2πωt + φ) * exp(-0.5*((t-t0)/δ)^2)
"""

import numpy as np
import matplotlib.pyplot as plt
from math import pi

# ========================
# 物理常数与单位转换
# ========================
RYDBERG_PER_BOHR_TO_V_PER_ANGSTROM = 36.36  # 1 Ry/bohr/e = 36.36 V/Å
INTENSITY_CONVERSION = 1.3271e13           # (V/Å)^2 → W/cm²
FS_TO_SECONDS = 1e-15                      # 飞秒到秒转换

# ========================
# 激光参数配置
# ========================
def configure_parameters():
    """返回包含所有模拟参数的字典"""
    return {
        # 时间参数
        'total_steps': 800,     # 总时间步数
        'time_step': 0.096756,   # 时间步长 (fs)
        't0': 30.0,              # 脉冲中心时间 (fs)
        'pulse_width': 7.0,      # 脉冲宽度δ (fs)
        
        # 电场参数 (Ry/bohr/e 单位)
        'E0_ry': 0.002,          # 电场振幅 0.002 Ry/bohr/e
        'frequency': 0.375,      # 载波频率 (fs⁻¹)
        'phase': 0.0,            # 相位 (弧度)
        
        # 介质参数
        'refractive_index': 2.0, # 折射率
        
        # 输出方向 (x,y,z 分量比例)
        'polarization': [0.0, 1.0, 0.0]  # Y方向极化
    }

# ========================
# 电场脉冲生成函数
# ========================
def generate_electric_field(t, params):
    """
    生成高斯包络调制的余弦电场脉冲
    
    参数:
        t (np.array): 时间序列 (fs)
        params (dict): 参数字典
    
    返回:
        np.array: 电场强度序列 (V/Å)
    """
    # 单位转换
    E0_v_per_angstrom = params['E0_ry'] * RYDBERG_PER_BOHR_TO_V_PER_ANGSTROM
    
    # 载波部分
    carrier = np.cos(2 * pi * params['frequency'] * t + params['phase'])
    
    # 高斯包络
    envelope = np.exp(-0.5 * ((t - params['t0']) / params['pulse_width'])**2)
    
    return E0_v_per_angstrom * carrier * envelope * 100000/36.36

# ========================
# 激光参数计算
# ========================
def calculate_laser_parameters(E_field, t, params):
    """计算并返回激光强度相关参数"""
    # 瞬时强度 (W/cm²)
    intensity = (
        params['refractive_index'] 
        * INTENSITY_CONVERSION 
        * (E_field ** 2)
    )
    
    # 脉冲能量密度
    energy_density = np.trapz(intensity, dx=params['time_step']*FS_TO_SECONDS)  # J/cm²
    
    return {
        'peak_intensity': np.max(intensity),
        'pulse_energy': energy_density * 1e3  # 转换为 mJ/cm²
    }

# ========================
# 主程序
# ========================
def main():
    # 初始化配置
    params = configure_parameters()
    time_points = np.arange(params['total_steps']) * params['time_step']
    
    # 生成电场
    E_field = generate_electric_field(time_points, params)
    
    # 转换为三维电场
    field_3d = np.outer(E_field, params['polarization'])
    
    # 保存文件（关键修改点）
    np.savetxt("TDEFIELD.in", field_3d, 
               fmt="    %.8f    %.8f    %.8f",
               header="Time-dependent Electric Field (V/A)",  # 移除非ASCII字符
               encoding='utf-8')  # 指定UTF-8编码
    
    # 计算激光参数
    laser_params = calculate_laser_parameters(E_field, time_points, params)
    
    # 打印结果
    print("\n=== 激光脉冲参数 ===")
    print(f"电场振幅:      {params['E0_ry']} Ry/bohr/e")
    print(f"峰值强度:      {laser_params['peak_intensity']:.8e} W/cm²")
    print(f"脉冲能量密度:  {laser_params['pulse_energy']:.8f} mJ/cm²")
    
    # 可视化
    plt.figure(figsize=(10, 6))
    plt.plot(time_points, E_field, lw=2)
    plt.xlabel('Time (fs)', fontsize=12)
    plt.ylabel('Electric Field (V/Å)', fontsize=12)
    plt.title(f'Gaussian Pulse (E0 = {params["E0_ry"]} Ry/bohr/e)', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()