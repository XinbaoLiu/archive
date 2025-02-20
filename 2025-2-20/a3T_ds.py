# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 23:23:40 2025

@author: iopca
"""

import numpy as np
import matplotlib.pyplot as plt

# 常量定义
a0 = 5.29177210903e-11  # 玻尔半径，单位：米 (m)

def calculate_lambda(n, a0):
    """计算无量纲参数 lambda"""
    return (3**(1/6)) * (np.pi**(5/6)) * np.sqrt(a0) * (n**(1/6))

def calculate_S3T(lambda_val):
    """计算无量纲参数 S3^T"""
    term1 = 5 - (lambda_val + 5/lambda_val) * np.arctan(lambda_val)
    
    arcsin_arg = lambda_val / np.sqrt(1 + lambda_val**2)
    term2 = (2 / lambda_val) * np.arcsin(arcsin_arg)
    
    denominator = lambda_val * np.sqrt(2 + lambda_val**2)
    arctan_arg = lambda_val / np.sqrt(2 + lambda_val**2)
    term3_part = (np.pi/2 - np.arctan(arctan_arg))
    term3 = (2 / denominator) * term3_part
    
    total = term1 - term2 + term3
    S3T = -total / (60 * np.pi)
    return S3T

def calculate_a3T(n, S3T):
    """计算 a3^T，单位：平方米 (m²)"""
    factor = 2 * (1/3)**(1/3) * np.pi**(-2/3)
    term = (3 / (4 * np.pi * n))**(2/3)
    a3T = factor * term * S3T
    return a3T*n

# 生成电荷密度范围（合理物理范围：1e24 到 1e30 m⁻³）
n_values = np.logspace(27, 31, 400)

# 计算 lambda、S3T、a3T
lambda_values = calculate_lambda(n_values, a0)
S3T_values = calculate_S3T(lambda_values)
a3T_values = calculate_a3T(n_values, S3T_values)

# 绘图
plt.figure(figsize=(10, 6))
plt.loglog(n_values, a3T_values, label='$a_3^T$')
plt.xlabel('Charge Density $n$ (m⁻³)', fontsize=12)
plt.ylabel('$a_3^T$ (m²)', fontsize=12)
plt.title('$a_3^T$ vs. Charge Density', fontsize=14)
plt.grid(True, which='both', linestyle='--')
plt.legend()
plt.show()
