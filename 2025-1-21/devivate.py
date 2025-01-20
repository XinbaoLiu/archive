# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 00:18:40 2025

@author: iopca
"""

import sympy as sp

# 定义符号变量
λ = sp.symbols('λ', real=True, positive=True)

# 原公式分解项
A = 5
B = - (λ + 5/λ) * sp.atan(λ)
C = - (2/λ) * sp.asin(λ / sp.sqrt(1 + λ**2))
D_term = (sp.pi/2 - sp.atan(λ / sp.sqrt(2 + λ**2)))
D = (2 / (λ * sp.sqrt(2 + λ**2))) * D_term

# 组合完整公式
S3T = -1/(60 * sp.pi) * (A + B + C + D)

# 对 λ 求导
dS3T_dλ = sp.diff(S3T, λ)

# 简化结果（可能需要较长时间）
dS3T_dλ_simplified = sp.simplify(dS3T_dλ)

# 输出结果
print("导数表达式（原始形式）：")
sp.pprint(dS3T_dλ)

print("\n导数表达式（简化后）：")
sp.pprint(dS3T_dλ_simplified)

# 输出为 LaTeX 格式（可直接复制到论文中）
print("\nLaTeX 格式:")
print(sp.latex(dS3T_dλ_simplified))