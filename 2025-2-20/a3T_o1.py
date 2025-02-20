import numpy as np
import matplotlib.pyplot as plt

# 玻尔半径（以米为单位）
a0 = 0.529e-10  # 约 0.529 Å

def S3_T(lambda_):
    """
    按照给定公式计算 S3^T 的值。
    公式:
      S3^T = -(1 / (60*pi)) * { 5 
                    - (lambda + 5/lambda)*arctan(lambda) 
                    - 2/lambda * arcsin(lambda / sqrt(1 + lambda^2)) 
                    + [2 / (lambda*sqrt(2 + lambda^2)) ] * [ pi/2 - arctan(lambda / sqrt(2 + lambda^2)) ]
                }
    其中 lambda = 3^(1/6)*pi^(5/6)*a0^(1/2)*n^(1/6)
    """
    from math import pi, atan, asin, sqrt
    
    term1 = 5
    term2 = (lambda_ + 5.0/lambda_) * np.arctan(lambda_)
    term3 = 2.0/lambda_ * np.arcsin(lambda_ / np.sqrt(1.0 + lambda_**2))
    
    # 后面那个大括号里的部分
    term4 = (2.0 / (lambda_ * np.sqrt(2.0 + lambda_**2))) * (
        pi/2 - np.arctan(lambda_ / np.sqrt(2.0 + lambda_**2))
    )
    
    # 把上面四项合并 (注意有正负号)
    inside_bracket = term1 - term2 - term3 + term4
    
    return -(1.0 / (60.0 * pi)) * inside_bracket

def a3_T(n):
    """
    计算 a3^T = 2 * (1/3)^(1/3) * pi^(-2/3) * (3/(4*pi*n))^(2/3) * S3^T
    """
    from math import pi
    
    # 先计算lambda
    lambda_ = (3.0**(1.0/6.0)) * (pi**(5.0/6.0)) * (a0**0.5) * (n**(1.0/6.0))
    
    # S3^T
    s3 = S3_T(lambda_)
    
    # 前面的常数因子
    prefactor = 2.0 * (1.0/3.0)**(1.0/3.0) * (pi**(-2.0/3.0))
    
    # (3 / (4*pi*n))^(2/3)
    bracket = (3.0 / (4.0 * pi * n))**(2.0/3.0)
    
    return prefactor * bracket * s3 * n

# === 主程序演示 ===

# 在这里选取一个电荷密度 n 的范围（单位 1/m^3 举例）
n_vals = np.logspace(27, 31, 200)  # 在 1e27 到 1e31 之间取若干点示范
a3_vals = [a3_T(nv) for nv in n_vals]

plt.figure(figsize=(6,4))
plt.loglog(n_vals, a3_vals, label=r"$a_3^T(n)$")  # 用对数坐标来展示

plt.xlabel(r"电荷密度 $n$ (m$^{-3}$)")
plt.ylabel(r"$a_3^T$ (m$^{2}$)")  # 量纲为面积
plt.title("a3^T 相对于电荷密度 n 的变化")
plt.legend()
plt.tight_layout()
plt.show()
