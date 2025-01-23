import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

def eta_integrand_term1(omega, t_diff, omega_pl, rho):
    x = omega/omega_pl
    a3 = 0.008191
    b3 = 0.1589
    term1 = a3 / (1 + b3 * (x/2)**2)**(5/4)
    return -0.5 * rho**2 * term1 * np.exp(-1j*omega*t_diff)

def eta_integrand_term2(omega, t_diff, omega_pl, rho):
    x = omega/omega_pl
    gamma3 = 1.380
    omega3 = -1.070
    term2 = (x/2)**2 * np.exp(-(np.abs(x/2) + np.abs(omega3))**2 / gamma3)
    return -0.5 * rho**2 * term2 * np.exp(-1j*omega*t_diff)

def calculate_eta0(t, omega_pl=1.0, rho=1.0, t_cutoff=10):
    eta0_term1 = integrate.dblquad(
        lambda omega, tp: eta_integrand_term1(omega, t-tp, omega_pl, rho).real,
        t-t_cutoff, t,
        lambda x: 0,
        lambda x: 10*omega_pl
    )[0] / (2*np.pi)
    
    eta0_term2 = integrate.dblquad(
        lambda omega, tp: eta_integrand_term2(omega, t-tp, omega_pl, rho).real,
        t-t_cutoff, t,
        lambda x: 0,
        lambda x: 10*omega_pl
    )[0] / (2*np.pi)
    
    return eta0_term1, eta0_term2

def eta_term1(t_diff, omega_pl, rho):
    omega = np.linspace(0, 10*omega_pl, 1000)
    x = omega/omega_pl
    a3 = 0.008191
    b3 = 0.1589
    integrand = -0.5 * rho**2 * a3 / (1 + b3 * (x/2)**2)**(5/4) * np.exp(-1j*omega*t_diff)
    return np.trapz(integrand.real, omega) / (2*np.pi)

def eta_term2(t_diff, omega_pl, rho):
    omega = np.linspace(0, 10*omega_pl, 1000)
    x = omega/omega_pl
    gamma3 = 1.380
    omega3 = -1.070
    integrand = -0.5 * rho**2 * (x/2)**2 * np.exp(-(np.abs(x/2) + np.abs(omega3))**2 / gamma3) * np.exp(-1j*omega*t_diff)
    return np.trapz(integrand.real, omega) / (2*np.pi)

# Plot η₀(t)
t_values = np.linspace(0, 10, 50)
eta0_term1_values = []
eta0_term2_values = []
eta0_total_values = []

for t in t_values:
    print(f"Calculating η₀ for t = {t:.1f}")
    term1, term2 = calculate_eta0(t)
    eta0_term1_values.append(term1)
    eta0_term2_values.append(term2)
    eta0_total_values.append(term1 + term2)

plt.figure(figsize=(12, 10))

# Plot η₀(t)
plt.subplot(2, 1, 1)
plt.plot(t_values, eta0_term1_values, 'b--', label='First term')
plt.plot(t_values, eta0_term2_values, 'r--', label='Second term')
plt.plot(t_values, eta0_total_values, 'k-', label='Total')
plt.xlabel('t')
plt.ylabel('η₀(t)')
plt.title('Time integral of memory function η₀ vs time')
plt.legend()
plt.grid(True)

# Plot η(t-t') for different t
plt.subplot(2, 1, 2)
t_diff_values = np.linspace(0, 10, 200)
for t in [1, 3, 8]:
    eta1_values = [eta_term1(td, 1.0, 1.0) for td in t_diff_values]
    eta2_values = [eta_term2(td, 1.0, 1.0) for td in t_diff_values]
    total_values = np.array(eta1_values) + np.array(eta2_values)
    
    plt.plot(t_diff_values, total_values, '-', label=f't = {t}')

plt.xlabel('t-t\'')
plt.ylabel('η(t-t\')')
plt.title('Memory function η vs t-t\' at different times')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
