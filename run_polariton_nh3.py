import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.special import factorial, eval_hermite
from numpy.polynomial.hermite import hermgauss
from dvr_engine import n, find_grid, build_U, build_Tg, build_Vg,second_derivative,x_math

cm_to_au=1/ 219474.63

m_q=4533.52
A0, A2, A4 = 9.249e-3, -3.289e-2, 2.923e-2
gamma, delta = 1.271, 0.8887

def V_nh3(q):
    return A0+A2*q**2+A4*q**4

def d_nh3(q):
    return -(gamma)*q*np.exp(-delta*q**2)

beta_q=0.16
x_phys_q = x_math * beta_q

def find_grid_well(V_func, m):
    res = minimize(V_func, x0=1.0)
    x_min = res.x[0]
    v_min = res.fun
    k = second_derivative(V_func, x_min)
    omega = np.sqrt(k / m)
    beta = np.sqrt(1 / (m * omega))
    x_phys = x_math * beta + x_min
    return x_phys, beta

x_phys, beta = find_grid_well(V_nh3, m_q)

U_matrix=build_U()
T_q = build_Tg(m_q, beta_q)
V_q = build_Vg(x_phys_q, V_nh3)
H_q = T_q + V_q
vals_q, vecs_q = np.linalg.eigh(H_q)

omega_c=vals_q[3]-vals_q[0]
D_q=np.diag(d_nh3(x_phys_q))
d_10 = np.abs(vecs_q[:, 3].T @ D_q @ vecs_q[:, 0])

m_c=1.0

beta_c=np.sqrt(1/(m_c*omega_c))
x_phys_c = x_math * beta_c
T_c = build_Tg(m_c, beta_c)

def V_cavity(xc):
    return 0.5 * omega_c**2 * xc**2

V_c = build_Vg(x_phys_c, V_cavity)
H_c = T_c + V_c
X_c = np.diag(x_phys_c) 

I_q=np.eye(n)
I_c=np.eye(n)

H0=np.kron(I_q,H_c)+np.kron(H_q,I_c)

D_2D=np.kron(D_q, I_c)

eta_list=[0.0,0.05,0.1,0.15,0.20]
colors=['b','g','r','c','m']
E_axis = np.linspace(500, 2500, 2000)
spectra_matrix = np.zeros((len(eta_list), len(E_axis)))

for idx,eta in enumerate(eta_list):
    g = eta * omega_c / d_10
    H_int = np.sqrt(2 * omega_c) * g * np.kron(D_q, X_c)
    H_DSE = (g**2 / omega_c) * np.kron(D_q @ D_q, I_c)
    H_tot = H0 + H_int + H_DSE

    vals_tot, vecs_tot = np.linalg.eigh(H_tot)
    psi_0 = vecs_tot[:, 0]  
    E_0 = vals_tot[0]       

    spectrum_y = np.zeros_like(E_axis)
    
    
    for i in range(1, 30):
        
        delta_E = (vals_tot[i] - E_0) / cm_to_au
        
        
        psi_i = vecs_tot[:, i]
        transition_dipole = psi_i.T @ D_2D @ psi_0
        intensity = np.abs(transition_dipole)**2
        
        
        gamma_width = 10.0
        peak = intensity * (gamma_width / ((E_axis - delta_E)**2 + gamma_width**2))
        
        spectrum_y += peak
        
    spectra_matrix[idx, :] = spectrum_y

plt.figure(figsize=(10, 6))

for idx, eta in enumerate(eta_list):
    
    plt.plot(E_axis, spectra_matrix[idx, :] * 1e5, label=f'$\eta$ = {eta}', color=colors[idx], lw=1.5)

plt.xlabel(r"$\omega / \mathrm{cm}^{-1}$", fontsize=14)
plt.ylabel("Intensity [arb. units]", fontsize=14)
plt.title(r"IR Spectrum", fontsize=16)
plt.xlim(500, 2500)
plt.legend(title=r"Coupling Strength", loc='upper right')
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()



 
