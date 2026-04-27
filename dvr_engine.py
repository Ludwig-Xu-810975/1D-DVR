import numpy as np
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import minimize
from scipy.special import factorial, eval_hermite
import math

n = 100

x_math, w = hermgauss(n)
U = np.zeros((n, n))
T = np.zeros((n, n))


def second_derivative(func, x, h=1e-6):
    return (func(x + h) - 2 * func(x) + func(x - h)) / (h**2)


def find_grid(V_func, m):
    res = minimize(V_func, x0=0)
    x_min = res.x[0]
    v_min = res.fun
    k = second_derivative(V_func, x_min)
    omega = np.sqrt(k / m)
    beta = np.sqrt(1 / (m * omega))
    x_phys = x_math * beta + x_min
    return x_phys, beta


def build_U():
    for j in range(n):
        factor = 1 / np.sqrt(np.sqrt(np.pi) * (2**j) * factorial(j))
        H_val = eval_hermite(j, x_math)
        U[j, :] = np.sqrt(w) * factor * H_val

    return U


def build_Tg(m_particle, beta):
    for j in range(n):
        T[j, j] = (1 / 4) * (2 * j + 1)
        if j + 2 < n:
            T[j, j + 2] = -(1 / 4) * math.sqrt((j + 1) * (j + 2))
            T[j + 2, j] = -(1 / 4) * math.sqrt((j + 1) * (j + 2))

    Tg = U.conj().T @ T @ U
    Tg = (1 / (m_particle * beta**2)) * Tg
    return Tg


def build_Vg(x_phys, V_func):
    v = V_func(x_phys)
    Vg = np.diag(v)
    return Vg

