import numpy as np
import matplotlib.pyplot as plt

from dvr_engine import n, find_grid, build_U, build_Tg, build_Vg

m_particle = 918


def V_morse(x):
    De = 0.174
    a = 1.03
    re = 1.40
    return De * (1 - np.exp(-a * (x - re))) ** 2


x_phys, beta = find_grid(V_morse, m_particle)

U_matrix = build_U()

T_grid = build_Tg(m_particle, beta)

V_grid = build_Vg(x_phys, V_morse)

H = T_grid + V_grid
vals, vecs = np.linalg.eigh(H)

print("前 5 个本征能级：")
print(vals[:5])

plt.plot(x_phys, V_morse(x_phys), "k-", lw=2, label="Morse Potential")

scale_factor = 0.05
for i in range(3):
    psi = vecs[:, i] * scale_factor + vals[i]
    plt.plot(x_phys, psi, label=f"v={i}")

plt.axhline(vals[0], color="gray", linestyle="--")

plt.ylim(0, 0.1)
plt.xlabel("r (Bohr or Angstrom)")
plt.ylabel("Energy")
plt.legend()
plt.title("1D-DVR: Morse Potential Wavefunctions")
plt.show()
