import numpy as np

asym = 1
B = 13
omega = 35 * np.sqrt(B) / 1000
gamma1 = 0.41
gamma4 = 0.15 * asym
# gamma3 = 0.3;
gamma3 = -0.3
gamma0 = 3
Delta = 0.018 * asym
Delta_ab = -0.001 * 0
Delta_td = (2 * gamma1 * gamma4 / gamma0 + Delta)
beta = (omega / gamma1) ** 2
Zm = 0.0178  # (57.9*(10^-6)*13)/k ,  57.9*(10^-6)*13=0.0007527eV=0.7527meV
# temperature: 4 K=4 * 0.0862=0.3448meV
# Zm = 10
ep = 0.5

nu = 8

clight = 299792458
hbar = (6.62607 * 10 ** (-34)) / (2 * np.pi)
el = 1.602176634 * 10 ** (-19)
dlayer = 3.35 * 10 ** (-10)
ep0 = 8.8541878128 * 10 ** (-12)

Lb = np.sqrt(hbar / (el * B))
x = dlayer / Lb
epr = 6
Eh = x / np.sqrt(2 * np.pi)

k = (np.sqrt(np.pi / 2) * el) / (4 * np.pi * ep0 * epr * Lb)
setH = [0, 1, -2, 2]

U0minD = 20e-3
U0maxD = 20.5e-3
dU0D = 0.25e-3


inter = 1
