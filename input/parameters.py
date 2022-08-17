import numpy as np

asym = 1
alpha_Zm = 1 #0.04227165829987071 # k / alpha_k
alpha_H_oct_int = 1
inter = 1
alpha_x = 1
alpha_k = 0.26
alpha_rho = 0
alpha_state = 1
alpha_rand=0.01
alpha_reg = 0
dens = 1

uz = 7e0 * 1e-3# * 0
uperp = -1.6e0 * 1e-3#  * 0
B = 13
omega = 35 * np.sqrt(B) / 1000
gamma0 = 3
gamma1 = 0.41
gamma4 = 0.15 * asym
eta4 = gamma4 / gamma0
gamma3 = -0.3 * asym # gamma3 = 0.3;
eta3 = gamma3 / gamma0
Delta = 0.018 * asym
Delta_ab = 0.001 * 0
Delta_td = (2 * gamma1 * gamma4 / gamma0 + Delta)
beta = (omega / gamma1) ** 2
Zm = 57.9e-6 * B * alpha_Zm   # 57.9e-6 is from Zhang2011PRB below equation 17 #temperature: 4 K=4 * 0.0862=0.3448meV
ep = 0.5
nu = 4

occupied_bands = nu + 8 # nu = 8 is the charge neutrality point in this schema
clight = 299792458
hbar = (6.62607 * 10 ** (-34)) / (2 * np.pi)
el = 1.602176634 * 10 ** (-19)
dlayer = 3.35 * 10 ** (-10)
ep0 = 8.8541878128 * 10 ** (-12)

Lb = np.sqrt(hbar / (el * B))
x = dlayer / Lb * alpha_x
epr = 6
Eh = x / np.sqrt(2 * np.pi)
k = (np.sqrt(np.pi / 2) * el) / (4 * np.pi * ep0 * epr * Lb) * alpha_k

U0minD = -2e-3
U0maxD = -1.25e-3
dU0D = 0.25e-3

