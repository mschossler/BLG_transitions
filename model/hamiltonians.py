from input.parameters import np, Delta_ab, beta, Delta_td, omega, gamma1, gamma3, Delta, eta3, eta4, Zm, Eh, uperp, uz
from model.exchange_integrals import Xskp, Xskm, Xdkpm, Xdkmp, Xzs, Xzd, Xos, Xod, Xfs, Xfd, Xsts, Xstd
from utils import tau_func
# from __main__ import Eh, deltatb, eigenvectorp, eigenvectorm, rho

def h0p(u):
    hamtx = [[u / 2 - Delta_ab, 0, 0, 0],
             [0, beta * (Delta_td - u) + u / 2 - Delta_ab, 0, 0],
             [0, 0, beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [0, 0, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6), 2 * beta * (Delta_td - u) + u / 2 - Delta_ab]]
    # hamtx=[[round(hamtx[i][j],4) for i in range(4)] for j in range(4)]
    return hamtx


def h0p2(u):
    hamtx = [[beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [-np.sqrt(2) * beta * (gamma1 + gamma3 / 6), 2 * beta * (Delta_td - u) + u / 2 - Delta_ab]]
    return hamtx


def h0m(u):
    hamtx = [[-u / 2 + Delta_ab, 0, 0, 0],
             [0, beta * (Delta_td + u) - u / 2 + Delta_ab, 0, 0],
             [0, 0, 2 * beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [0, 0, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6), beta * (Delta_td - u) + u / 2 - Delta_ab]]
    # hamtx=[[round(hamtx[i][j],4) for i in range(4)] for j in range(4)]
    return hamtx


def h0m2(u):
    hamtx = [[2 * beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [-np.sqrt(2) * beta * (gamma1 + gamma3 / 6), beta * (Delta_td - u) + u / 2 - Delta_ab]]
    return hamtx

def H0(u):
    return [m + [0, 0, 0, 0] for m in h0p(u)] + [[0, 0, 0, 0] + m for m in h0m(u)]

def hAp(u,Delta_ab):
    "LL: -2, 2,...  - plot color: Black"
    hamtx = [[-u / 2 + Delta_ab, -omega, eta4 * omega, 0],
             [-omega, - u / 2 + Delta - Delta_ab, gamma1, np.sqrt(2) * eta4 * omega],
             [eta4 * omega, gamma1, u / 2 + Delta + Delta_ab, -np.sqrt(2) * omega],
             [0, np.sqrt(2) * eta4 * omega, -np.sqrt(2) * omega, u / 2 - Delta_ab]]
    return hamtx

def hCp(u,Delta_ab):
    "LL: 0, -3, 3,...  - plot color: Blue"
    hamtx = [[u / 2 - Delta_ab, eta3 * omega, 0, 0, 0],
             [eta3 * omega, -u / 2 + Delta_ab, -np.sqrt(2) * omega, np.sqrt(2) * eta4 * omega, 0],
             [0, -np.sqrt(2) * omega, - u / 2 + Delta - Delta_ab, gamma1, np.sqrt(3) * eta4 * omega],
             [0, np.sqrt(2) * eta4 * omega, gamma1, u / 2 + Delta + Delta_ab, -np.sqrt(3) * omega],
             [0, 0, np.sqrt(3) * eta4 * omega, -np.sqrt(3) * omega, u / 2 - Delta_ab]]
    return hamtx

def hBp(u,Delta_ab):
    "LL: 1, -4, 4...  - plot color: Red"
    hamtx = [[- u / 2 + Delta - Delta_ab, gamma1, eta4 * omega, 0, 0, 0, 0],
             [gamma1, u / 2 + Delta + Delta_ab, -omega, 0, 0, 0, 0],
             [eta4 * omega, - omega, u / 2 - Delta_ab, np.sqrt(2) * eta3 * omega, 0, 0, 0],
             [0, 0, np.sqrt(2) * eta3 * omega, -u / 2 + Delta_ab, -np.sqrt(3) * omega, np.sqrt(3) * eta4 * omega, 0],
             [0, 0, 0, -np.sqrt(3) * omega, - u / 2 + Delta - Delta_ab, gamma1, 2 * eta4 * omega],
             [0, 0, 0, np.sqrt(3) * eta4 * omega, gamma1, u / 2 + Delta + Delta_ab, -2 * omega],
             [0, 0, 0, 0, 2 * eta4 * omega, - 2 * omega, u / 2 - Delta_ab]]

    return hamtx

idp_dic = {0:0,1:1,-2:2,2:3}
def idp(n):
    return idp_dic[n]

idm_dic = {0:4,1:5,-2:6,2:7}
def idm(n):
    return idm_dic[n]

idps_dic = {0:8,1:9,-2:10,2:11}
def idps(n):
    return idps_dic[n]

idms_dic = {0:12,1:13,-2:14,2:15}
def idms(n):
    return idms_dic[n]


def capacitance_energy_func(n, nprime, s1, s2, Eh, deltatb):
    if (n == nprime) and (s1 == s2):
        return Eh * deltatb
    return 0

def id_i(n,s,idx,idxs):
    idxfsu = round((1 + s) / 2)
    idxfsd = round((1 - s) / 2)
    return idxfsu * idx(n) + idxfsd * idxs(n)

def full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho):
    electrostatic_energy = capacitance_energy_func(n, nprime, s1, s2, Eh, deltatb)

    def id1(n):
        return id_i(n,s1,idp,idps)
    def id2(n):
        return id_i(n,s2,idp,idps)

    res = -sum(
        [Xskp(n2, nprime, n, n1, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
            Xskp(n2, nprime, n, n1, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in
            [-2, 2] for n2 in [-2, 2]]) - electrostatic_energy
    return res


def full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho):
    electrostatic_energy = capacitance_energy_func(n, nprime, s1, s2, Eh, deltatb)

    def id1(n):
        return id_i(n,s1,idm,idms)
    def id2(n):
        return id_i(n,s2,idm,idms)

    res = -sum(
        [Xskm(n2, nprime, n, n1, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
            Xskm(n2, nprime, n, n1, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in
            [-2, 2] for n2 in [-2, 2]]) + electrostatic_energy
    return res


def full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho):
    def id1(n):
        return id_i(n,s1,idp,idps)
    def id2(n):
        return id_i(n,s2,idm,idms)

    res = -sum(
        [Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in
         [0, 1]] + [
            Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [-2, 2] for n2 in
            [-2, 2]])
    return res


def full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho):
    def id1(n):
        return id_i(n,s1,idm,idms)
    def id2(n):
        return id_i(n,s2,idp,idps)

    res = -sum(
        [Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in
         [0, 1]] + [
            Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [-2, 2] for n2 in
            [-2, 2]])
    return res

mZm = np.diag([-Zm for i in range(8)] + [Zm for i in range(8)])

def asymmetric_h( rho, tau, u):
    first = u * np.trace(tau @ rho) * (tau @ rho)
    second = - u * tau @ rho @ tau @ rho
    return first + second



taux = tau_func([[0, 1], [1, 0]])
tauy = tau_func([[0, -1], [1, 0]]) * 1j
tauz = tau_func([[1, 0], [0, -1]])

def H_asym(rho):
    H_asym_tmp = asymmetric_h(rho, taux, uperp) + asymmetric_h(rho, tauy, uperp) + asymmetric_h(rho, tauz, uz)
    return H_asym_tmp



base_oct = ['0m-', '0p-', '1m-', '1p-','0m+', '0p+', '1m+', '1p+']
base_full = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']

base_dict = {}
for base in base_oct:
#     print(base_oct.index(base)+1,base_full.index(base))
    base_dict[base_oct.index(base)+1] = base_full.index(base)

base_dict_inverse = {}
for key in base_dict:
    base_dict_inverse[base_dict[key]] = key

# base_dict_inverse
def Hint_oct(rhotmp):
    def rho(i, j):
        return rhotmp[base_dict[i]][base_dict[j]]

    deltatb = (rho(1, 1) - rho(2, 2) + rho(3, 3) - rho(4, 4) + rho(5, 5) - rho(6, 6) + rho(7, 7) - rho(8, 8))
    hamtx = [[- Xzs * rho(1, 1) - Xfs * rho(3, 3) + Eh * deltatb, -Xzd * rho(1, 2) - Xfd * rho(3, 4), -Xsts * rho(1, 3), -Xstd * rho(1, 4), -Xzs * rho(1, 5) - Xfs * rho(3, 7), -Xzd * rho(1, 6) - Xfd * rho(3, 8), -Xsts * rho(1, 7), -Xstd * rho(1, 8)],
             [-Xzd * rho(1, 2) - Xfd * rho(3, 4), - Xzs * rho(2, 2) - Xfs * rho(4, 4) - Eh * deltatb, -Xstd * rho(2, 3), -Xsts * rho(2, 4), -Xzd * rho(2, 5) - Xfd * rho(4, 7), -Xzs * rho(2, 6) - Xfs * rho(4, 8), -Xstd * rho(2, 7), -Xsts * rho(2, 8)],
             [-Xsts * rho(1, 3), -Xstd * rho(2, 3), - Xfs * rho(1, 1) - Xos * rho(3, 3) + Eh * deltatb, -Xfd * rho(1, 2) - Xod * rho(3, 4), - Xsts * rho(3, 5), -Xstd * rho(3, 6), -Xfs * rho(1, 5) - Xos * rho(3, 7), -Xfd * rho(1, 6) - Xod * rho(3, 8)],
             [-Xstd * rho(1, 4), -Xsts * rho(2, 4), -Xfd * rho(1, 2) - Xod * rho(3, 4), - Xfs * rho(2, 2) - Xos * rho(4, 4) - Eh * deltatb, -Xstd * rho(4, 5), -Xsts * rho(4, 6), -Xfd * rho(2, 5) - Xod * rho(4, 7), -Xfs * rho(2, 6) - Xos * rho(4, 8)],
             [-Xzs * rho(1, 5) - Xfs * rho(3, 7), -Xzd * rho(2, 5) - Xfd * rho(4, 7), -Xsts * rho(3, 5), -Xstd * rho(4, 5), - Xzs * rho(5, 5) - Xfs * rho(7, 7) + Eh * deltatb, -Xfd * rho(3, 5) - Xzd * rho(5, 6), -Xsts * rho(5, 7), -Xstd * rho(5, 8)],
             [-Xzd * rho(1, 6) - Xfd * rho(3, 8), -Xzs * rho(2, 6) - Xfs * rho(4, 8), -Xstd * rho(3, 6), -Xsts * rho(4, 6), -Xfd * rho(3, 5) - Xzd * rho(5, 6), - Xzs * rho(6, 6) - Eh * deltatb - Xfs * rho(8, 8), -Xstd * rho(6, 7), -Xsts * rho(6, 8)],
             [-Xsts * rho(1, 7), -Xstd * rho(2, 7), -Xfs * rho(1, 5) - Xos * rho(3, 7), -Xfd * rho(2, 5) - Xod * rho(4, 7), -Xsts * rho(5, 7), -Xstd * rho(6, 7), - Xfs * rho(5, 5) - Xos * rho(7, 7) + Eh * deltatb, -Xfd * rho(5, 6) - Xod * rho(7, 8)],
             [-Xstd * rho(1, 8), -Xsts * rho(2, 8), -Xfd * rho(1, 6) - Xod * rho(3, 8), -Xfs * rho(2, 6) - Xos * rho(4, 8), -Xstd * rho(5, 8), -Xsts * rho(6, 8), -Xfd * rho(5, 6) - Xod * rho(7, 8), - Xfs * rho(6, 6) - Eh * deltatb - Xos * rho(8, 8)]]
    Hint_16 = np.zeros((16, 16)).tolist()

    for key_i in base_dict_inverse:
        i = key_i
        a = base_dict_inverse[key_i]-1
        for key_j in base_dict_inverse:
            j = key_j
            b = base_dict_inverse[key_j]-1
            Hint_16[i][j] = hamtx[a][b]

    return np.array(Hint_16)

########################################################################