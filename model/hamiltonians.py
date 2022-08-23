from input.parameters import np, Delta_ab, beta, Delta_td, omega, gamma1, gamma3, Delta, eta3, eta4, Zm
from model.exchange_integrals import Xskp, Xskm, Xdkpm, Xdkmp

# Zeeman energy
mZm = np.diag([-Zm for i in range(8)] + [Zm for i in range(8)])


def h0p(u):
    hamtx = [[u / 2 - Delta_ab, 0, 0, 0], [0, beta * (Delta_td - u) + u / 2 - Delta_ab, 0, 0],
             [0, 0, beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [0, 0, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6), 2 * beta * (Delta_td - u) + u / 2 - Delta_ab]]
    return hamtx


def h0p2(u):
    hamtx = [[beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [-np.sqrt(2) * beta * (gamma1 + gamma3 / 6), 2 * beta * (Delta_td - u) + u / 2 - Delta_ab]]
    return hamtx


def h0m(u):
    hamtx = [[-u / 2 + Delta_ab, 0, 0, 0], [0, beta * (Delta_td + u) - u / 2 + Delta_ab, 0, 0],
             [0, 0, 2 * beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [0, 0, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6), beta * (Delta_td - u) + u / 2 - Delta_ab]]
    return hamtx


def h0m2(u):
    hamtx = [[2 * beta * (Delta_td + u) - u / 2 + Delta_ab, -np.sqrt(2) * beta * (gamma1 + gamma3 / 6)],
             [-np.sqrt(2) * beta * (gamma1 + gamma3 / 6), beta * (Delta_td - u) + u / 2 - Delta_ab]]
    return hamtx


def H0(u):
    return [m + [0, 0, 0, 0] for m in h0p(u)] + [[0, 0, 0, 0] + m for m in h0m(u)]



def hAp(u):
    "LL: -2, 2,...  - plot color: Black"
    if u < 0:
        Delta_ab_local = - Delta_ab
    else:
        Delta_ab_local = Delta_ab

    hamtx = [[-u / 2 + Delta_ab_local, -omega, eta4 * omega, 0],
             [-omega, - u / 2 + Delta - Delta_ab_local, gamma1, np.sqrt(2) * eta4 * omega],
             [eta4 * omega, gamma1, u / 2 + Delta + Delta_ab_local, -np.sqrt(2) * omega],
             [0, np.sqrt(2) * eta4 * omega, -np.sqrt(2) * omega, u / 2 - Delta_ab_local]]
    return hamtx


def hCp(u):
    "LL: 0, -3, 3,...  - plot color: Blue"
    if u < 0:
        Delta_ab_local = - Delta_ab
    else:
        Delta_ab_local = Delta_ab

    hamtx = [[u / 2 - Delta_ab_local, eta3 * omega, 0, 0, 0],
             [eta3 * omega, -u / 2 + Delta_ab_local, -np.sqrt(2) * omega, np.sqrt(2) * eta4 * omega, 0],
             [0, -np.sqrt(2) * omega, - u / 2 + Delta - Delta_ab_local, gamma1, np.sqrt(3) * eta4 * omega],
             [0, np.sqrt(2) * eta4 * omega, gamma1, u / 2 + Delta + Delta_ab_local, -np.sqrt(3) * omega],
             [0, 0, np.sqrt(3) * eta4 * omega, -np.sqrt(3) * omega, u / 2 - Delta_ab_local]]
    return hamtx


def hBp(u):
    "LL: 1, -4, 4...  - plot color: Red"
    if u < 0:
        Delta_ab_local = - Delta_ab
    else:
        Delta_ab_local = Delta_ab

    hamtx = [[- u / 2 + Delta - Delta_ab_local, gamma1, eta4 * omega, 0, 0, 0, 0],
             [gamma1, u / 2 + Delta + Delta_ab_local, -omega, 0, 0, 0, 0],
             [eta4 * omega, - omega, u / 2 - Delta_ab_local, np.sqrt(2) * eta3 * omega, 0, 0, 0],
             [0, 0, np.sqrt(2) * eta3 * omega, -u / 2 + Delta_ab_local, -np.sqrt(3) * omega, np.sqrt(3) * eta4 * omega, 0],
             [0, 0, 0, -np.sqrt(3) * omega, - u / 2 + Delta - Delta_ab_local, gamma1, 2 * eta4 * omega],
             [0, 0, 0, np.sqrt(3) * eta4 * omega, gamma1, u / 2 + Delta + Delta_ab_local, -2 * omega],
             [0, 0, 0, 0, 2 * eta4 * omega, - 2 * omega, u / 2 - Delta_ab_local]]

    return hamtx


idp_dic = {0: 0, 1: 1, -2: 2, 2: 3}


def idp(n):
    return idp_dic[n]


idm_dic = {0: 4, 1: 5, -2: 6, 2: 7}


def idm(n):
    return idm_dic[n]


idps_dic = {0: 8, 1: 9, -2: 10, 2: 11}


def idps(n):
    return idps_dic[n]


idms_dic = {0: 12, 1: 13, -2: 14, 2: 15}


def idms(n):
    return idms_dic[n]


def capacitance_energy_func(n, nprime, s1, s2, Eh, deltatb):
    if (n == nprime) and (s1 == s2):
        return Eh * deltatb
    return 0


def id_i(n, s, idx, idxs):
    idxfsu = round((1 + s) / 2)
    idxfsd = round((1 - s) / 2)
    return idxfsu * idx(n) + idxfsd * idxs(n)


def full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho):
    electrostatic_energy = capacitance_energy_func(n, nprime, s1, s2, Eh, deltatb)

    def id1(n):
        return id_i(n, s1, idp, idps)

    def id2(n):
        return id_i(n, s2, idp, idps)

    res = -sum(
        [Xskp(n2, nprime, n, n1, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [Xskp(n2, nprime, n, n1, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in
                                                                                                             [-2, 2] for n2 in [-2, 2]]) - electrostatic_energy
    return res


def full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho):
    electrostatic_energy = capacitance_energy_func(n, nprime, s1, s2, Eh, deltatb)

    def id1(n):
        return id_i(n, s1, idm, idms)

    def id2(n):
        return id_i(n, s2, idm, idms)

    res = -sum(
        [Xskm(n2, nprime, n, n1, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [Xskm(n2, nprime, n, n1, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in
                                                                                                             [-2, 2] for n2 in [-2, 2]]) + electrostatic_energy
    return res


def full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho):
    def id1(n):
        return id_i(n, s1, idp, idps)

    def id2(n):
        return id_i(n, s2, idm, idms)

    res = -sum([Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
        Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [-2, 2] for n2 in [-2, 2]])
    return res


def full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho):
    def id1(n):
        return id_i(n, s1, idm, idms)

    def id2(n):
        return id_i(n, s2, idp, idps)

    res = -sum([Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
        Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [-2, 2] for n2 in [-2, 2]])
    return res
