from config import itmax
from input.parameters import *
from model.densities_small_U import rho0constUp, rho0constUm
from model.exchange_integrals import Xzs, Xzd, Xos, Xod, Xfs, Xfd, Xsts, Xstd
from model.hamiltonians import mZm, hAp, hBp, hCp  # , asymmetric_h, taux, tauy, tauz
from utils import eigen, nonedimmerp, nonedimmerm, tau_func

base_oct = ['0m-', '0p-', '1m-', '1p-', '0m+', '0p+', '1m+', '1p+']
base_full = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']

base_dict = {}
for base in base_oct:
    base_dict[base_oct.index(base) + 1] = base_full.index(base)

base_dict_inverse = {}
for key in base_dict:
    base_dict_inverse[base_dict[key]] = key


def Hint_oct(rhotmp):
    def rho(i, j):
        return rhotmp[base_dict[i]][base_dict[j]]

    deltatb = (rho(1, 1) - rho(2, 2) + rho(3, 3) - rho(4, 4) + rho(5, 5) - rho(6, 6) + rho(7, 7) - rho(8, 8))
    hamtx = [[- Xzs * rho(1, 1) - Xfs * rho(3, 3) + Eh * deltatb, -Xzd * rho(1, 2) - Xfd * rho(3, 4), -Xsts * rho(1, 3), -Xstd * rho(1, 4), -Xzs * rho(1, 5) - Xfs * rho(3, 7),
              -Xzd * rho(1, 6) - Xfd * rho(3, 8), -Xsts * rho(1, 7), -Xstd * rho(1, 8)],
             [-Xzd * rho(1, 2) - Xfd * rho(3, 4), - Xzs * rho(2, 2) - Xfs * rho(4, 4) - Eh * deltatb, -Xstd * rho(2, 3), -Xsts * rho(2, 4), -Xzd * rho(2, 5) - Xfd * rho(4, 7),
              -Xzs * rho(2, 6) - Xfs * rho(4, 8), -Xstd * rho(2, 7), -Xsts * rho(2, 8)],
             [-Xsts * rho(1, 3), -Xstd * rho(2, 3), - Xfs * rho(1, 1) - Xos * rho(3, 3) + Eh * deltatb, -Xfd * rho(1, 2) - Xod * rho(3, 4), - Xsts * rho(3, 5), -Xstd * rho(3, 6),
              -Xfs * rho(1, 5) - Xos * rho(3, 7), -Xfd * rho(1, 6) - Xod * rho(3, 8)],
             [-Xstd * rho(1, 4), -Xsts * rho(2, 4), -Xfd * rho(1, 2) - Xod * rho(3, 4), - Xfs * rho(2, 2) - Xos * rho(4, 4) - Eh * deltatb, -Xstd * rho(4, 5), -Xsts * rho(4, 6),
              -Xfd * rho(2, 5) - Xod * rho(4, 7), -Xfs * rho(2, 6) - Xos * rho(4, 8)],
             [-Xzs * rho(1, 5) - Xfs * rho(3, 7), -Xzd * rho(2, 5) - Xfd * rho(4, 7), -Xsts * rho(3, 5), -Xstd * rho(4, 5), - Xzs * rho(5, 5) - Xfs * rho(7, 7) + Eh * deltatb,
              -Xfd * rho(3, 5) - Xzd * rho(5, 6), -Xsts * rho(5, 7), -Xstd * rho(5, 8)],
             [-Xzd * rho(1, 6) - Xfd * rho(3, 8), -Xzs * rho(2, 6) - Xfs * rho(4, 8), -Xstd * rho(3, 6), -Xsts * rho(4, 6), -Xfd * rho(3, 5) - Xzd * rho(5, 6),
              - Xzs * rho(6, 6) - Eh * deltatb - Xfs * rho(8, 8), -Xstd * rho(6, 7), -Xsts * rho(6, 8)],
             [-Xsts * rho(1, 7), -Xstd * rho(2, 7), -Xfs * rho(1, 5) - Xos * rho(3, 7), -Xfd * rho(2, 5) - Xod * rho(4, 7), -Xsts * rho(5, 7), -Xstd * rho(6, 7),
              - Xfs * rho(5, 5) - Xos * rho(7, 7) + Eh * deltatb, -Xfd * rho(5, 6) - Xod * rho(7, 8)],
             [-Xstd * rho(1, 8), -Xsts * rho(2, 8), -Xfd * rho(1, 6) - Xod * rho(3, 8), -Xfs * rho(2, 6) - Xos * rho(4, 8), -Xstd * rho(5, 8), -Xsts * rho(6, 8),
              -Xfd * rho(5, 6) - Xod * rho(7, 8), - Xfs * rho(6, 6) - Eh * deltatb - Xos * rho(8, 8)]]
    Hint_16 = np.zeros((16, 16)).tolist()
    for key_i in base_dict_inverse:
        i = key_i
        a = base_dict_inverse[key_i] - 1
        for key_j in base_dict_inverse:
            j = key_j
            b = base_dict_inverse[key_j] - 1
            Hint_16[i][j] = hamtx[a][b]
    return np.array(Hint_16)


#########################################################################################################
def asymmetric_h(rho, tau, u):
    first = u * np.trace(tau @ rho) * (tau @ rho)
    second = - u * tau @ rho @ tau @ rho
    return first + second


taux = tau_func([[0, 1], [1, 0]])
tauy = tau_func([[0, -1], [1, 0]]) * 1j
tauz = tau_func([[1, 0], [0, -1]])


def H_asym(rho):
    H_asym_tmp = asymmetric_h(rho, taux, uperp) + asymmetric_h(rho, tauy, uperp) + asymmetric_h(rho, tauz, uz)
    return H_asym_tmp


##########################################################################################################


def loopU0(u):
    if u >= 0:
        rho0 = rho0constUp  # rhodiagUp
    else:
        rho0 = rho0constUm  # rhodiagUm
    print('u=%.2f meV' % (u * 1000))
    rho = rho0

    eigenvaluep2, eigenvectorp2 = eigen(hAp(u))[0][1:3], eigen(hAp(u))[1][1:3]
    eigenvectorp2 = nonedimmerp(eigenvectorp2)
    # print(eigenvectorp2)

    eigenvaluem2, eigenvectorm2 = eigen(hAp(-u))[0][1:3], eigen(hAp(-u))[1][1:3]
    eigenvectorm2 = nonedimmerm(eigenvectorm2)

    # print(eigenvectorm2)

    eigenvaluep2 = eigen(hAp(u))[0][1:3]
    eigenvaluem2 = eigen(hAp(-u))[0][1:3]
    # print(eigenvaluep2,eigenvaluem2)

    # Delta_ab = dab
    eigenvaluep1 = eigen(hBp(u))[0][3]
    eigenvaluem1 = eigen(hBp(-u))[0][3]
    # print(eigenvaluep1,eigenvaluem1)

    eigenvaluep0 = eigen(hCp(u))[0][2]
    eigenvaluem0 = eigen(hCp(-u))[0][2]
    # print(eigenvaluep0,eigenvaluem0)

    eigenvaluep = [eigenvaluep0] + [eigenvaluep1] + (eigenvaluep2).tolist()
    eigenvaluem = [eigenvaluem0] + [eigenvaluem1] + (eigenvaluem2).tolist()
    # bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
    h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)  # follows notation order

    it = 1;
    while it < itmax:
        Hint_longrange = k * alpha_H_oct_int * Hint_oct(rho)
        H_asym = asymmetric_h(taux, rho, uperp) + asymmetric_h(tauy, rho, uperp) + asymmetric_h(tauz, rho, uz)
        H = h0 + Hint_longrange + H_asym + mZm
        eigenvalue_loop, eigenvector_loop = eigen(H);
        rhotemp = rho
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(nu));
        rho = (1 - alpha_rho) * rho + alpha_rho * rhotemp
        it += 1

    # eigenvalue, eigenvector = npla.eig(H) # follows notation (bands list) order
    eigenvalue, eigenvector = eigen(H)
    # eigenvector = np.transpose(eigenvector)# follows notation (bands list) order
    # eigenvalue_h0, eigenvector = npla.eig(h0+mZm)
    # eigenvalue_H_asym, eigenvector = npla.eig(H_asym)
    #
    #
    # print(sorted(np.real(eigenvalue_h0)*1e3))
    # print(sorted(np.real(eigenvalue_H_asym)*1e3))
    # print(sorted(np.real(eigenvalue)*1e3))

    # return [u * 1e3]+ (np.real(eigenvalue) * 1e3).tolist()
    return u * 1e3, np.real(eigenvalue) * 1e3, eigenvector, rho
