from config import setH
from input.parameters import *
from model.densities import rho0constUp, rho0constUm
from model.exchange_integrals import Xskm, Xskp
from model.hamiltonians import hAp
from utils import eigen, df_round, nonedimmerp, nonedimmerm


##########################################################################################################
# regularization (self energy) U dependent

def delta_e_kp(n, nu0, nu1, num2, nu2, eigenvectorp):
    pzm = (nu1 - 1 / 2) * Xskp(n, 1, 1, n, eigenvectorp) + (nu0 - 1 / 2) * Xskp(n, 0, 0, n, eigenvectorp)

    m2and2 = nu2 * Xskp(n, 2, 2, n, eigenvectorp) - (1 - num2) * Xskp(n, -2, -2, n, eigenvectorp)

    F = (Xskp(n, 2, 2, n, eigenvectorp) - Xskp(n, -2, -2, n, eigenvectorp))

    #     F0=Xskp(0, 2, 2, 0, eigenvectorp) - Xskp(0, -2, -2, 0, eigenvectorp)

    #     F=F-F0

    res = (F / 2 - pzm - m2and2) * k
    return res


def delta_e_km(n, nu0, nu1, num2, nu2, eigenvectorm):
    pzm = (nu1 - 1 / 2) * Xskm(n, 1, 1, n, eigenvectorm) + (nu0 - 1 / 2) * Xskm(n, 0, 0, n, eigenvectorm)

    m2and2 = nu2 * Xskm(n, 2, 2, n, eigenvectorm) - (1 - num2) * Xskm(n, -2, -2, n, eigenvectorm)

    F = (Xskm(n, 2, 2, n, eigenvectorm) - Xskm(n, -2, -2, n, eigenvectorm))

    #     F0=Xskm(0, 2, 2, 0, eigenvectorm) - Xskm(0, -2, -2, 0, eigenvectorm)

    #     F=F-F0

    res = (F / 2 - pzm - m2and2) * k
    return res


def delta_e_regmatrix(rho0const, eigenvectorp, eigenvectorm):
    # print('here3 ', np.diag(rho0const))
    nu0kp, nu1kp, num2kp, nu2kp, nu0km, nu1km, num2km, nu2km = tuple(np.diag(rho0const)[:8])
    # print('here4')
    regmatrixUmspindown = [delta_e_kp(n, nu0kp, nu1kp, num2kp, nu2kp, eigenvectorp) for n in setH] + [delta_e_km(n, nu0km, nu1km, num2km, nu2km, eigenvectorm) for n in setH]

    nu0kp, nu1kp, num2kp, nu2kp, nu0km, nu1km, num2km, nu2km = tuple(np.diag(rho0const)[8:16])
    regmatrixUmspinup = [delta_e_kp(n, nu0kp, nu1kp, num2kp, nu2kp, eigenvectorp) for n in setH] + [delta_e_km(n, nu0km, nu1km, num2km, nu2km, eigenvectorm) for n in setH]
    return np.diag(np.array(regmatrixUmspindown + regmatrixUmspinup))


##########################################################################################################

def exciton_j_to_n_km(n, j, eigenvectorm):
    A_nj = Xskm(n, n, j, j, eigenvectorm) * k * alpha_k
    nu_j = 1
    nu_n = 0
    return -(nu_j - nu_n) * A_nj


def exciton_j_to_n_kp(n, j, eigenvectorp):
    A_nj = Xskp(n, n, j, j, eigenvectorp) * k * alpha_k
    nu_j = 1
    nu_n = 0
    return -(nu_j - nu_n) * A_nj


def regmatrix_u(u):
    global eigenvectorp, eigenvectorm, deltatb, rho
    if u >= 0:
        rho0 = rho0constUp
    else:
        rho0 = rho0constUm
    print('running nu=%(nu)i u=%(u).2fmeV' % {'u': (u * 1000), 'nu': nu})
    rho = rho0

    ################### warping #############################################################
    eigenvaluep2, eigenvectorp2 = eigen(hAp(u))[0][1:3], eigen(hAp(u))[1][1:3]
    eigenvectorp2 = nonedimmerp(eigenvectorp2)

    eigenvaluem2, eigenvectorm2 = eigen(hAp(-u))[0][1:3], eigen(hAp(-u))[1][1:3]
    eigenvectorm2 = nonedimmerm(eigenvectorm2)

    # eigenvaluep1 = eigen(hBp(u))[0][3]
    # eigenvaluem1 = eigen(hBp(-u))[0][3]
    # print(eigenvaluep1,eigenvaluem1)

    # eigenvaluep0 = eigen(hCp(u))[0][2]
    # eigenvaluem0 = eigen(hCp(-u))[0][2]
    # print(eigenvaluep0,eigenvaluem0)

    # eigenvaluep = [eigenvaluep0] + [eigenvaluep1] + eigenvaluep2.tolist()
    # eigenvaluem = [eigenvaluem0] + [eigenvaluem1] + eigenvaluem2.tolist()

    #########################################################################################

    # ######################################
    #     eigenvaluep2, eigenvectorp2 = eigen(h0p2(u))
    #     eigenvaluem2, eigenvectorm2 = eigen(h0m2(u))
    #     eigenvaluep = [h0p(u)[0][0]] + [h0p(u)[1][1]] + eigenvaluep2.tolist()
    #     eigenvaluem = [h0m(u)[0][0]] + [h0m(u)[1][1]] + eigenvaluem2.tolist()
    # ######################################
    # h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)  # none interacting matrix

    # eigenvectorp = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()]))
    # eigenvectorm = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()]))

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    eigenvectorp_none_interact = eigenvectorp
    eigenvectorm_none_interact = eigenvectorm
    ###### regularization (self energy) U dependent
    # print('here1')

    # print('here2')
    ######

    ###### regularization (self energy) U dependent
    regmatrix = delta_e_regmatrix(rho0, eigenvectorp, eigenvectorm) * alpha_reg
    ######

    dict_quantities_u = {'u': u * 1e3,
                         # 'eigenvalue': 1e3 * eigenvalue,
                         # 'eigenvector': eigenvector,
                         # 'Et': 1e3 * Et,
                         # 'h0': df_round(1e3 * h0),
                         'rhoU': df_round(rho),
                         'Eh_deltaU': 1e3 * k * Eh * deltatb,
                         # 'Hint': df_round(1e3 * Hint),
                         'regmatrix': 1e3 * regmatrix
                         }
    ############################# exciton ####################
    # eigenvaluep2, eigenvectorp2 = eigen(hAp(u))[0][1:3], eigen(hAp(u))[1][1:3]
    # eigenvectorp2 = nonedimmerp(eigenvectorp2)
    #
    # eigenvaluem2, eigenvectorm2 = eigen(hAp(-u))[0][1:3], eigen(hAp(-u))[1][1:3]
    # eigenvectorm2 = nonedimmerm(eigenvectorm2)
    #
    # eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    # eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    exciton = np.array([exciton_j_to_n_km(-2, 1, eigenvectorm_none_interact), exciton_j_to_n_kp(-2, 1, eigenvectorp_none_interact),
                        exciton_j_to_n_km(1, 2, eigenvectorm_none_interact), exciton_j_to_n_kp(1, 2, eigenvectorp_none_interact)])

    dict_quantities_u['exciton_energy'] = 1e3 * exciton
    # dict_quantities_u['exciton'] = exciton

    # for t in allowed_transitions:
    #     transition_energy()

    return dict_quantities_u
