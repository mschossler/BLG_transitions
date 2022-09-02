from config import itmax, setH
from input.parameters import *
from model.densities_small_U import rho0constUp, rho0constUm
# from model.density_test import rho0constUp, rho0constUm
from model.exchange_integrals import Xskm, Xskp
from model.hamiltonians import full_hp, full_hm, full_hpm, full_hmp, idp, idps, idm, idms, mZm, hAp, hBp, hCp
from utils import eigen, df_round, nonedimmerp, nonedimmerm


def hp(n, nprime, s1, s2):
    return full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho)


def hm(n, nprime, s1, s2):
    return full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho)


def hpm(n, nprime, s1, s2):
    return full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


def hmp(n, nprime, s1, s2):
    return full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


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


def loopU(u):
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

    eigenvaluep1 = eigen(hBp(u))[0][3]
    eigenvaluem1 = eigen(hBp(-u))[0][3]
    # print(eigenvaluep1,eigenvaluem1)

    eigenvaluep0 = eigen(hCp(u))[0][2]
    eigenvaluem0 = eigen(hCp(-u))[0][2]
    # print(eigenvaluep0,eigenvaluem0)

    eigenvaluep = [eigenvaluep0] + [eigenvaluep1] + eigenvaluep2.tolist()
    eigenvaluem = [eigenvaluem0] + [eigenvaluem1] + eigenvaluem2.tolist()

    #########################################################################################

    # ######################################
    #     eigenvaluep2, eigenvectorp2 = eigen(h0p2(u))
    #     eigenvaluem2, eigenvectorm2 = eigen(h0m2(u))
    #     eigenvaluep = [h0p(u)[0][0]] + [h0p(u)[1][1]] + eigenvaluep2.tolist()
    #     eigenvaluem = [h0m(u)[0][0]] + [h0m(u)[1][1]] + eigenvaluem2.tolist()
    # ######################################
    h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)  # none interacting matrix

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

    it = 1
    while it < itmax:
        deltatb = sum([rho[idp(n)][idp(n)] + rho[idps(n)][idps(n)] - rho[idm(n)][idm(n)] - rho[idms(n)][idms(n)] for n in setH])
        # hphpm = [[hp(n, nprime) for n in setH] + [hpm(n, nprime) for n in setH] for nprime in setH]
        # hmphm = [[hmp(n, nprime) for n in setH] + [hm(n, nprime) for n in setH] for nprime in setH]
        hphpm = [[hp(n, nprime, 1, 1) for nprime in setH] + [hpm(n, nprime, 1, 1) for nprime in setH] for n in setH]
        hmphm = [[hmp(n, nprime, 1, 1) for nprime in setH] + [hm(n, nprime, 1, 1) for nprime in setH] for n in setH]

        Hintup = np.vstack((hphpm, hmphm))

        hphpms = [[hp(n, nprime, -1, -1) for nprime in setH] + [hpm(n, nprime, -1, -1) for nprime in setH] for n in setH]
        hmphms = [[hmp(n, nprime, -1, -1) for nprime in setH] + [hm(n, nprime, -1, -1) for nprime in setH] for n in setH]

        Hintdown = np.vstack((hphpms, hmphms))

        hphpmsud = [[hp(n, nprime, 1, -1) for nprime in setH] + [hpm(n, nprime, 1, -1) for nprime in setH] for n in setH]
        hmphmsud = [[hmp(n, nprime, 1, -1) for nprime in setH] + [hm(n, nprime, 1, -1) for nprime in setH] for n in setH]

        Hintupdown = np.vstack((hphpmsud, hmphmsud))

        hphpmsdu = [[hp(n, nprime, -1, 1) for nprime in setH] + [hpm(n, nprime, -1, 1) for nprime in setH] for n in setH]
        hmphmsdu = [[hmp(n, nprime, -1, 1) for nprime in setH] + [hm(n, nprime, -1, 1) for nprime in setH] for n in setH]

        Hintdownup = np.vstack((hphpmsdu, hmphmsdu))

        Hint = k * np.vstack((np.hstack((Hintup, Hintupdown)), np.hstack((Hintdownup, Hintdown))))
        H = Hint + h0 + mZm + regmatrix  # np.add(Hint, h0)
        eigenvalue_loop, eigenvector_loop = eigen(H)
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(number_occupied_bands))

        # regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg # we should not update regmatrix here, make results_old for nu=0 bad and won't affect nu=4

        it += 1

    eigenvalue, eigenvector = eigen(H)
    # eigenvalue, eigenvector = npla.eig(H)
    # idxfunc = np.argsort(eigenvalue)
    #
    # eigenvalue = eigenvalue[idxfunc]
    # eigenvector = eigenvector[:, idxfunc]

    ehf = - sum([Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in setH for nprime in setH] +
                [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in setH for nprime in setH] +
                [Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in setH for nprime in setH] +
                [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in setH for nprime in setH] +
                [2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in setH for nprime in setH] +
                [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in setH for nprime in setH]
                )
    Et = sum([eigenvalue[i] for i in range(number_occupied_bands)]) + ehf

    dict_quantities_u = {'u': u * 1e3,
                         'eigenvalue': 1e3 * eigenvalue,
                         'eigenvector': eigenvector,
                         'Et': 1e3 * Et,
                         'h0': df_round(1e3 * h0),
                         'rhoU': df_round(rho),
                         'Eh_deltaU': 1e3 * k * Eh * deltatb,
                         'Hint': df_round(1e3 * Hint),
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
