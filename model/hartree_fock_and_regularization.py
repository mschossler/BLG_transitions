from config import itmax, setH
from input.parameters import *
from model.densities_small_U import rho0constUp, rho0constUm
from model.exchange_integrals import Xskm, Xskp
from model.hamiltonians import full_hp, full_hm, full_hpm, full_hmp, idp, idps, idm, idms, delta_e_regmatrix, mZm, h0p2, h0m2, h0p, h0m, hAp
from utils import eigen, df_round, nonedimmerp, nonedimmerm


def hp(n, nprime, s1, s2):
    return full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho)


def hm(n, nprime, s1, s2):
    return full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho)


def hpm(n, nprime, s1, s2):
    return full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


def hmp(n, nprime, s1, s2):
    return full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


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
    print('u=%.2f' % (u * 1000))
    rho = rho0

    eigenvaluep2, eigenvectorp2 = eigen(h0p2(u))
    eigenvaluem2, eigenvectorm2 = eigen(h0m2(u))
    eigenvaluep = [h0p(u)[0][0]] + [h0p(u)[1][1]] + eigenvaluep2.tolist()
    eigenvaluem = [h0m(u)[0][0]] + [h0m(u)[1][1]] + eigenvaluem2.tolist()

    h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)

    # eigenvectorp = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()]))
    # eigenvectorm = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()]))

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    ###### regularization (self energy) U dependent
    # print('here1')

    # print('here2')
    ######

    ###### regularization (self energy) U dependent
    regmatrix = delta_e_regmatrix(rho0, eigenvectorp, eigenvectorm) * k * alpha_reg
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
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(occupied_bands))

        regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * k * alpha_reg

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
    Et = sum([eigenvalue[i] for i in range(nu)]) + ehf

    dict_quantities_u = {'u': u * 1e3,
                         'eigenvalue': 1e3 * eigenvalue,
                         'eigenvector': eigenvector,
                         'Et': 1e3 * Et,
                         'h0': df_round(1e3 * h0),
                         'rhoU': df_round(rho),
                         'Eh_deltaU': 1e3 * k * Eh * deltatb,
                         'Hint': df_round(1e3 * Hint),
                         }
    ############################# exciton ####################
    eigenvaluep2, eigenvectorp2 = eigen(hAp(u))[0][1:3], eigen(hAp(u))[1][1:3]
    eigenvectorp2 = nonedimmerp(eigenvectorp2)

    eigenvaluem2, eigenvectorm2 = eigen(hAp(-u))[0][1:3], eigen(hAp(-u))[1][1:3]
    eigenvectorm2 = nonedimmerm(eigenvectorm2)

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    exciton = np.array([exciton_j_to_n_km(-2, 1, eigenvectorm), exciton_j_to_n_kp(-2, 1, eigenvectorp),
                        exciton_j_to_n_km(1, 2, eigenvectorm), exciton_j_to_n_kp(1, 2, eigenvectorp)])

    dict_quantities_u['regmatrix'] = 1e3 * regmatrix
    dict_quantities_u['exciton_energy'] = 1e3 * exciton
    # dict_quantities_u['exciton'] = exciton

    return dict_quantities_u
