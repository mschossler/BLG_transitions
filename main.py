import multiprocessing
import time

import pandas as pd
from numpy import linalg as npla

from config import aux_dir_path, namecsv, itmax, nprocesses, t0, bands
from input.parameters import *
from model.densities import rho0
from model.hamiltonians import full_hp, full_hm, full_hpm, full_hmp, h0p2, h0m2, h0p, h0m, mZm, idp, idps, idm, idms
from utils import eigen, frange, df_round, sort_dict, observable_to_csv, idxcalc


def hp(n, nprime, s1, s2):
    return full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho)


def hm(n, nprime, s1, s2):
    return full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho)


def hpm(n, nprime, s1, s2):
    return full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


def hmp(n, nprime, s1, s2):
    return full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


def loopU(u):
    global eigenvectorp, eigenvectorm, deltatb, rho
    print('u=%.2f' % (u * 1000))
    rho = rho0

    eigenvaluep2, eigenvectorp2 = eigen(h0p2(u))
    eigenvaluem2, eigenvectorm2 = eigen(h0m2(u))
    eigenvaluep = [h0p(u)[0][0]] + [h0p(u)[1][1]] + eigenvaluep2.tolist()
    eigenvaluem = [h0m(u)[0][0]] + [h0m(u)[1][1]] + eigenvaluem2.tolist()

    h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    it = 1;
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
        H = Hint + h0 + mZm  # np.add(Hint, h0)
        eigenvalue_loop, eigenvector_loop = eigen(H);
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(nu));

        it += 1

    eigenvalue, eigenvector = npla.eig(H)
    idxfunc = np.argsort(eigenvalue)

    eigenvalue = eigenvalue[idxfunc]
    eigenvector = eigenvector[:, idxfunc]

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
    return dict_quantities_u


a_pool = multiprocessing.Pool(processes=nprocesses)
quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))

quantities_dict = {}
for dict in quantities:
    quantities_dict[dict['u']] = dict

quantities_dict = sort_dict(quantities_dict)

energies = []
for k, v in quantities_dict.items():
    u_temp, eigenvalue_temp, eigenvector_temp = v['u'], v['eigenvalue'], v['eigenvector']
    idx = idxcalc(eigenvector_temp)
    energies.append([u_temp] + eigenvalue_temp[idx].tolist())

energies_df = pd.DataFrame(energies, columns=['u'] + bands)
energies_df.to_csv(aux_dir_path + namecsv, index=False)

for quantity in ['h0', 'rhoU', 'Eh_deltaU', 'Hint', 'Et']:
    observable_to_csv(quantities_dict, quantity)

print(time.time() - t0)
print('file ' + namecsv + ' saved')
print(" \n done in ", time.time() - t0)