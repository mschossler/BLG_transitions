from config import itmax, bands, bands_oct
from input.parameters import *
from model.densities import rho0constUp, rho0constUm
from model.exchange_integrals import Xzs, Xzd, Xos, Xod, Xfs, Xfd, Xsts, Xstd
from model.hamiltonians import mZm, hAp, hBp, hCp  # , asymmetric_h, taux, tauy, tauz
from utils import eigen, nonedimmerp, nonedimmerm, tau_func, df_round

# itmax = 10000

# base_oct = ['0m-', '0p-', '1m-', '1p-', '0m+', '0p+', '1m+', '1p+']
# bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
# 
# map_octet_to_all_LL_index = {}
# for base in base_oct:
#     map_octet_to_all_LL_index[base_oct.index(base) + 1] = bands.index(base)
#
# map_octet_to_all_LL_index_inverse = {}
# for key in map_octet_to_all_LL_index:
#     map_octet_to_all_LL_index_inverse[map_octet_to_all_LL_index[key]] = key
#
# def Hint_oct(rhotmp):
#     def rho(i, j):
#         return rhotmp[map_octet_to_all_LL_index[i]][map_octet_to_all_LL_index[j]]
#
#     deltatb = (rho(1, 1) - rho(2, 2) + rho(3, 3) - rho(4, 4) + rho(5, 5) - rho(6, 6) + rho(7, 7) - rho(8, 8))
#     hamtx = [[- Xzs * rho(1, 1) - Xfs * rho(3, 3) + Eh * deltatb, -Xzd * rho(1, 2) - Xfd * rho(3, 4), -Xsts * rho(1, 3), -Xstd * rho(1, 4), -Xzs * rho(1, 5) - Xfs * rho(3, 7),
#               -Xzd * rho(1, 6) - Xfd * rho(3, 8), -Xsts * rho(1, 7), -Xstd * rho(1, 8)],
#              [-Xzd * rho(1, 2) - Xfd * rho(3, 4), - Xzs * rho(2, 2) - Xfs * rho(4, 4) - Eh * deltatb, -Xstd * rho(2, 3), -Xsts * rho(2, 4), -Xzd * rho(2, 5) - Xfd * rho(4, 7),
#               -Xzs * rho(2, 6) - Xfs * rho(4, 8), -Xstd * rho(2, 7), -Xsts * rho(2, 8)],
#              [-Xsts * rho(1, 3), -Xstd * rho(2, 3), - Xfs * rho(1, 1) - Xos * rho(3, 3) + Eh * deltatb, -Xfd * rho(1, 2) - Xod * rho(3, 4), - Xsts * rho(3, 5), -Xstd * rho(3, 6),
#               -Xfs * rho(1, 5) - Xos * rho(3, 7), -Xfd * rho(1, 6) - Xod * rho(3, 8)],
#              [-Xstd * rho(1, 4), -Xsts * rho(2, 4), -Xfd * rho(1, 2) - Xod * rho(3, 4), - Xfs * rho(2, 2) - Xos * rho(4, 4) - Eh * deltatb, -Xstd * rho(4, 5), -Xsts * rho(4, 6),
#               -Xfd * rho(2, 5) - Xod * rho(4, 7), -Xfs * rho(2, 6) - Xos * rho(4, 8)],
#              [-Xzs * rho(1, 5) - Xfs * rho(3, 7), -Xzd * rho(2, 5) - Xfd * rho(4, 7), -Xsts * rho(3, 5), -Xstd * rho(4, 5), - Xzs * rho(5, 5) - Xfs * rho(7, 7) + Eh * deltatb,
#               -Xfd * rho(3, 5) - Xzd * rho(5, 6), -Xsts * rho(5, 7), -Xstd * rho(5, 8)],
#              [-Xzd * rho(1, 6) - Xfd * rho(3, 8), -Xzs * rho(2, 6) - Xfs * rho(4, 8), -Xstd * rho(3, 6), -Xsts * rho(4, 6), -Xfd * rho(3, 5) - Xzd * rho(5, 6),
#               - Xzs * rho(6, 6) - Eh * deltatb - Xfs * rho(8, 8), -Xstd * rho(6, 7), -Xsts * rho(6, 8)],
#              [-Xsts * rho(1, 7), -Xstd * rho(2, 7), -Xfs * rho(1, 5) - Xos * rho(3, 7), -Xfd * rho(2, 5) - Xod * rho(4, 7), -Xsts * rho(5, 7), -Xstd * rho(6, 7),
#               - Xfs * rho(5, 5) - Xos * rho(7, 7) + Eh * deltatb, -Xfd * rho(5, 6) - Xod * rho(7, 8)],
#              [-Xstd * rho(1, 8), -Xsts * rho(2, 8), -Xfd * rho(1, 6) - Xod * rho(3, 8), -Xfs * rho(2, 6) - Xos * rho(4, 8), -Xstd * rho(5, 8), -Xsts * rho(6, 8),
#               -Xfd * rho(5, 6) - Xod * rho(7, 8), - Xfs * rho(6, 6) - Eh * deltatb - Xos * rho(8, 8)]]
#     Hint_16 = np.zeros((16, 16)).tolist()
#     for key_i in map_octet_to_all_LL_index_inverse:
#         i = key_i
#         a = map_octet_to_all_LL_index_inverse[key_i] - 1
#         for key_j in map_octet_to_all_LL_index_inverse:
#             j = key_j
#             b = map_octet_to_all_LL_index_inverse[key_j] - 1
#             Hint_16[i][j] = hamtx[a][b]
#     return np.array(Hint_16)
#
#
# ##########################################################################################################
# # regularization (self energy) U dependent
#
# def delta_e_kp(n, nu0, nu1, num2, nu2, eigenvectorp):
#     pzm = (nu1 - 1 / 2) * Xskp(n, 1, 1, n, eigenvectorp) + (nu0 - 1 / 2) * Xskp(n, 0, 0, n, eigenvectorp)
#
#     m2and2 = nu2 * Xskp(n, 2, 2, n, eigenvectorp) - (1 - num2) * Xskp(n, -2, -2, n, eigenvectorp)
#
#     F = (Xskp(n, 2, 2, n, eigenvectorp) - Xskp(n, -2, -2, n, eigenvectorp))
#
#     #     F0=Xskp(0, 2, 2, 0, eigenvectorp) - Xskp(0, -2, -2, 0, eigenvectorp)
#
#     #     F=F-F0
#
#     res = (F / 2 - pzm - m2and2) * k
#     return res
#
#
# def delta_e_km(n, nu0, nu1, num2, nu2, eigenvectorm):
#     pzm = (nu1 - 1 / 2) * Xskm(n, 1, 1, n, eigenvectorm) + (nu0 - 1 / 2) * Xskm(n, 0, 0, n, eigenvectorm)
#
#     m2and2 = nu2 * Xskm(n, 2, 2, n, eigenvectorm) - (1 - num2) * Xskm(n, -2, -2, n, eigenvectorm)
#
#     F = (Xskm(n, 2, 2, n, eigenvectorm) - Xskm(n, -2, -2, n, eigenvectorm))
#
#     #     F0=Xskm(0, 2, 2, 0, eigenvectorm) - Xskm(0, -2, -2, 0, eigenvectorm)
#
#     #     F=F-F0
#
#     res = (F / 2 - pzm - m2and2) * k
#     return res
#
#
# def delta_e_regmatrix(rho0const, eigenvectorp, eigenvectorm):
#     # print('here3 ', np.diag(rho0const))
#     nu0kp, nu1kp, num2kp, nu2kp, nu0km, nu1km, num2km, nu2km = tuple(np.diag(rho0const)[:8])
#     # print('here4')
#     regmatrixUmspindown = [delta_e_kp(n, nu0kp, nu1kp, num2kp, nu2kp, eigenvectorp) for n in setH] + [delta_e_km(n, nu0km, nu1km, num2km, nu2km, eigenvectorm) for n in setH]
#
#     nu0kp, nu1kp, num2kp, nu2kp, nu0km, nu1km, num2km, nu2km = tuple(np.diag(rho0const)[8:16])
#     regmatrixUmspinup = [delta_e_kp(n, nu0kp, nu1kp, num2kp, nu2kp, eigenvectorp) for n in setH] + [delta_e_km(n, nu0km, nu1km, num2km, nu2km, eigenvectorm) for n in setH]
#     return np.diag(np.array(regmatrixUmspindown + regmatrixUmspinup))
#
#
# ##########################################################################################################
#
# def exciton_j_to_n_km(n, j, eigenvectorm):
#     A_nj = Xskm(n, n, j, j, eigenvectorm) * k * alpha_k
#     nu_j = 1
#     nu_n = 0
#     return -(nu_j - nu_n) * A_nj
#
#
# def exciton_j_to_n_kp(n, j, eigenvectorp):
#     A_nj = Xskp(n, n, j, j, eigenvectorp) * k * alpha_k
#     nu_j = 1
#     nu_n = 0
#     return -(nu_j - nu_n) * A_nj
#
#
# #########################################################################################################
# def asymmetric_h(rho, tau, u):
#     first = u * np.trace(tau @ rho) * (tau @ rho)
#     second = - u * tau @ rho @ tau @ rho
#     return first + second
#
#
# taux = tau_func([[0, 1], [1, 0]])
# tauy = tau_func([[0, -1], [1, 0]]) * 1j
# tauz = tau_func([[1, 0], [0, -1]])
#
#
# def H_asym(rho):
#     H_asym_tmp = asymmetric_h(rho, taux, uperp) + asymmetric_h(rho, tauy, uperp) + asymmetric_h(rho, tauz, uz)
#     return H_asym_tmp
#
#
# ##########################################################################################################
#
#
# def loopU0(u):
#     if u >= 0:
#         rho0 = rho0constUp  # rhodiagUp
#     else:
#         rho0 = rho0constUm  # rhodiagUm
#     print('running asymmetric calcs for  nu=%(nu)i u=%(u).2fmeV' % {'u':(u * 1e3),'nu':nu})
#     rho = rho0
#
#     eigenvaluep2, eigenvectorp2 = eigen(hAp(u))[0][1:3], eigen(hAp(u))[1][1:3]
#     eigenvectorp2 = nonedimmerp(eigenvectorp2)
#     # print(eigenvectorp2)
#
#     eigenvaluem2, eigenvectorm2 = eigen(hAp(-u))[0][1:3], eigen(hAp(-u))[1][1:3]
#     eigenvectorm2 = nonedimmerm(eigenvectorm2)
#
#     # print(eigenvectorm2)
#
#     eigenvaluep2 = eigen(hAp(u))[0][1:3]
#     eigenvaluem2 = eigen(hAp(-u))[0][1:3]
#     # print(eigenvaluep2,eigenvaluem2)
#
#     # Delta_ab = dab
#     eigenvaluep1 = eigen(hBp(u))[0][3]
#     eigenvaluem1 = eigen(hBp(-u))[0][3]
#     # print(eigenvaluep1,eigenvaluem1)
#
#     eigenvaluep0 = eigen(hCp(u))[0][2]
#     eigenvaluem0 = eigen(hCp(-u))[0][2]
#     # print(eigenvaluep0,eigenvaluem0)
#
#     eigenvaluep = [eigenvaluep0] + [eigenvaluep1] + (eigenvaluep2).tolist()
#     eigenvaluem = [eigenvaluem0] + [eigenvaluem1] + (eigenvaluem2).tolist()
#     # bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
#     h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)  # follows notation order
#
#     eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
#     eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])
#
#     eigenvectorp_none_interact = eigenvectorp
#     eigenvectorm_none_interact = eigenvectorm
#
#     regmatrix = delta_e_regmatrix(rho0, eigenvectorp, eigenvectorm) * alpha_reg
#
#     it = 1;
#     while it < itmax:
#         Hint_longrange = k * alpha_H_oct_int * Hint_oct(rho)
#         H_asym = asymmetric_h(taux, rho, uperp) + asymmetric_h(tauy, rho, uperp) + asymmetric_h(tauz, rho, uz)
#         H = h0 + Hint_longrange + H_asym + mZm # + regmatrix
#         eigenvalue_loop, eigenvector_loop = eigen(H);
#         rhotemp = rho
#         rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(nu));
#         rho = (1 - alpha_rho) * rho + alpha_rho * rhotemp
#         it += 1
#
#     # eigenvalue, eigenvector = npla.eig(H) # follows notation (bands list) order
#     eigenvalue, eigenvector = eigen(H)
#     # eigenvector = np.transpose(eigenvector)# follows notation (bands list) order
#     # eigenvalue_h0, eigenvector = npla.eig(h0+mZm)
#     # eigenvalue_H_asym, eigenvector = npla.eig(H_asym)
#
#     dict_quantities_u = {'u': u * 1e3,
#                          'eigenvalue': 1e3 * np.real(eigenvalue),
#                          'eigenvector': eigenvector,
#                          # 'Et': 1e3 * Et,
#                          # 'h0': df_round(1e3 * h0),
#                          'rhoU': df_round(rho),
#                          # 'Eh_deltaU': 1e3 * k * Eh * deltatb,
#                          # 'Hint': df_round(1e3 * Hint),
#                          'regmatrix': 1e3 * regmatrix
#                          }
#     exciton = np.array([exciton_j_to_n_km(-2, 1, eigenvectorm_none_interact), exciton_j_to_n_kp(-2, 1, eigenvectorp_none_interact),
#                         exciton_j_to_n_km(1, 2, eigenvectorm_none_interact), exciton_j_to_n_kp(1, 2, eigenvectorp_none_interact)])
#
#     dict_quantities_u['exciton_energy'] = 1e3 * exciton
#     # print(sorted(np.real(eigenvalue_h0)*1e3))
#     # print(sorted(np.real(eigenvalue_H_asym)*1e3))
#     # print(sorted(np.real(eigenvalue)*1e3))
#
#     # return [u * 1e3]+ (np.real(eigenvalue) * 1e3).tolist()
#     return dict_quantities_u

# print(rho0constUp)
# print(bands_oct)
# print(bands_LL2)
# base_oct = ['0m-', '0p-', '1m-', '1p-', '0m+', '0p+', '1m+', '1p+']
# bands_oct = ['0m-', '0p-', '1m-', '1p-', '0m+', '0p+', '1m+', '1p+']
# bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
# bands_oct = [ 'LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL1_Km_Sdown', 'LL1_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sup', 'LL1_Kp_Sup']
# bands = bands

map_octet_to_all_LL_index = {}
for band in bands_oct:
    #     print(base_oct.index(base)+1,bands.index(base))
    map_octet_to_all_LL_index[bands_oct.index(band) + 1] = bands.index(band)

map_octet_to_all_LL_index_inverse = {}
for key in map_octet_to_all_LL_index:
    map_octet_to_all_LL_index_inverse[map_octet_to_all_LL_index[key]] = key


# def rho(i, j):
#     return rhotmp[map_octet_to_all_LL_index[i]][map_octet_to_all_LL_index[j]]

# map_octet_to_all_LL_index_inverse
def Hint_oct(rhotmp):
    # global rhotmp
    def rho(i, j):
        return rhotmp[map_octet_to_all_LL_index[i]][map_octet_to_all_LL_index[j]]

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

    for key_i in map_octet_to_all_LL_index_inverse:
        i = key_i
        a = map_octet_to_all_LL_index_inverse[key_i] - 1
        for key_j in map_octet_to_all_LL_index_inverse:
            j = key_j
            b = map_octet_to_all_LL_index_inverse[key_j] - 1
            Hint_16[i][j] = hamtx[a][b]

    return np.array(Hint_16)


########################################################################
# ##########################################################################################################
# # regularization (self energy) U dependent
#
# def delta_e_kp(n, nu0, nu1, num2, nu2, eigenvectorp):
#     pzm = (nu1 - 1 / 2) * Xskp(n, 1, 1, n, eigenvectorp) + (nu0 - 1 / 2) * Xskp(n, 0, 0, n, eigenvectorp)
#
#     m2and2 = nu2 * Xskp(n, 2, 2, n, eigenvectorp) - (1 - num2) * Xskp(n, -2, -2, n, eigenvectorp)
#
#     F = (Xskp(n, 2, 2, n, eigenvectorp) - Xskp(n, -2, -2, n, eigenvectorp))
#
#     #     F0=Xskp(0, 2, 2, 0, eigenvectorp) - Xskp(0, -2, -2, 0, eigenvectorp)
#
#     #     F=F-F0
#
#     res = (F / 2 - pzm - m2and2) * k
#     return res
#
#
# def delta_e_km(n, nu0, nu1, num2, nu2, eigenvectorm):
#     pzm = (nu1 - 1 / 2) * Xskm(n, 1, 1, n, eigenvectorm) + (nu0 - 1 / 2) * Xskm(n, 0, 0, n, eigenvectorm)
#
#     m2and2 = nu2 * Xskm(n, 2, 2, n, eigenvectorm) - (1 - num2) * Xskm(n, -2, -2, n, eigenvectorm)
#
#     F = (Xskm(n, 2, 2, n, eigenvectorm) - Xskm(n, -2, -2, n, eigenvectorm))
#
#     #     F0=Xskm(0, 2, 2, 0, eigenvectorm) - Xskm(0, -2, -2, 0, eigenvectorm)
#
#     #     F=F-F0
#
#     res = (F / 2 - pzm - m2and2) * k
#     return res
#
#
# def delta_e_regmatrix(rho0const, eigenvectorp, eigenvectorm):
#     # print('here3 ', np.diag(rho0const))
#     nu0kp, nu1kp, num2kp, nu2kp, nu0km, nu1km, num2km, nu2km = tuple(np.diag(rho0const)[:8])
#     # print('here4')
#     regmatrixUmspindown = [delta_e_kp(n, nu0kp, nu1kp, num2kp, nu2kp, eigenvectorp) for n in setH] + [delta_e_km(n, nu0km, nu1km, num2km, nu2km, eigenvectorm) for n in setH]
#
#     nu0kp, nu1kp, num2kp, nu2kp, nu0km, nu1km, num2km, nu2km = tuple(np.diag(rho0const)[8:16])
#     regmatrixUmspinup = [delta_e_kp(n, nu0kp, nu1kp, num2kp, nu2kp, eigenvectorp) for n in setH] + [delta_e_km(n, nu0km, nu1km, num2km, nu2km, eigenvectorm) for n in setH]
#     return np.diag(np.array(regmatrixUmspindown + regmatrixUmspinup))
#
#
# ##########################################################################################################
#
# def exciton_j_to_n_km(n, j, eigenvectorm):
#     A_nj = Xskm(n, n, j, j, eigenvectorm) * k * alpha_k
#     nu_j = 1
#     nu_n = 0
#     return -(nu_j - nu_n) * A_nj
#
#
# def exciton_j_to_n_kp(n, j, eigenvectorp):
#     A_nj = Xskp(n, n, j, j, eigenvectorp) * k * alpha_k
#     nu_j = 1
#     nu_n = 0
#     return -(nu_j - nu_n) * A_nj
#
#
# #########################################################################################################
#########################################################################################################
# def tau_func(tau):
#     tau = np.array(tau)
#     temp = np.zeros((8, 8), dtype='float_')  # ,dtype = 'complex_')
#
#     # for i in range(4):
#     for i in range(2):
#         temp[i, i] = tau[0, 0]
#         temp[i, i + 4] = tau[0, 1]
#         temp[i + 4, i] = tau[1, 0]
#         temp[i + 4, i + 4] = tau[1, 1]
#
#     return np.kron(np.eye(2, dtype=int), temp)

taux = tau_func([[0, 1], [1, 0]])
tauy = tau_func([[0, -1], [1, 0]]) * 1j
tauz = tau_func([[1, 0], [0, -1]])

def asymmetric_h(tau, rho, u):
    first = u * np.trace(tau @ rho) * (tau @ rho)
    second = - u * tau @ rho @ tau @ rho
    return first + second


# def H_asym(rho):
#     H_asym_tmp = asymmetric_h(rho, taux, uperp) + asymmetric_h(rho, tauy, uperp) + asymmetric_h(rho, tauz, uz)
#     return H_asym_tmp


##########################################################################################################


def loopU0(u):
    # global rhotemp
    if u >= 0:
        rho0 = rho0constUp  # rhodiagUp
    else:
        rho0 = rho0constUm  # rhodiagUm
    print('running asymmetric calcs for  nu=%(nu)i u=%(u).2fmeV' % {'u': (u * 1e3), 'nu': nu})
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

    # print(h0)

    # regmatrix = delta_e_regmatrix(rho0, eigenvectorp, eigenvectorm) * alpha_reg
    # print(np.diag(rho.real))
    # print('here')
    # print(k)
    it = 1
    while it < itmax:
        Hint_longrange = k * alpha_H_oct_int * Hint_oct(rho)
        H_asym = asymmetric_h(taux, rho, uperp) + asymmetric_h(tauy, rho, uperp) + asymmetric_h(tauz, rho, uz)
        H = h0 + Hint_longrange + H_asym + mZm
        eigenvalue_loop, eigenvector_loop = eigen(H)
        rhotemp = rho
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(number_occupied_bands))
        rho = (1 - alpha_rho) * rho + alpha_rho * rhotemp
        it += 1

    # print(Hint_longrange)
    # eigenvalue, eigenvector = npla.eig(H) # follows notation (bands list) order
    eigenvalue, eigenvector = eigen(H)
    # print(eigenvalue*1e3)
    # eigenvector = np.transpose(eigenvector)# follows notation (bands list) order
    # eigenvalue_h0, eigenvector = npla.eig(h0+mZm)
    # eigenvalue_H_asym, eigenvector = npla.eig(H_asym)
    #
    #
    # print(sorted(np.real(eigenvalue_h0)*1e3))
    # print(sorted(np.real(eigenvalue_H_asym)*1e3))
    # print(sorted(np.real(eigenvalue)*1e3))

    # return [u * 1e3]+ (np.real(eigenvalue) * 1e3).tolist()
    dict_quantities_u = {'u': u * 1e3,
                         'eigenvalue': 1e3 * np.real(eigenvalue),
                         'eigenvector': eigenvector,
                         # 'Et': 1e3 * Et,
                         # 'h0': df_round(1e3 * h0),
                         'rhoU': df_round(rho),
                         # 'Eh_deltaU': 1e3 * k * Eh * deltatb,
                         # 'Hint': df_round(1e3 * Hint),
                         # 'regmatrix': 1e3 * regmatrix
                         }
    # exciton = np.array([exciton_j_to_n_km(-2, 1, eigenvectorm_none_interact), exciton_j_to_n_kp(-2, 1, eigenvectorp_none_interact),
    #                     exciton_j_to_n_km(1, 2, eigenvectorm_none_interact), exciton_j_to_n_kp(1, 2, eigenvectorp_none_interact)])

    # dict_quantities_u['exciton_energy'] = 1e3 * exciton
    # print(sorted(np.real(eigenvalue_h0)*1e3))
    # print(sorted(np.real(eigenvalue_H_asym)*1e3))
    # print(sorted(np.real(eigenvalue)*1e3))

    # return [u * 1e3]+ (np.real(eigenvalue) * 1e3).tolist()
    return dict_quantities_u