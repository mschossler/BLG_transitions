from config import bands, base_octet, indexes_octet_on_bands, index_octet_on_bands_oct, setH  # , itmax_asymmetric_calcs, alpha_rho
from input.parameters import *
from model.densities import density_by_model_regime
from model.exchange_integrals import Xzs, Xzd, Xos, Xod, Xfs, Xfd, Xsts, Xstd, Xskm, Xskp
# from model.exchange_integrals import Xskm, Xskp
from model.hamiltonians import mZm, hAp, hBp, hCp, tau_func, idp, idps, idm, idms  # , asymmetric_h, taux, tauy, tauz
from utils import eigen, nonedimmerp, nonedimmerm, df_round, remove_small_imag, check_if_complex, occupation_band

model_regime = 'no_LL2_mixing_and_asym'

# print('executing hartree_fock_with_asymmetric_interactions for no_LL2_mixing_and_asym regime')

# print('here_asymmetric_calcs')
# import time
# time.sleep(.2)
map_octet_to_all_LL_index = {}
for band in base_octet:
    #     print(base_oct.index(base)+1,bands.index(base))
    map_octet_to_all_LL_index[base_octet.index(band) + 1] = bands.index(band)

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
##########################################################################################################
# regularization (self energy) U dependent
# def occupation_band(x):
#     if abs(x) < 0.4:
#     # if round(x)==0:
#         return abs(x)
#     # if x == 0:
#     # return abs(x)
#     else:
#         return 1


def delta_e_kp(n, ll0, ll1, llm2, ll2, eigenvectorp):
    # Shizuya, PRB 2012 eq. 24 / Shizuya, PRB 2020 eq. 32  / and my notes
    state_occupancy = {0: occupation_band(ll0),
                       1: occupation_band(ll1),
                       -2: occupation_band(llm2),
                       2: occupation_band(ll2)}
    # if nu==0: state_occupancy = {0: 1, 1: 0, -2: 1, 2: 0} #occupancy as of the LL porlarized state
    res = sum([(1 / 2 - state_occupancy[m]) * Xskp(n, m, m, n, eigenvectorp) for m in setH]) * k  # boxed equation on regularization.pdf

    return res


def delta_e_km(n, ll0, ll1, llm2, ll2, eigenvectorm):
    # Shizuya, PRB 2012 eq. 24 / Shizuya, PRB 2020 eq. 32 / and my regularization notes
    state_occupancy = {0: occupation_band(ll0),
                       1: occupation_band(ll1),
                       -2: occupation_band(llm2),
                       2: occupation_band(ll2)}
    # if nu==0: state_occupancy = {0: 1, 1: 0, -2: 1, 2: 0} #occupancy as of the LL porlarized state
    res = sum([(1 / 2 - state_occupancy[m]) * Xskm(n, m, m, n, eigenvectorm) for m in setH]) * k  # boxed equation on regularization.pdf

    return res


def delta_e_regmatrix(rho0const, eigenvectorp, eigenvectorm):
    # Shizuya, PRB 2012 eq. 24 / Shizuya, PRB 2020 eq. 32  / and my notes
    # print('here3 ', np.diag(rho0const))
    ll0kp, ll1kp, llm2kp, ll2kp, ll0km, ll1km, llm2km, ll2km = tuple(np.diag(rho0const)[:8])
    # print('here4')
    regmatrixUmspindown = [delta_e_kp(n, ll0kp, ll1kp, llm2kp, ll2kp, eigenvectorp) for n in setH] + [delta_e_km(n, ll0km, ll1km, llm2km, ll2km, eigenvectorm) for n in setH]

    ll0kp, ll1kp, llm2kp, ll2kp, ll0km, ll1km, llm2km, ll2km = tuple(np.diag(rho0const)[8:16])
    regmatrixUmspinup = [delta_e_kp(n, ll0kp, ll1kp, llm2kp, ll2kp, eigenvectorp) for n in setH] + [delta_e_km(n, ll0km, ll1km, llm2km, ll2km, eigenvectorm) for n in setH]
    return np.diag(np.array(regmatrixUmspindown + regmatrixUmspinup))


##########################################################################################################
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
    first = u * np.trace(tau @ rho) * (tau @ rho) / 2
    second = - u * tau @ rho @ tau @ rho / 2
    return first + second


# def H_asym(rho):
#     H_asym_tmp = asymmetric_h(rho, taux, uperp) + asymmetric_h(rho, tauy, uperp) + asymmetric_h(rho, tauz, uz)
#     return H_asym_tmp


##########################################################################################################

if same_rhoRandom:
    rho0const_small_u = np.array(density_by_model_regime(model_regime)['rho0const_small_u'])


def loopU0(u):
    global rhotemp
    # if u >= 0:
    #     rho0 = density_by_model_regime(model_regime)['rho0constUp']  # rhodiagUp
    # else:
    #     rho0 = density_by_model_regime(model_regime)['rho0constUm']  # rhodiagUm
    rho0 = np.array(density_by_model_regime(model_regime)['rho0const_small_u'])
    print('running hartree_fock_with_asymmetric_interactions for  nu=%(nu)i u=%(u).2fmeV ' % {'u': (u * 1e3), 'nu': nu})

    if same_rhoRandom:
        rho0 = rho0const_small_u

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

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    # regmatrix = delta_e_regmatrix(rho0, eigenvectorp, eigenvectorm) * alpha_reg_asym_calcs
    # print(np.diag(rho.real))
    # print('here')
    # print(k)
    it = 1
    while it < itmax_asymmetric_calcs:
        Hint_longrange = k * alpha_H_oct_int * Hint_oct(rho)
        H_asym = asymmetric_h(taux, rho, uperp) + asymmetric_h(tauy, rho, uperp) + asymmetric_h(tauz, rho, uz)
        Hint = Hint_longrange + H_asym * apha_H_asym
        # regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg
        H = h0 + Hint + mZm  # + regmatrix
        eigenvalue_loop, eigenvector_loop = eigen(H)
        rhotemp = rho
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(number_occupied_bands))
        rho = (1 - alpha_rho) * rho + alpha_rho * rhotemp
        # regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg_asym_calcs # we should not update regmatrix here, won't affect nu=4 besides a small imag part
        it += 1

    # print(Hint_longrange)
    # eigenvalue, eigenvector = npla.eig(H) # follows notation (bands list) order
    regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg_asym_calcs
    H = H + regmatrix
    eigenvalue, eigenvector = eigen(H)

    ehf = - sum([Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in setH for nprime in setH] +
                [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in setH for nprime in setH] +
                [Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in setH for nprime in setH] +
                [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in setH for nprime in setH]
                # [2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in setH for nprime in setH] +
                # [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in setH for nprime in setH]
                )
    Et = sum([eigenvalue[i] for i in range(number_occupied_bands)]) + ehf
    eigenvector_octet = eigenvector[4:12, index_octet_on_bands_oct]
    # print(eigenvector_octet.shape)
    if check_if_complex(eigenvalue, u, nu):
        # print('here')
        # print('here0', check_if_complex(energy_u, ind, nu))
        eigenvalue = np.real(eigenvalue)
        # print('here0', check_if_complex(energy_u, ind, nu))
        # print('here2')
    # eigenvalue = remove_small_imag(eigenvalue)
    # eigenvector_octet_norms = [np.linalg.norm(one_eigenvector_octet) for one_eigenvector_octet in  ]
    eigenvector_octet_norms = np.linalg.norm(eigenvector_octet, axis=1)

    rho_diag = np.diag(rho)
    rho_diag_octet = np.diag(rho)[indexes_octet_on_bands]
    trace_rho_diag_octet = round(sum(rho_diag_octet), 3)
    sum_off_diag_rho_octet = round(np.sum(rho) - trace_rho_diag_octet - 4, 3)
    rho_diag_octet_real_part = np.real(np.diag(rho)[indexes_octet_on_bands])
    trace_rho_diag_octet_real_part = round(sum(rho_diag_octet_real_part), 3)
    # print(eigenvalue*1e3)
    # eigenvector = np.transpose(eigenvector)# follows notation (bands list) order
    # eigenvalue_h0, eigenvector = npla.eig(h0+mZm)
    # eigenvalue_H_asym, eigenvector = npla.eig(H_asym)
    #
    #
    # print(sorted(remove_small_imag(eigenvalue_h0)*1e3))
    # print(sorted(remove_small_imag(eigenvalue_H_asym)*1e3))
    # print(sorted(remove_small_imag(eigenvalue)*1e3))

    # return [u * 1e3]+ (remove_small_imag(eigenvalue) * 1e3).tolist()
    dict_quantities_u = {'u': u * 1e3,
                         'eigenvalue': 1e3 * remove_small_imag(eigenvalue),
                         'eigenvector': eigenvector,
                         'eigenvector_octet': eigenvector_octet,
                         'eigenvector_octet_norms': eigenvector_octet_norms,
                         'regmatrix': 1e3 * np.diag(regmatrix),
                         'Et': 1e3 * Et,
                         # 'h0': df_round(1e3 * h0),
                         'rho(density)': df_round(rho),
                         'rho_diag': rho_diag,
                         'rho_diag_octet': rho_diag_octet,
                         'trace_rho_diag_octet': trace_rho_diag_octet,
                         'sum_off_diag_rho_octet': sum_off_diag_rho_octet,
                         'rho_diag_octet_real_part': rho_diag_octet_real_part,
                         'trace_rho_diag_octet_real_part': trace_rho_diag_octet_real_part
                         # 'Eh_deltaU': 1e3 * k * Eh * deltatb,
                         # 'Hint': df_round(1e3 * Hint),
                         # 'regmatrix': 1e3 * regmatrix
                         }
    # exciton = np.array([exciton_j_to_n_km(-2, 1, eigenvectorm_none_interact), exciton_j_to_n_kp(-2, 1, eigenvectorp_none_interact),
    #                     exciton_j_to_n_km(1, 2, eigenvectorm_none_interact), exciton_j_to_n_kp(1, 2, eigenvectorp_none_interact)])

    # dict_quantities_u['exciton_energy'] = 1e3 * exciton
    # print(sorted(remove_small_imag(eigenvalue_h0)*1e3))
    # print(sorted(remove_small_imag(eigenvalue_H_asym)*1e3))
    # print(sorted(remove_small_imag(eigenvalue)*1e3))

    # return [u * 1e3]+ (remove_small_imag(eigenvalue) * 1e3).tolist()
    return dict_quantities_u
