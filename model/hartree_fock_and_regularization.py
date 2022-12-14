from config import setH, indexes_octet_on_bands, index_octet_on_bands_oct  # itmax_full_range,
from input.parameters import *
from model.densities import density_by_model_regime
# from model.density_test import rho0constUp, rho0constUm
from model.exchange_integrals import Xskm, Xskp
from model.hamiltonians import h0p, h0p2, h0m, h0m2, full_hp, full_hm, full_hpm, full_hmp, idp, idps, idm, idms, mZm, hAp, hBp, hCp, tau_func
from utils import eigen, df_round, nonedimmerp, nonedimmerm, remove_small_imag, check_if_complex, occupation_band_regularization

# if model_regime == 'no_LL2_mixing_and_asym':
#     print('executing hartree_fock_and_regularization to return LL2 and LLm2 for no_LL2_mixing_and_asym regime')
# elif model_regime == 'full_range':
#     print('executing hartree_fock_and_regularization in full_range regime')

model_regime = 'full_range'


# print('here_hartree_fock_only_calcs')
# # import time
# # time.sleep(.2)

def hp(n, nprime, s1, s2):
    """ interacting hamiltonian for valley k^p """
    return full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho)


def hm(n, nprime, s1, s2):
    """ interacting hamiltonian for valley k^m """
    return full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho)


def hpm(n, nprime, s1, s2):
    """ interacting hamiltonian for valley k^p x k^m """
    if potential_asym_layers:
        return full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)
    else:
        return 0


def hmp(n, nprime, s1, s2):
    """ interacting hamiltonian for valley k^m x k^p """
    if potential_asym_layers:
        return full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)
    else:
        return 0


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
    state_occupancy = {0: occupation_band_regularization(ll0),
                       1: occupation_band_regularization(ll1),
                       -2: occupation_band_regularization(llm2),
                       2: occupation_band_regularization(ll2)}
    # if nu==0: state_occupancy = {0: 1, 1: 0, -2: 1, 2: 0} #occupancy as of the LL porlarized state
    res = sum([(1 / 2 - state_occupancy[m]) * Xskp(n, m, m, n, eigenvectorp) for m in setH]) * k  # boxed equation on regularization.pdf

    return res


def delta_e_km(n, ll0, ll1, llm2, ll2, eigenvectorm):
    # Shizuya, PRB 2012 eq. 24 / Shizuya, PRB 2020 eq. 32 / and my regularization notes
    state_occupancy = {0: occupation_band_regularization(ll0),
                       1: occupation_band_regularization(ll1),
                       -2: occupation_band_regularization(llm2),
                       2: occupation_band_regularization(ll2)}
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

def exciton_j_to_n_km(n, j, eigenvectorm):
    A_nj = Xskm(n, n, j, j, eigenvectorm) * k
    nu_j = 1
    nu_n = 0
    return -(nu_j - nu_n) * A_nj


def exciton_j_to_n_kp(n, j, eigenvectorp):
    A_nj = Xskp(n, n, j, j, eigenvectorp) * k
    nu_j = 1
    nu_n = 0
    return -(nu_j - nu_n) * A_nj


########################################################################################
taux = tau_func([[0, 1], [1, 0]])
tauy = tau_func([[0, -1], [1, 0]]) * 1j
tauz = tau_func([[1, 0], [0, -1]])


def asymmetric_h(tau, rho, u):
    first = u * np.trace(tau @ rho) * (tau @ rho) / 2
    second = - u * tau @ rho @ tau @ rho / 2
    return first + second

if same_rhoRandom:
    rho0constUm = np.array(density_by_model_regime(model_regime)['rho0constUm'])
    rho0constUp = np.array(density_by_model_regime(model_regime)['rho0constUp'])
    rho0const_small_u = np.array(density_by_model_regime(model_regime)['rho0const_small_u'])


def loopU(u):
    global eigenvectorp, eigenvectorm, deltatb, rho
    if same_rhoRandom:
        if u <= -u_critical:
            # print('negative u')
            rho0 = rho0constUm
        elif (u > -u_critical) and (u < u_critical):
            # print('small u')
            rho0 = rho0const_small_u

        elif u >= u_critical:
            # print('positive u')
            rho0 = rho0constUp
    else:
        if u <= -u_critical:
            # print('negative u')
            rho0 = density_by_model_regime(model_regime)['rho0constUm']
        elif (u > -u_critical) and (u < u_critical):
            # print('small u')
            rho0 = density_by_model_regime(model_regime)['rho0const_small_u']
        elif u >= u_critical:
            # print('positive u')
            rho0 = density_by_model_regime(model_regime)['rho0constUp']

    print('running hartree_fock_and_regularization with nu=%(nu)i u=%(u).2fmeV ' % {'u': (u * 1e3), 'nu': nu})
    rho = rho0

    ################### warping #############################################################
    if projected_four_band_H0:
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
    if effective_H0:
        eigenvaluep2, eigenvectorp2 = eigen(h0p2(u))
        eigenvaluem2, eigenvectorm2 = eigen(h0m2(u))
        eigenvaluep = [h0p(u)[0][0]] + [h0p(u)[1][1]] + eigenvaluep2.tolist()
        eigenvaluem = [h0m(u)[0][0]] + [h0m(u)[1][1]] + eigenvaluem2.tolist()

    # ######################################
    h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)  # none interacting matrix

    # eigenvectorp = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()]))
    # eigenvectorm = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()]))

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    # eigenvectorp_none_interact = eigenvectorp
    # eigenvectorm_none_interact = eigenvectorm
    ###### regularization (self energy) U dependent
    # print('here1')

    # print('here2')
    ######

    it = 1
    while it < itmax_full_range:
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

        H_asym = asymmetric_h(taux, rho, uperp) + asymmetric_h(tauy, rho, uperp) + asymmetric_h(tauz, rho, uz)
        Hint = k * alpha_int_H * np.vstack((np.hstack((Hintup, Hintupdown)), np.hstack((Hintdownup, Hintdown)))) + H_asym * apha_H_asym
        # regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg
        H = Hint + h0 + mZm  # + regmatrix  # np.add(Hint, h0)
        eigenvalue_loop, eigenvector_loop = eigen(H)
        # rhotemp = rho
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(number_occupied_bands))
        # rho = (1 - alpha_rho) * rho + alpha_rho * rhotemp

        # regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg # we should not update regmatrix here, won't affect nu=4

        it += 1
    ###### regularization (self energy) U dependent
    regmatrix = delta_e_regmatrix(rho, eigenvectorp, eigenvectorm) * alpha_reg
    # constant_matrix = np.diag([])
    ######
    H = H + regmatrix  # + constant_matrix
    eigenvalue, eigenvector = eigen(H)
    eigenvector_octet = eigenvector[4:12, index_octet_on_bands_oct]
    # print(eigenvector_octet.shape)
    if check_if_complex(eigenvalue, u, nu):
        # print('here')
        # print('here0', check_if_complex(energy_u, ind, nu))
        eigenvalue = np.real(eigenvalue)
        # print('here0', check_if_complex(energy_u, ind, nu))
        # print('here2')

    eigenvector_octet_norms = np.linalg.norm(eigenvector_octet, axis=1)
    rho_diag = remove_small_imag(np.diag(rho))
    rho_diag_octet = remove_small_imag(np.diag(rho)[indexes_octet_on_bands])
    trace_rho_diag_octet = round(sum(rho_diag_octet), 3)
    sum_off_diag_rho_octet = round(np.sum(rho) - trace_rho_diag_octet - 4, 3)
    rho_diag_octet_real_part = np.real(np.diag(rho)[indexes_octet_on_bands])
    trace_rho_diag_octet_real_part = round(sum(rho_diag_octet_real_part), 3)
    # eigenvalue, eigenvector = npla.eig(H)
    # idxfunc = np.argsort(eigenvalue)
    #
    # eigenvalue = eigenvalue[idxfunc]
    # eigenvector = eigenvector[:, idxfunc]
    # if check_if_complex(eigenvalue, u, nu):
    #     # print('here')
    #     # print('here0', check_if_complex(energy_u, ind, nu))
    #     eigenvalue = np.real(eigenvalue)
    ehf = - sum([Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in setH for nprime in setH] +
                [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in setH for nprime in setH] +
                [Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in setH for nprime in setH] +
                [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in setH for nprime in setH]
                # [2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in setH for nprime in setH] +
                # [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in setH for nprime in setH]
                )

    Et = sum([eigenvalue[i] for i in range(number_occupied_bands)]) + ehf

    # ehf = - sum([Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in setH for nprime in setH] +
    #             [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in setH for nprime in setH] +
    #             [Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in setH for nprime in setH] +
    #             [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in setH for nprime in setH] +
    #             [2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in setH for nprime in setH] +
    #             [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in setH for nprime in setH]
    #             )
    # ehf = - sum([Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in [0,1] for nprime in setH] +
    #             [Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in [-2, 2] for nprime in setH] +
    #             [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in [0,1] for nprime in setH] +
    #             [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in [-2, 2] for nprime in setH] +
    #             [Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in [0,1] for nprime in setH] +
    #             [Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in [-2, 2] for nprime in setH] +
    #             [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in [0,1] for nprime in setH] +
    #             [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in [-2, 2] for nprime in setH] +
    #             [2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in [0,1] for nprime in setH] +
    #             [2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in [-2, 2] for nprime in setH] +
    #             [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in [0,1] for nprime in setH] +
    #             [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in [-2, 2] for nprime in setH]
    #             )
    # Et = sum(np.sort(np.diag(h0+mZm+regmatrix))[:number_occupied_bands]) - ehf

    dict_quantities_u = {'u': u * 1e3,
                         'eigenvalue': 1e3 * remove_small_imag(eigenvalue),  # np.real due to numerical fluctuations
                         'eigenvector': eigenvector,
                         'eigenvector_octet': eigenvector_octet,
                         'eigenvector_octet_norms': eigenvector_octet_norms,
                         'rho(density)': df_round(rho),
                         'rho_diag': rho_diag,
                         'rho_diag_octet': rho_diag_octet,
                         'trace_rho_diag_octet': trace_rho_diag_octet,
                         'sum_off_diag_rho_octet': sum_off_diag_rho_octet,
                         'rho_diag_octet_real_part': rho_diag_octet_real_part,
                         'trace_rho_diag_octet_real_part': trace_rho_diag_octet_real_part,
                         'Et': 1e3 * remove_small_imag(Et),
                         'h0': df_round(1e3 * h0),
                         'rhoU': df_round(rho),
                         'Eh_deltaU': 1e3 * k * Eh * deltatb,
                         'Hint': df_round(1e3 * Hint),
                         'regmatrix': 1e3 * np.diag(regmatrix),
                         'regmatrix_oct': 1e3 * np.diag(regmatrix)[index_octet_on_bands_oct]
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

    exciton = np.array([exciton_j_to_n_km(-2, 1, eigenvectorm), exciton_j_to_n_kp(-2, 1, eigenvectorp),
                        exciton_j_to_n_km(1, 2, eigenvectorm), exciton_j_to_n_kp(1, 2, eigenvectorp),
                        exciton_j_to_n_km(2, 1, eigenvectorm), exciton_j_to_n_kp(2, 1, eigenvectorp),
                        exciton_j_to_n_km(1, -2, eigenvectorm), exciton_j_to_n_kp(1, -2, eigenvectorp)])

    dict_quantities_u['exciton_energy'] = 1e3 * exciton
    # dict_quantities_u['exciton'] = exciton

    # for t in allowed_transitions:
    #     transition_energy()

    return dict_quantities_u
