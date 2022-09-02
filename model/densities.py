import numpy as np
import pandas as pd

input_dir = 'input/'

if __name__ == "__main__":
    # setting path
    import sys

    sys.path.append('../')
    input_dir = '../input/'

from input.parameters import nu, number_occupied_bands, model_regime
from config import bands, bands_oct, alpha_rand_full_range, alpha_rand_asymmetric_calcs
from utils import eigen

# filling_order_Upositive = ['-2p-', '-2p+', '-2m-', '-2m+', '0m-', '0p-', '0m+', '0p+', '1m-', '1p-', '1m+', '1p+', '2p-', '2p+']  # spin_polarz -> 0LL_polarz -> 1LL_spin_polarz
# filling_order_Unegative = ['-2m-', '-2m+', '-2p-', '-2p+', '0p-', '0m-', '0p+', '0m+', '1p-', '1m-', '1p+', '1m+', '2m-', '2m-']  # spin_polarz -> 0LL_polarz -> 1LL_spin_polarz

filling_order_Upositive = ['LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sdown',
                           'LL1_Kp_Sdown', 'LL1_Km_Sup', 'LL1_Kp_Sup', 'LL2_Kp_Sdown', 'LL2_Kp_Sup']

filling_order_Unegative = ['LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LL0_Kp_Sdown', 'LL0_Km_Sdown', 'LL0_Kp_Sup', 'LL0_Km_Sup', 'LL1_Kp_Sdown',
                           'LL1_Km_Sdown', 'LL1_Kp_Sup', 'LL1_Km_Sup', 'LL2_Km_Sdown', 'LL2_Km_Sup']

# bands_oct = ['LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL1_Km_Sdown', 'LL1_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sup', 'LL1_Kp_Sup']
seed_oct_dict = {-4: (0, 0, 0, 0, 0, 0, 0, 0),
                 -3: (1, 0, 0, 0, 0, 0, 0, 0),
                 -2: (1, 1, 0, 0, 0, 0, 0, 0),
                 -1: (1, 1, 1, 0, 0, 0, 0, 0),
                 0: (1, 1, 1, 1, 0, 0, 0, 0),
                 1: (1, 1, 1, 1, 1, 0, 0, 0),
                 2: (1, 1, 1, 1, 1, 1, 0, 0),
                 3: (1, 1, 1, 1, 1, 1, 1, 0),
                 4: (1, 1, 1, 1, 1, 1, 1, 1)
                 }


def ramdom_16x16_density(numb_of_filled_states):
    if (model_regime == 'full_range') and (nu == 0):
        rhorand8 = pd.read_csv('input/' + 'rho0phbroken.csv', header=None).values.tolist()

        rho0file8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
        zero_rho0file8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

        rhorand16 = rho0file8_zero + zero_rho0file8
    else:
        rhorand8 = np.random.rand(8, 8)
        eigenvalue, eigenvector = eigen(rhorand8)
        # rhorand8 = sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(int(number_occupied_bands / 2)))
        rhorand8 = sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(numb_of_filled_states // 2))

        rhorand8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
        zero_rhorand8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

        rhorand16 = rhorand8_zero + zero_rhorand8

    return rhorand16


def diag_full_regime(u_signal, number_occupied_bands):
    if u_signal >= 0:
        filling_order = filling_order_Upositive
    if u_signal < 0:
        filling_order = filling_order_Unegative
    diag = [(1 if band in filling_order[0:number_occupied_bands] else 0) for band in bands]
    if sum(diag) != number_occupied_bands:
        print('is the filling factor right?:', sum(diag) == number_occupied_bands, 'u_signal:', u_signal)
        exit()

    # if u_signal >= 0:
    #     rho0const = rho0constUp
    # if u_signal < 0:
    #     filling_order = filling_order_Unegative
    # seed_dict = {'rho0constUp': rho0constUp, 'rho0constUm': rho0constUm}
    rhorand16 = ramdom_16x16_density(number_occupied_bands)
    return (1 - alpha_rand_full_range) * np.diag(diag) + alpha_rand_full_range * rhorand16


# for band in bands:
#     diag_Up = [1 if band in filling_order_Upositive[0:number_occupied_bands] else 0]
# rho0constUp = np.diag(diag_Up)

# print(diag_Um)

# for band in bands:
#     diag_Um = [1 if band in filling_order_Upositive[0:number_occupied_bands] else 0]
# rho0constUm = np.diag(diag_Um)
# for state in filling_order_Upositive[0:number_occupied_bands]:
#     print(state)
#
#     print('is the filling factor right (U>0)?: ', sum(diag) == number_occupied_bands)
#     rho0constUp = np.diag(diag)
# if sum(diag_Up) != number_occupied_bands:
#     print('is the filling factor right (U>0)?: ', sum(diag_Up) == number_occupied_bands)
#     exit()

# if nu == 0:
#     parameters_density = {'alpha_rand': alpha_rand,
#                           'dens': dens,
#                           'alpha_state': alpha_state
#                           }
#     print('parameters for seed at nu=0: alpha_rand=%(alpha_rand).3f, dens=%(dens)i, alpha_state=%(alpha_state).2f' % parameters_density)
#
#     # rhorand8 = pd.read_csv(input_dir + 'rho0phbroken.csv', header=None).values.tolist()
#     #
#     # rho0file8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
#     # zero_rho0file8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')
#     #
#     # rho0file16 = rho0file8_zero + zero_rho0file8
#     # rhorand = rho0file16
#
#     rhorand8 = np.random.rand(8, 8)
#     eigenvalue, eigenvector = eigen(rhorand8)
#     rhorand8 = sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(int(number_occupied_bands / 2)))
#
#     rhorand8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
#     zero_rhorand8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')
#
#     rhorand16 = rhorand8_zero + zero_rhorand8
#     rhorand = rhorand16
#
#     if dens == 1:
#         spindown = [('0p-', alpha_state), ('1p-', alpha_state), ('-2p-', 1), ('2p-', 0), ('0m-', 1 - alpha_state), ('1m-', 1 - alpha_state), ('-2m-', 1), ('2m-', 0)]
#         spinup = [('0p+', alpha_state), ('1p+', alpha_state), ('-2p+', 1), ('2p+', 0), ('0m+', 1 - alpha_state), ('1m+', 1 - alpha_state), ('-2m+', 1), ('2m+', 0)]
#
#     if dens == 2:  # nu8density_CAF - 11000011
#         spindown = [('0p-', alpha_state), ('1p-', alpha_state), ('-2p-', 1), ('2p-', 0), ('0m-', 1 - alpha_state), ('1m-', 1 - alpha_state), ('-2m-', 1), ('2m-', 0)]
#         spinup = [('0p+', 1 - alpha_state), ('1p+', 1 - alpha_state), ('-2p+', 1), ('2p+', 0), ('0m+', alpha_state), ('1m+', alpha_state), ('-2m+', 1), ('2m+', 0)]
#
#     if dens == 3:  # alpha_state =1 ferromagnetic state, spin polarized 11110000
#         spindown = [('0p-', alpha_state), ('1p-', alpha_state), ('-2p-', 1), ('2p-', 0), ('0m-', alpha_state), ('1m-', alpha_state), ('-2m-', 1), ('2m-', 0)]
#         spinup = [('0p+', 1 - alpha_state), ('1p+', 1 - alpha_state), ('-2p+', 1), ('2p+', 0), ('0m+', 1 - alpha_state), ('1m+', 1 - alpha_state), ('-2m+', 1), ('2m+', 0)]
#
#     if dens == 4:  # alpha_state =1 LL polarz
#         spindown = [('0p-', alpha_state), ('1p-', 1 - alpha_state), ('-2p-', 1), ('2p-', 0), ('0m-', alpha_state), ('1m-', 1 - alpha_state), ('-2m-', 1), ('2m-', 0)]
#         spinup = [('0p+', alpha_state), ('1p+', 1 - alpha_state), ('-2p+', 1), ('2p+', 0), ('0m+', alpha_state), ('1m+', 1 - alpha_state), ('-2m+', 1), ('2m+', 0)]
#
#     diag_caf = [x[1] for x in spindown + spinup]
#     # print(diag_caf)
#     # spindownUp = [('0p-', 1), ('1p-', 0), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 0), ('-2m-', 1), ('2m-', 0)]
#     # spinupUp = [('0p+', 1), ('1p+', 0), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 0), ('-2m+', 1), ('2m+', 0)]
#     #
#     # spindownUm = [('0p-', 1), ('1p-', 0), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 0), ('-2m-', 1), ('2m-', 0)]
#     # spinupUm = [('0p+', 1), ('1p+', 0), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 0), ('-2m+', 1), ('2m+', 0)]
#     #
#     # diag = [x[1] for x in spindownUp + spinupUp]
#     # if sum(diag) != number_occupied_bands:
#     #     print('is the filling factor right (U>0)?: ', sum(diag) == number_occupied_bands)
#     #     exit()
#     rho0constUp = (1 - alpha_rand) * np.diag(diag_caf) + alpha_rand * rhorand
#     # diag = [x[1] for x in spindownUm + spinupUm]
#     # if sum(diag) != number_occupied_bands:
#     #     print('is the filling factor right (U<0)?:', sum(diag) == number_occupied_bands)
#     #     exit()
#     rho0constUm = (1 - alpha_rand) * np.diag(diag_caf) + alpha_rand * rhorand

# if nu == 4:
#     spindownUp = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
#     spinupUp = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]
#
#     spindownUm = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
#     spinupUm = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]
#
# if nu != 0:
#     diag = [occupation[1] for occupation in spindownUp + spinupUp]
#     if sum(diag) != number_occupied_bands:
#         print('is the filling factor right (U>0)?: ', sum(diag) == number_occupied_bands)
#         exit()
#     rho0constUp = np.diag(diag)
#
#     diag = [occupation[1] for occupation in spindownUm + spinupUm]
#     if sum(diag) != number_occupied_bands:
#         print('is the filling factor right (U<0)?:', sum(diag) == number_occupied_bands)
#         exit()
#     rho0constUm = np.diag(diag)

# print(np.diag(rho0constUp))
# print(np.diag(rho0constUm))

# bands_oct = ['LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL1_Km_Sdown', 'LL1_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sup', 'LL1_Kp_Sup']
# seed_oct = (1, 1, 0, 1, 0, 1, 0, 1)
# print(sum(seed_oct))
# #
# # occupied_octet_states = []
# # for i in range(8):
# #     if seed_oct[i]:
# #         occupied_octet_states.append(bands_oct[i])
# #
# # print(occupied_octet_states)
# occupied_octet_states = [bands_oct[i] for i in range(8) if seed_oct[i]]
# #
# #
# print(bands_oct)
# print(occupied_octet_states)

# print(filling_order_Upositive[-2::1])

def seed_asymmetric_calcs(nu):
    seed_oct = seed_oct_dict[nu]
    number_occupied_bands_local = sum(seed_oct) + 4
    # number_occupied_bands = nu + 8
    occupied_octet_states = [bands_oct[i] for i in range(8) if seed_oct[i]]

    filling_order = filling_order_Upositive[0:4] + occupied_octet_states + filling_order_Upositive[-2::1]

    diag = [(1 if band in filling_order[0:number_occupied_bands_local] else 0) for band in bands]
    if sum(diag) != number_occupied_bands:
        print('Wrong filling factor for nu %i', nu)
        exit()

    rhorand16 = ramdom_16x16_density(number_occupied_bands_local)

    rho0const = (1 - alpha_rand_asymmetric_calcs) * np.diag(diag) + alpha_rand_asymmetric_calcs * rhorand16

    return rho0const


########################################################################################################################################################
# def seed_asymmetric_calcs(seed_oct, u_signal=1):
#     number_occupied_bands_local = sum(seed_oct) + 4
#     # number_occupied_bands = nu + 8
#     occupied_octet_states = [bands_oct[i] for i in range(8) if seed_oct[i]]
#     if u_signal >= 0:
#         filling_order = filling_order_Upositive[0:4] + occupied_octet_states + filling_order_Upositive[-2::1]
#     elif u_signal < 0:
#         filling_order = filling_order_Unegative[0:4] + occupied_octet_states + filling_order_Unegative[-2::1]
#
#     rhorand8 = np.random.rand(8, 8)
#     eigenvalue, eigenvector = eigen(rhorand8)
#     # rhorand8 = sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(int(number_occupied_bands / 2)))
#     rhorand8 = sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(number_occupied_bands_local // 2))
#
#     rhorand8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
#     zero_rhorand8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')
#
#     rhorand16 = rhorand8_zero + zero_rhorand8
#
#     diag = [(1 if band in filling_order[0:number_occupied_bands_local] else 0) for band in bands]
#     if sum(diag) != number_occupied_bands:
#         print('Wrong filling factor for u_signal:', u_signal)
#         exit()
#
#     rho0constUp = (1 - alpha_rand) * np.diag(diag_caf) + alpha_rand * rhorand16
#     rho0constUm = (1 - alpha_rand) * np.diag(diag_caf) + alpha_rand * rhorand16
#
#     seed_dict = {'rho0constUp': rho0constUp, 'rho0constUm': rho0constUm}
#     return seed_dict
########################################################################################################################################################

if model_regime == 'full_range':
    rho0constUp = diag_full_regime(+1, number_occupied_bands)
    rho0constUm = diag_full_regime(-1, number_occupied_bands)
elif model_regime == 'near_zero_dielectric_field':
    rho0constUp = seed_asymmetric_calcs(nu)
    rho0constUm = seed_asymmetric_calcs(nu)

# print(seed_asymmetric_calcs(nu))
