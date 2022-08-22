import numpy as np
import pandas as pd

if __name__ == "__main__":
    # setting path
    import sys
    sys.path.append('../')

from input.parameters import nu, number_occupied_bands
from config import bands

# filling_order_Upositive = ['-2p-', '-2p+', '-2m-', '-2m+', '0m-', '0p-', '0m+', '0p+', '1m-', '1p-', '1m+', '1p+', '2p-', '2p+']  # spin_polarz -> 0LL_polarz -> 1LL_spin_polarz
# filling_order_Unegative = ['-2m-', '-2m+', '-2p-', '-2p+', '0p-', '0m-', '0p+', '0m+', '1p-', '1m-', '1p+', '1m+', '2m-', '2m-']  # spin_polarz -> 0LL_polarz -> 1LL_spin_polarz

filling_order_Upositive = ['LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sdown',
                           'LL1_Kp_Sdown', 'LL1_Km_Sup', 'LL1_Kp_Sup', 'LL2_Kp_Sdown', 'LL2_Kp_Sup']
filling_order_Unegative = ['LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LL0_Kp_Sdown', 'LL0_Km_Sdown', 'LL0_Kp_Sup', 'LL0_Km_Sup', 'LL1_Kp_Sdown',
                           'LL1_Km_Sdown', 'LL1_Kp_Sup', 'LL1_Km_Sup', 'LL2_Km_Sdown', 'LL2_Km_Sup']


def diag(u_signal, number_occupied_bands):
    if u_signal >= 0:
        filling_order = filling_order_Upositive
    if u_signal < 0:
        filling_order = filling_order_Unegative
    diag = [(1 if band in filling_order[0:number_occupied_bands] else 0) for band in bands]
    if sum(diag) != number_occupied_bands:
        print('is the filling factor right?:', sum(diag) == number_occupied_bands, 'u_signal:', u_signal)
        exit()
    return np.diag(diag)

rho0constUp = diag(+1, number_occupied_bands)
rho0constUm = diag(-1, number_occupied_bands)

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

if nu == 0:
    alpha = 0.6

    rhorand8 = pd.read_csv('input/' + 'rho0phbroken.csv', header=None).values.tolist()

    rho0file8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
    zero_rho0file8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

    rho0file16 = rho0file8_zero + zero_rho0file8
    rhorand = rho0file16

    # spindownUp = [('0p-', 1), ('1p-', 0), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 0), ('-2m-', 1), ('2m-', 0)]
    # spinupUp = [('0p+', 1), ('1p+', 0), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 0), ('-2m+', 1), ('2m+', 0)]
    #
    # spindownUm = [('0p-', 1), ('1p-', 0), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 0), ('-2m-', 1), ('2m-', 0)]
    # spinupUm = [('0p+', 1), ('1p+', 0), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 0), ('-2m+', 1), ('2m+', 0)]
    #
    # diag = [x[1] for x in spindownUp + spinupUp]
    # if sum(diag) != number_occupied_bands:
    #     print('is the filling factor right (U>0)?: ', sum(diag) == number_occupied_bands)
    #     exit()
    rho0constUp = (1 - alpha) * rho0constUp + alpha * rhorand
    #
    # diag = [x[1] for x in spindownUm + spinupUm]
    # if sum(diag) != number_occupied_bands:
    #     print('is the filling factor right (U<0)?:', sum(diag) == number_occupied_bands)
    #     exit()
    rho0constUm = (1 - alpha) * rho0constUm + alpha * rhorand

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
