import numpy as np
import pandas as pd

from input.parameters import nu, occupied_bands

filling_order_Upositive = ['-2p-', '-2p+', '-2m-', '-2m+', '0m-', '0p-', '0m+', '0p+', '1m-', '1p-', '1m+', '1p+', '2p-', '2p+']  # spin_polarz -> 0LL_polarz -> 1LL_spin_polarz
filling_order_Unegative = ['-2m-', '-2m+', '-2p-', '-2p+', '0p-', '0m-', '0p+', '0m+', '1p-', '1m-', '1p+', '1m+', '2m-', '2m-']  # spin_polarz -> 0LL_polarz -> 1LL_spin_polarz

if nu == 0:
    alpha = 0.6

    rhorand8 = pd.read_csv('input/' + 'rho0phbroken.csv', header=None).values.tolist()

    rho0file8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
    zero_rho0file8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

    rho0file16 = rho0file8_zero + zero_rho0file8
    rhorand = rho0file16

    spindownUp = [('0p-', 1), ('1p-', 0), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 0), ('-2m-', 1), ('2m-', 0)]
    spinupUp = [('0p+', 1), ('1p+', 0), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 0), ('-2m+', 1), ('2m+', 0)]

    spindownUm = [('0p-', 1), ('1p-', 0), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 0), ('-2m-', 1), ('2m-', 0)]
    spinupUm = [('0p+', 1), ('1p+', 0), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 0), ('-2m+', 1), ('2m+', 0)]

    diag = [x[1] for x in spindownUp + spinupUp]
    print('is the filling factor right (U>0)?: ', sum(diag) == occupied_bands)
    rho0constUp = (1 - alpha) * np.diag(diag) + alpha * rhorand

    diag = [x[1] for x in spindownUm + spinupUm]
    print('is the filling factor right?:  (U<0)', sum(diag) == occupied_bands)
    rho0constUm = (1 - alpha) * np.diag(diag) + alpha * rhorand

if nu == 4:
    spindownUp = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUp = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]

    spindownUm = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUm = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]

if nu != 0:
    diag = [occupation[1] for occupation in spindownUp + spinupUp]
    print('is the filling factor right (U>0)?: ', sum(diag) == occupied_bands)
    rho0constUp = np.diag(diag)

    diag = [occupation[1] for occupation in spindownUm + spinupUm]
    print('is the filling factor right (U<0)?:', sum(diag) == occupied_bands)
    rho0constUm = np.diag(diag)
