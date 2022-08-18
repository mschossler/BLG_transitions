import numpy as np
import pandas as pd

from input.parameters import nu, occupied_bands

if nu == 4:
    spindownUp = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUp = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]

    spindownUm = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUm = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]

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

diag = [occupation[1] for occupation in spindownUp + spinupUp]
print('is the filling factor right (U>0)?: ', sum(diag) == occupied_bands)
rho0constUp = np.diag(diag)

diag = [occupation[1] for occupation in spindownUm + spinupUm]
print('is the filling factor right (U<0)?:', sum(diag) == occupied_bands)
rho0constUm = np.diag(diag)
