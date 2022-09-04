import numpy as np
import pandas as pd

from config import results_dir_path, tol, alpha_rand, file_name_csv
from input.parameters import number_occupied_bands, nu
from utils import eigen, check_hermitian, check_real

if nu == 4:
    spindownUp = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUp = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]

    spindownUm = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUm = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]


diag = [x[1] for x in spindownUp + spinupUp]
print('is the filling factor right (U>0)?: ', sum(diag) == number_occupied_bands)
rho0constUp = np.diag(diag)

diag = [x[1] for x in spindownUm + spinupUm]
print('is the filling factor right (U<0)?:', sum(diag) == number_occupied_bands)
rho0constUm = np.diag(diag)

# rho0rand = pd.read_csv(  results_dir_path + 'rhoconstphbroken.csv', header=None).values.tolist()
# # print(rho0rand)
# rho = rho0rand

# 0p-, 1p-, -2p-, 2p-, 0m-, 1m-, -2m-, 2m-,   0p+, 1p+, -2p+, 2p+, 0m+, 1m+, -2m+, 2m+
rho0phbroken = np.diag([1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]).tolist()

# 0p-, 1p-, -2p-, 2p-, 0m-, 1m-, -2m-, 2m-,   0p+, 1p+, -2p+, 2p+, 0m+, 1m+, -2m+, 2m+
rho0ph = np.diag([0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0]).tolist()

# randvec = [np.random.uniform(0, 1) for i in range(16)]
randvec = pd.read_csv('input/' + 'randvec.csv', header=None).values.tolist()[0]
randmatrix = np.outer(randvec, randvec)
randeigenvalue, randeigenvector = eigen(randmatrix)
rho0rand = sum(np.outer(randeigenvector[i, :], randeigenvector[i, :]) for i in range(number_occupied_bands))
rho0 = (1 - alpha_rand) * np.array(rho0rand) + alpha_rand * np.array(rho0phbroken)

res = check_hermitian(rho0, tol)
resreal = check_real(rho0, tol)

if res:
    print('rho0 is hermitian')
    if resreal:
        print('rho0 is real')
        rho0 = rho0.real
    else:
        print('rho0 is not real')
        exit()

else:
    print('rho0rand is not hermitian')
    pd.DataFrame(rho0).to_csv(results_dir_path + 'rho0_nsymmetric_' + file_name_csv, mode='w', index=False, header=False)
    exit()

pd.DataFrame(rho0).to_csv(results_dir_path + 'rho0' + file_name_csv, mode='w', index=False, header=False)
#
#
# rhoUdf = pd.DataFrame([])
# rhoUdf.to_csv( results_dir_path + 'rhoU_tmp_' + namecsv, mode='w', index=False, header=False)
