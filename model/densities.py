import numpy as np
import pandas as pd

from config import aux_dir_path, namecsv, tol
from input.parameters import ep, nu
from utils import eigen, check_hermitian, check_real

# rho0rand = pd.read_csv(  aux_dir_path + 'rhoconstphbroken.csv', header=None).values.tolist()
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
rho0rand = sum(np.outer(randeigenvector[i, :], randeigenvector[i, :]) for i in range(nu))
rho0 = (1 - ep) * np.array(rho0rand) + ep * np.array(rho0phbroken)

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
    pd.DataFrame(rho0).to_csv('aux/' + 'rho0_nsymmetric_' + namecsv, mode='w', index=False, header=False)
    exit()

pd.DataFrame(rho0).to_csv(aux_dir_path + 'rho0' + namecsv, mode='w', index=False, header=False)
#
#
# rhoUdf = pd.DataFrame([])
# rhoUdf.to_csv( aux_dir_path + 'rhoU_tmp_' + namecsv, mode='w', index=False, header=False)