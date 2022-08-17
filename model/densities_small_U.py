import numpy as np
from input.parameters import nu, occupied_bands


if nu == 4:
    spindownUp = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUp = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]

    spindownUm = [('0p-', 1), ('1p-', 1), ('-2p-', 1), ('2p-', 0), ('0m-', 1), ('1m-', 1), ('-2m-', 1), ('2m-', 0)]
    spinupUm = [('0p+', 1), ('1p+', 1), ('-2p+', 1), ('2p+', 0), ('0m+', 1), ('1m+', 1), ('-2m+', 1), ('2m+', 0)]


diag = [occupation[1] for occupation in spindownUp + spinupUp]
print('is the filling factor right (U>0)?: ', sum(diag) == occupied_bands)
rho0constUp = np.diag(diag)

diag = [occupation[1] for occupation in spindownUm + spinupUm]
print('is the filling factor right (U<0)?:', sum(diag) == occupied_bands)
rho0constUm = np.diag(diag)
