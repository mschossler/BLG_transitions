import multiprocessing
import time

from config import aux_dir_path, namecsv, nprocesses, t0, bands
import pandas as pd
from input.parameters import U0minD, U0maxD, dU0D
from utils import frange, sort_dict, observable_to_csv, idxcalc, idxcalc_v0

from model.hartree_fock_regularization import loopU

a_pool = multiprocessing.Pool(processes=nprocesses)
quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))

quantities_dict = {}
for dict in quantities:
    quantities_dict[dict['u']] = dict

quantities_dict = sort_dict(quantities_dict)

energies = []
for k, v in quantities_dict.items():
    u_temp, eigenvalue_temp, eigenvector_temp = v['u'], v['eigenvalue'], v['eigenvector']
    idx = idxcalc(eigenvector_temp)
    v['eigenvector'] = eigenvector_temp[:, idx]
    energies.append([u_temp] + eigenvalue_temp[idx].tolist())

energies_df = pd.DataFrame(energies, columns=['u'] + bands)
energies_df.to_csv(aux_dir_path + namecsv, index=False)

for quantity in ['h0', 'rhoU', 'Eh_deltaU', 'Hint', 'Et','eigenvector']:
    observable_to_csv(quantities_dict, quantity)

print(time.time() - t0)
print('file ' + namecsv + ' saved')
print(" \n done in ", time.time() - t0)

from visualization.plots import plot_energies

print(energies_df)
plot_energies(energies_df)