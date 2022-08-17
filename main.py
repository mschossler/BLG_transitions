import multiprocessing
import time

import pandas as pd

from config import aux_dir_path, namecsv, nprocesses, t0, bands
from input.parameters import U0minD, U0maxD, dU0D
from model.hartree_fock_and_regularization import loopU
from utils import frange, sort_dict, observable_to_csv, idxcalc, transition_energy, allowed_transitions

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
    v['eigenvalue'] = eigenvalue_temp[idx]
    v['eigenvector'] = eigenvector_temp[:, idx]
    energies.append([u_temp] + v['eigenvalue'].tolist())

energies_df = pd.DataFrame(energies, columns=['u'] + bands)
energies_df.to_csv(aux_dir_path + namecsv, index=False)

for quantity in ['h0', 'rhoU', 'Eh_deltaU', 'Hint', 'Et', 'eigenvector', 'exciton_energy', 'regmatrix']:
    observable_to_csv(quantities_dict, quantity)

from visualization.plots import plot_energies

transition_energy_dic = {}
for t in allowed_transitions:
    transition_label, energy = transition_energy(energies_df, t)
    transition_energy_dic[transition_label] = energy

print(transition_energy_dic)

print(energies_df)
plot_energies(energies_df)

print(time.time() - t0)
print('file ' + namecsv + ' saved')
print(" \n done in ", time.time() - t0)
