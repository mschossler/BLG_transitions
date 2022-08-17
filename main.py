import multiprocessing
import time

import pandas as pd

from config import aux_dir_path, namecsv, nprocesses, t0, bands
from input.parameters import U0minD, U0maxD, dU0D
from model.hartree_fock_and_regularization import loopU
from utils import frange, sort_dict, observable_to_csv, idxcalc, transition_fermi_energy

a_pool = multiprocessing.Pool(processes=nprocesses)
quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))

quantities_dict = {}
for dict_u in quantities:
    quantities_dict[dict_u['u']] = dict_u

quantities_dict = sort_dict(quantities_dict)


energies = []
# allowed_transitions_nu = {}
for k, v in quantities_dict.items():
    u_temp, eigenvalue_temp, eigenvector_temp = v['u'], v['eigenvalue'], v['eigenvector']
    idx = idxcalc(eigenvector_temp)
    v['eigenvalue'] = eigenvalue_temp[idx]
    v['eigenvector'] = eigenvector_temp[:, idx]
    # allowed_transitions_nu[v['u']] = allowed_transitions(v['u'],v['eigenvalue'])
    energies.append([u_temp] + v['eigenvalue'].tolist())

energies_df = pd.DataFrame(energies, columns=['u'] + bands)
energies_df.to_csv(aux_dir_path + namecsv, index=False)

for quantity in ['h0', 'rhoU', 'Eh_deltaU', 'Hint', 'Et', 'eigenvector', 'exciton_energy', 'regmatrix']:
    observable_to_csv(quantities_dict, quantity)

energies_df, transition_energy_df = transition_fermi_energy(energies_df)
transition_energy_df.to_csv(aux_dir_path + 'trans_' + namecsv, index=False)

from visualization.plots import plot_energies, plot_transitions

# transition_energy_all = pd.DataFrame([])
# for t in allowed_transitions:
#     transition_energy_df = transition_energy(energies_df, t)
#     transition_energy_all = pd.concat([transition_energy_all, transition_energy_df], axis=1)
#     # print(transition_energy_df)
#
# transition_energy_all['u'] = transition_energy_all.index


print(energies_df)
print(transition_energy_df)
plot_energies(energies_df)
plot_transitions(transition_energy_df)

print(time.time() - t0)
print('file ' + namecsv + ' saved')
print(" \n done in ", time.time() - t0)
