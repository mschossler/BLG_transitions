import multiprocessing
import time

t0 = time.time()

import pandas as pd

from config import aux_dir_path, file_name_csv, nprocesses, bands
from input.parameters import U0minD, U0maxD, dU0D, nu
from model.hartree_fock_and_regularization import loopU
from utils import frange, sort_dict, observable_to_csv, idxcalc, transitions_energy_fermi_energy

a_pool = multiprocessing.Pool(processes=nprocesses)
quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))

quantities_dict = {}
for dict_u in quantities:
    quantities_dict[dict_u['u']] = dict_u

quantities_dict = sort_dict(quantities_dict)

energies = []
# allowed_transitions_nu = {}
for k, v in quantities_dict.items():
    u_temp, eigenvalue_temp, eigenvector_temp = round(v['u'], 4), v['eigenvalue'], v['eigenvector']
    idx = idxcalc(eigenvector_temp)
    v['eigenvalue'] = eigenvalue_temp[idx]
    v['eigenvector'] = eigenvector_temp[:, idx]
    # allowed_transitions_nu[v['u']] = allowed_transitions(v['u'],v['eigenvalue'])
    energies.append([u_temp] + v['eigenvalue'].tolist())

energies_df = pd.DataFrame(energies, columns=['u'] + bands)

for quantity in ['h0', 'rhoU', 'Eh_deltaU', 'Hint', 'Et', 'eigenvector', 'exciton_energy', 'regmatrix']:
    observable_to_csv(quantities_dict, quantity)

energies_df, transition_energy_df = transitions_energy_fermi_energy(energies_df)  # add fermi_energy to energies_df

energies_df.to_csv(aux_dir_path + 'energies_' + file_name_csv, index=False)
transition_energy_df.to_csv(aux_dir_path + 'transitions_' + file_name_csv, index=False)

from visualization.plots import plot_energies, plot_transitions

print(energies_df)
print(transition_energy_df)
plot_energies(energies_df)
plot_transitions(transition_energy_df)

print('files ' + file_name_csv + ' saved')
print('working duration for nu=%(nu)i: %(t).1fs' % {'t': time.time() - t0, 'nu': nu})
