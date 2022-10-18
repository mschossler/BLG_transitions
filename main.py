import multiprocessing
import time
from datetime import datetime

import numpy as np
import pandas as pd

t0 = time.time()

from config import model_regime, results_dir_path, file_name_csv, nprocesses, bands, bands_LLm2_LL2, tests_mode
from input.parameters import U0minD, U0maxD, dU0D, nu, u_critical, parameters_to_save, replace_LLm2_LL2_low_u, save_folder_name
from utils import frange, sort_dict, observable_to_csv, idxcalc, transitions_energy_fermi_energy

a_pool = multiprocessing.Pool(processes=nprocesses)

if model_regime == 'full_range':
    from model.hartree_fock_and_regularization import loopU

    quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
    # quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
elif model_regime == 'no_LL2_mixing_and_asym':
    # print('here_condition_near_zero_calcs')
    # import time
    # time.sleep(.1)
    from model.hartree_fock_with_asymmetric_interactions import loopU0

    U0minD_tmp = max(U0minD, -u_critical)
    U0maxD_tmp = min(U0maxD, u_critical)
    if U0maxD_tmp > U0minD_tmp:
        quantities = a_pool.map(loopU0, frange(U0minD_tmp, U0maxD_tmp, dU0D))
    else:
        print('range of u not compatible with ' + model_regime + '. Change for a valid range (check u_critical)')
        exit()
    # quantities = a_pool.map(loopU0, [1e-3,2e-3])



def select_quantities_and_save_to_file(quantities_dict, model_regime_local):
    list_of_us = list(quantities_dict.keys())
    list_of_observables = list(quantities_dict[list_of_us[0]].keys())
    list_of_observables.remove('u')

    if bool(model_regime_local):
        # loop over u
        for dict_quantities in quantities_dict.values():
            # loop over keys of each dictionary of u
            for sub_key in list_of_observables:
                dict_quantities[sub_key + '_' + model_regime_local] = dict_quantities.pop(sub_key)
            # print(dict_quantities)

        list_of_observables = [quantity + '_' + model_regime_local for quantity in list_of_observables]

    # print(list_of_observables)
    for quantity in list_of_observables:
        observable_to_csv(quantities_dict, quantity)


def energies_and_observable_to_csv(quantities, model_regime_local=''):
    '''
    identify the highest weight base element in the eigenvector,
    and orders eigensystem for u propagation.
    calls function to export observables to csv
    '''

    quantities_dict = {}
    for dict_u in quantities:
        quantities_dict[dict_u['u']] = dict_u

    quantities_dict = sort_dict(quantities_dict)  # returns a sorted dictionary of {'u':dict_u}, dict_u are the dictionaries from loopU or loopU0 from the main calculations

    energies = []
    # allowed_transitions_nu = {}
    for k, v in quantities_dict.items():
        u_temp, eigenvalue_temp, eigenvector_temp = round(v['u'], 4), v['eigenvalue'], v['eigenvector']
        idx = idxcalc(eigenvector_temp)
        v['eigenvalue'] = eigenvalue_temp[idx]
        v['eigenvector'] = eigenvector_temp[:, idx]
        energies.append([u_temp] + v['eigenvalue'].tolist())

    select_quantities_and_save_to_file(quantities_dict, model_regime_local)

    energies_df = pd.DataFrame(energies, columns=['u'] + bands)
    return energies_df


energies_df = energies_and_observable_to_csv(quantities)

if model_regime == 'no_LL2_mixing_and_asym' and ((U0minD < -u_critical) or (U0maxD > u_critical) or replace_LLm2_LL2_low_u):

    def assign_HighFieldRange(energies, energies_high_u):
        # print(frange(U0minD, U0minD_tmp, dU0D) * 1e3)
        # energies.loc[energies['u'] > u_critical*1e3] = energies_high_u[np.isin(energies_high_u['u'],frange(U0maxD_tmp,U0maxD,dU0D)*1e3)].values
        # energies.loc[energies['u'] < -u_critical*1e3] = energies_high_u[np.isin(energies_high_u['u'],frange(U0minD,U0minD_tmp,dU0D)*1e3)].values

        energies_high_u_positive = energies_high_u[np.isin(energies_high_u['u'], frange(U0maxD_tmp, U0maxD, dU0D) * 1e3)]
        energies_high_u_negative = energies_high_u[np.isin(energies_high_u['u'], frange(U0minD, U0minD_tmp, dU0D) * 1e3)]
        # print(energies_high_u_negative)
        energies = pd.concat([energies_high_u_negative, energies, energies_high_u_positive]).reset_index(drop=True)
        # print(energies)
        # exit()
        return energies


    from model.hartree_fock_and_regularization import loopU

    a_pool = multiprocessing.Pool(processes=nprocesses)
    quantities_full_range = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
    model_regime_local = 'full_regime'
    energies_df_full_range = energies_and_observable_to_csv(quantities_full_range, model_regime_local)
    energies_df_full_range_small_u = energies_df_full_range[np.isin(energies_df_full_range['u'], frange(U0minD_tmp, U0maxD_tmp, dU0D) * 1e3)].reset_index(drop=True)

    if replace_LLm2_LL2_low_u:
        # energies_df_full_range = energies_df_full_range[np.isin(energies_df_full_range['u'], frange(U0minD_tmp, U0maxD_tmp, dU0D) * 1e3)]
        energies_df[bands_LLm2_LL2] = energies_df_full_range_small_u[bands_LLm2_LL2]
    energies_df = assign_HighFieldRange(energies_df, energies_df_full_range)


energies_df, transition_energy_df = transitions_energy_fermi_energy(energies_df, nu)  # add fermi_energy to energies_df
energies_df.round(8).to_csv(results_dir_path + 'energies_' + file_name_csv, index=False)
transition_energy_df.round(8).to_csv(results_dir_path + 'transitions_' + file_name_csv, index=False)

from visualization.plots import plot_energies, plot_transitions

# print and plot energies and transition energies
print(energies_df)
print(transition_energy_df)
plot_energies(energies_df, nu)
plot_transitions(transition_energy_df, nu)

# save parameters used in this simulation
parameters_to_save['script working duration'] = time.time() - t0
with open(results_dir_path + 'parameters.txt', 'w') as parameters_file:
    # with open('parameters.txt', 'w') as parameters_file:
    # for key in sorted(parameters_to_save.keys()):
    parameters_file.write('created on (Y-m-d): %s \n' % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for key in sorted(list(parameters_to_save.keys()), key=str.lower):
        parameters_file.write('%s = %s\n' % (key, parameters_to_save[key]))

print('files for nu=' + str(nu) + ' saved')

# save folder name for later ref in total energy comparations
if save_folder_name:
    with open(results_dir_path + '/../folder_list.txt', 'a') as f:
        print(tests_mode, file=f)

# #compare with previews results
# file1 = 'results/results_11092022/occupation_' + str(nu + 8) + '/energies_nu_' + str(nu) + '_to_compare.csv'
# file2 = 'results/results_' + current_date + '/occupation_' + str(nu + 8) + tests_mode + 'energies_' + file_name_csv
# print('same as previews results : %s' % filecmp.cmp(file1, file2))

print('running time for nu=%(nu)i: %(t).1fs' % {'t': time.time() - t0, 'nu': nu})
