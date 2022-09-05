import multiprocessing
import time
from datetime import datetime

t0 = time.time()

import pandas as pd

from config import model_regime, results_dir_path, file_name_csv, nprocesses, bands, bands_LLm2_LL2, bands_LL2, bands_LLm2
from input.parameters import U0minD, U0maxD, dU0D, nu, u_critical, parameters_to_save, mode
from utils import frange, sort_dict, observable_to_csv, idxcalc, transitions_energy_fermi_energy

a_pool = multiprocessing.Pool(processes=nprocesses)
if model_regime == 'full_range':
    from model.hartree_fock_and_regularization import loopU

    quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
elif model_regime == 'near_zero_dielectric_field':
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


def energies_and_observable_to_csv(quantities, export_to_file=True):
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
    #
    # for quantity in ['h0', 'rhoU', 'Eh_deltaU', 'Hint', 'Et', 'eigenvector', 'exciton_energy', 'regmatrix']:
    #     observable_to_csv(quantities_dict, quantity)
    list_of_u = list(quantities_dict.keys())
    list_of_observables = list(quantities_dict[list_of_u[0]].keys())
    list_of_observables.remove('u')
    # print(list_of_observables)
    if export_to_file:
        for quantity in list_of_observables:
            observable_to_csv(quantities_dict, quantity)

    energies_df = pd.DataFrame(energies, columns=['u'] + bands)
    return energies_df


energies_df = energies_and_observable_to_csv(quantities)

# print('before \n', energies_df)
# energies_df_from_file = pd.read_csv('input/' + 'energies_nu_0_for_LLm2_LL2_HighFieldRange.csv')
# energies_df_from_file = energies_df_from_file[(energies_df_from_file['u'] >= U0minD * 1e3) & (energies_df_from_file['u'] < U0maxD * 1e3)].reset_index()
# print('file \n', energies_df_from_file[(energies_df_from_file['u'] >= U0minD * 1e3) & (energies_df_from_file['u'] < U0maxD * 1e3)])
# print('file_full \n', energies_df_from_file)
#
# print(type(energies_df), type(energies_df_from_file), bands_LLm2_LL2)
# energies_df[bands_LLm2_LL2] = energies_df_from_file[bands_LLm2_LL2]
# # test = energies_df[bands_LLm2_LL2]
# print('after \n', energies_df)
# exit()


if model_regime == 'near_zero_dielectric_field':

    print('approximation mode for LL2 and LLm2: %s' % mode)

    energies_df_from_file = pd.read_csv('input/' + 'energies_nu_0_for_LLm2_LL2_HighFieldRange.csv')
    energies_df_from_file = energies_df_from_file[(energies_df_from_file['u'] >= U0minD * 1e3) & (energies_df_from_file['u'] < U0maxD * 1e3)].reset_index(drop=True)


    # print('file \n', energies_df_from_file)

    def assign_HighFieldRange(energies, energies_high_u):
        energies.loc[energies['u'] > u_critical] = energies_high_u.loc[energies_high_u['u'] > u_critical].values
        energies.loc[energies['u'] < -u_critical] = energies_high_u.loc[energies_high_u['u'] < -u_critical].values
        return energies


    if mode == 'hartree_fock_and_regularization_calcs':
        from model.hartree_fock_and_regularization import loopU

        a_pool = multiprocessing.Pool(processes=nprocesses)
        quantities_full_range = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
        energies_df_full_range = energies_and_observable_to_csv(quantities_full_range, False)

        energies_df[bands_LLm2_LL2] = energies_df_full_range[bands_LLm2_LL2]
        energies_df = assign_HighFieldRange(energies_df, energies_df_full_range)

    elif mode == 'fast_none_interact':
        # energies_df_source_LLm2_LL2 = energies_df
        energies_df = assign_HighFieldRange(energies_df, energies_df_from_file)


    elif mode == 'fast_from_file':
        # energies_df[energies_df['u'] > u_critical] = energies_df_from_file[energies_df_from_file['u'] > u_critical]
        # energies_df[energies_df['u'] < -u_critical] = energies_df_from_file[energies_df_from_file['u'] < -u_critical]
        # print('before \n',energies_df)
        # print(type(energies_df),type(energies_df_from_file),bands_LLm2_LL2)
        energies_df[bands_LLm2_LL2] = energies_df_from_file[bands_LLm2_LL2]
        energies_df = assign_HighFieldRange(energies_df, energies_df_from_file)
        # test = energies_df[bands_LLm2_LL2]
        # print('after \n',test)
        # energies_df = assign_HighFieldRange(energies_df, energies_df_from_file)

    elif mode == 'fast_from_constant':
        constant_LLm2 = -14.14  # -(-56.701644+42.5611020) ##-56.701644 is average of LL-2 with full interactions
        constant_LL2 = -5.84  # 49.672084-55.512116000000006 ##49.672084 is average of LL2 with full interactions
        # energies_df_LLm2_LL2 = pd.DataFrame([])
        # energies_df_cosntant = energies_df
        energies_df[bands_LLm2] = energies_df[bands_LLm2] + constant_LLm2
        energies_df[bands_LL2] = energies_df[bands_LL2] + constant_LL2
        # print('before_assign\n', energies_df)
        energies_df = assign_HighFieldRange(energies_df, energies_df_from_file)
        # energies_df[energies_df['u'] > u_critical] = energies_df_from_file[energies_df_from_file['u'] > u_critical].values
        # energies_df[energies_df['u'] < -u_critical] = energies_df_from_file[energies_df_from_file['u'] < -u_critical].values

    # energies_df = assign_HighFieldRange(energies_df, energies_df_cosntant)
    # energies_df[bands_LLm2_LL2] = energies_df_cosntant[bands_LLm2_LL2]
# print('after')
energies_df, transition_energy_df = transitions_energy_fermi_energy(energies_df, nu)  # add fermi_energy to energies_df

energies_df.to_csv(results_dir_path + 'energies_' + file_name_csv, index=False)
transition_energy_df.to_csv(results_dir_path + 'transitions_' + file_name_csv, index=False)

from visualization.plots import plot_energies, plot_transitions

print(energies_df)
print(transition_energy_df)
plot_energies(energies_df, nu)
plot_transitions(transition_energy_df, nu)

with open(results_dir_path + 'parameters.txt', 'w') as parameters_file:
    # with open('parameters.txt', 'w') as parameters_file:
    # for key in sorted(parameters_to_save.keys()):
    parameters_file.write('created on (Y-m-d): %s \n' % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    for key in sorted(list(parameters_to_save.keys()), key=str.lower):
        parameters_file.write('%s = %s\n' % (key, parameters_to_save[key]))

print('files ' + file_name_csv + ' saved')
print('working duration for nu=%(nu)i: %(t).1fs' % {'t': time.time() - t0, 'nu': nu})
