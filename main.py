import multiprocessing
import time

t0 = time.time()

import pandas as pd

from config import model_regime, aux_dir_path, file_name_csv, nprocesses, bands, bands_LL2
from input.parameters import U0minD, U0maxD, dU0D, nu
from utils import frange, sort_dict, observable_to_csv, idxcalc, transitions_energy_fermi_energy

# from model.hartree_fock_and_regularization import loopU
# from model.hartree_fock_with_asymmetric_interactions import loopU0

a_pool = multiprocessing.Pool(processes=nprocesses)

if model_regime == 'full_range':
    from model.hartree_fock_and_regularization import loopU

    quantities = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
elif model_regime == 'near_zero_dielectric_field':
    # print('here_condition_near_zero_calcs')
    # import time
    # time.sleep(.1)
    from model.hartree_fock_with_asymmetric_interactions import loopU0

    quantities = a_pool.map(loopU0, frange(U0minD, U0maxD, dU0D))
    # quantities = a_pool.map(loopU0, [1e-3,2e-3])


def energies_and_observable_to_csv(quantities):
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
    # print(list_of_observables)
    for quantity in list_of_observables:
        observable_to_csv(quantities_dict, quantity)

    energies_df = pd.DataFrame(energies, columns=['u'] + bands)
    return energies_df


energies_df = energies_and_observable_to_csv(quantities)
if model_regime == 'near_zero_dielectric_field':
    # mode = 'hartree_fock_and_regularization_calcs'
    # mode = 'fast_from_file'
    mode = 'fast_from_constant'

    if mode == 'hartree_fock_and_regularization_calcs':
        from model.hartree_fock_and_regularization import loopU

        a_pool = multiprocessing.Pool(processes=nprocesses)
        quantities_full_range = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
        energies_df_full_range = energies_and_observable_to_csv(quantities_full_range)

    elif mode == 'fast_from_file':
        energies_df_full_range = pd.read_csv('input/' + 'energies_nu_0_for_LLm2_and_LL2.csv')

    elif mode == 'fast_from_constant':
        constant_LLm2 = -14.14  # -(-56.701644+42.5611020) ##-56.701644 is average of LL-2 with full interactions
        constant_LL2 = -5.84  # 49.672084-55.512116000000006 ##49.672084 is average of LL2 with full interactions
        energies_df_full_range = pd.DataFrame([])
        energies_df_full_range[['LLm2_Kp_Sdown', 'LLm2_Km_Sdown', 'LLm2_Kp_Sup', 'LLm2_Km_Sup']] = energies_df[['LLm2_Kp_Sdown', 'LLm2_Km_Sdown', 'LLm2_Kp_Sup',
                                                                                                                'LLm2_Km_Sup']] + constant_LLm2
        energies_df_full_range[['LL2_Kp_Sdown', 'LL2_Km_Sdown', 'LL2_Kp_Sup', 'LL2_Km_Sup']] = energies_df[
                                                                                                   ['LL2_Kp_Sdown', 'LL2_Km_Sdown', 'LL2_Kp_Sup', 'LL2_Km_Sup']] + constant_LL2

    energies_df[bands_LL2] = energies_df_full_range[bands_LL2]

energies_df, transition_energy_df = transitions_energy_fermi_energy(energies_df, nu)  # add fermi_energy to energies_df

energies_df.to_csv(aux_dir_path + 'energies_' + file_name_csv, index=False)
transition_energy_df.to_csv(aux_dir_path + 'transitions_' + file_name_csv, index=False)

from visualization.plots import plot_energies, plot_transitions

print(energies_df)
print(transition_energy_df)
plot_energies(energies_df, nu)
plot_transitions(transition_energy_df, nu)

print('files ' + file_name_csv + ' saved')
print('working duration for nu=%(nu)i: %(t).1fs' % {'t': time.time() - t0, 'nu': nu})
