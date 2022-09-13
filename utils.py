import numpy as np
import pandas as pd
from numpy import linalg as npla

from config import results_dir_path, file_name_csv, bands, tol


# import globals


def frange(start, end, inc):
    return np.round(np.arange(start, end, inc).tolist(), 6)


# print(tol_machine_epsilon)
def remove_small_imag(scalar):
    # print(scalar)
    return np.real_if_close(scalar, tol_machine_epsilon)


energy_imag_max = 0


def check_if_complex(energy_iterable, u, nu):
    global energy_imag_max
    complex_part = False
    # print(u)
    for energy in energy_iterable:
        if abs(np.imag(energy)) > 0:
            energy_imag_max = max(energy_imag_max, abs(energy.imag * 1e3))
            complex_part = True
            # print(u)
            warning_dict = {'u': u * 1e3, 'nu': nu, 'imag': energy.imag * 1e3, 'max_imag': energy_imag_max}
            warning = '-> warning: none-zero imaginary part energy %(imag).4fmeV for nu=%(nu)i, u=%(u).2fmeV. max abs imaginary part: %(max_imag).4fmeV' % warning_dict
            print(warning)
            with open(results_dir_path + 'warnings.txt', 'a') as f:
                print(warning, file=f)
            # return complex_part
    return complex_part


def eigen(A):
    'returns eigenvalues and respective eigenvectors ordered by np.argsort'
    eigenValues, eigenVectors = npla.eig(A)
    idxfunc = np.argsort(eigenValues)
    eigenValues = eigenValues[idxfunc]
    eigenVectors = eigenVectors[:, idxfunc]
    # return eigenValues.real, np.transpose(eigenVectors)
    return remove_small_imag(eigenValues), remove_small_imag(np.transpose(eigenVectors))


def check_hermitian(a, tol):
    return np.all(np.abs(a - np.conjugate(a).T) < tol)


def check_real(a, tol):
    return np.all(np.abs(a - np.conjugate(a)) < tol)


def df_round(observable):
    return pd.DataFrame(observable.round(decimals=3))


def sort_dict(dict):
    sorted_dict = {k: v for k, v in sorted(dict.items())}
    return sorted_dict


# returns real part if imag<tol
tol_machine_epsilon = tol / np.finfo(float).eps


def observable_to_csv(obeservables_dict, obeservable):
    "takes a dictionary of shape {'u',dict_u} for various u's and exports file of dict_u[observable] to csv"
    obeservable_list = []
    for k, v in obeservables_dict.items():
        obeservable_list.append([round(k, 2), v[obeservable]])
    obeservable_df = pd.DataFrame(obeservable_list)
    obeservable_df = obeservable_df.applymap(remove_small_imag)
    obeservable_df.round(8).to_csv(results_dir_path + obeservable + '_' + file_name_csv, index=False, header=False)


# ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
base = np.identity(16)


def idxcalc(vecs1):
    lenth = len(base)
    idx = []
    for i in range(lenth):
        overlap = 0
        for j in range(lenth):
            overlaptmp = np.abs(np.dot(base[i, :], vecs1[j, :]))
            # print('here')
            if overlaptmp > overlap:
                idxtemp = j
                overlap = overlaptmp
        idx.append(idxtemp)
        # print('u=%.3f' % (u), i, j)
    # if idx != idxtmp:
    #     print('crossing here', u * 10 ** 3)
    return idx


def nonedimmerp(vectors):
    newvectors = []
    for x in vectors:
        tmp = [x[0], x[3]]
        newvectors.append((np.array(tmp) / npla.norm(tmp)).tolist())
    return np.array(newvectors)


def nonedimmerm(vectors):
    newvectors = []
    for x in vectors:
        tmp = [-x[3], x[0]]
        newvectors.append((np.array(tmp) / npla.norm(tmp)).tolist())
    return np.array(newvectors)


bands_by_sector = [list(t) for t in zip(*[iter(bands)] * 4)]

allowed_transitions = []
for sector in bands_by_sector:
    for band1 in sector:
        for band2 in sector:
            if (('LL1' in band1) and ('LL2' in band2)) or (('LLm2' in band1) and ('LL1' in band2)):  # only states with LL1 and LL2, or LL1 and LLm2
                allowed_transitions.append((band1, band2))


# def remove_imaginary_energy(energy_u_series,u,nu):
#     def check_if_real(scalar):
#         if np.imag(scalar)>tol:
#             # print(energy_u_series, u)
#             # exit()
#             print('none-zero imaginary part for nu=%(nu)i, u=%(u).2fmeV '% {'u': (u * 1e3), 'nu': nu} )
#             # scalar = np.real(scalar)
#             # return True
#         # else:
#             # return False
#         # return scalar.real
#     energy_u_series.map(check_if_real, na_action='ignore')


# scalar = np.real(scalar)
# return True
# else:
# return False
# return scalar.real
# energy_u_series.map(check_if_real, na_action='ignore')
# def check_if_real(scalar,u,nu):
#     if np.imag(scalar)>tol:
#         print('none-zero imaginary part for nu=%(nu)i, u=%(u).2fmeV '% {'u': (u * 1e3), 'nu': nu} )
#         # scalar = np.real(scalar)
#     return scalar.real

# .apply(check_if_real, result_type='ignore', args=(u,nu))


def transitions_energy_and_fermi_energy_u(energy_u, nu):
    u = energy_u['u']
    # print(u)
    energy_u.drop('u', axis=0, inplace=True)
    # if check_if_complex(energy_u, u, nu):
    #     # print('here')
    #     # print('here0', check_if_complex(energy_u, ind, nu))
    #     energy_u = energy_u.map(np.real)
    #     # print('here0', check_if_complex(energy_u, ind, nu))
    #     # print('here2')
    # energy_u = energy_u.map(remove_small_imag)
    # energy_u
    # if check_if_complex(energy_u,u,nu):
    #     print('here')
    #     print('here0',check_if_complex(energy_u,u,nu))
    #     energy_u = energy_u.map(np.real)
    #     print('here0',check_if_complex(energy_u,u,nu))
    #     print('here2')
    # print('here3', check_if_complex(energy_u, u, nu))
    # print(u,energy_u)
    # energy_u = energy_u.apply(check_if_real, args=(u, nu)).copy()
    number_occupied_bands = nu + 8
    occupied_states = energy_u.nsmallest(number_occupied_bands, keep='all')
    unoccupied_states = energy_u.nlargest(len(bands) - number_occupied_bands, keep='all')[::-1]  # must be reversed for right fermi_energy
    # print(occupied_states)
    # print(unoccupied_states)

    fermi_energy = (unoccupied_states.iloc[0] + occupied_states.iloc[-1]) / 2

    allowed_transitions_nu = []
    for two_allowed_bands in allowed_transitions:
        band1, band2 = two_allowed_bands
        if ('LL1' in band1) and (band1 in occupied_states.index) and (band2 in unoccupied_states.index):
            allowed_transitions_nu.append(two_allowed_bands)
        elif ('LL1' in band2) and (band2 in unoccupied_states.index) and (band1 in occupied_states.index):
            allowed_transitions_nu.append(two_allowed_bands)

    transition_energy_u_df = pd.DataFrame([])
    for allowed_transition_nu in allowed_transitions_nu:
        from_band, to_band = allowed_transition_nu
        # if (u >= 0) and (from_band in filling_order_Upositive[0:number_occupied_bands]): #doesn't work for nu = 5 and 6
        transition_energy_u_df[from_band + '_to_' + to_band] = [energy_u[to_band] - energy_u[from_band]]
        # if (u < 0) and (from_band in filling_order_Unegative[0:number_occupied_bands]): #doesn't work for nu = 5 and 6
        transition_energy_u_df[from_band + '_to_' + to_band] = [energy_u[to_band] - energy_u[from_band]]

    transition_energy_u_df.index = [u]

    transitions_energy_and_fermi_energy_u_dict = {'u': u,
                                                  'fermi_energy': fermi_energy,
                                                  'allowed_transitions_u': transition_energy_u_df.columns,
                                                  'transitions_energy_u_df': transition_energy_u_df,
                                                  }

    return transitions_energy_and_fermi_energy_u_dict


# pd.options.mode.chained_assignment = None  # default='warn'
# def check_if_complex(energy_u_series, u, nu):
#     complex_part = False
#     print(u)
#     for energy in energy_u_series.values:
#         if abs(energy.imag) > 0:
#             complex_part = True
#             # print(u)
#             print('-> warning: none-zero imaginary part energy %(imag).4fmeV for nu=%(nu)i, u=%(u).2fmeV ' % {'u': u, 'nu': nu, 'imag': energy.imag})
#             # return complex_part
#     return complex_part


def transitions_energy_fermi_energy(energies, nu):
    transitions_energy = pd.DataFrame([])
    fermi_energy = []
    # for ind in energies.index:
    for _, energy_u in energies.iterrows():
        # print(index, energy_u['u'])
        # energy_u = energies.loc[ind]
        # u = np.real(energy_u['u'])
        # if check_if_complex(energy_u, u, nu):
        #     # print('here')
        #     # print('here0', check_if_complex(energy_u, ind, nu))
        #     energy_u = energy_u.map(np.real)
        #     # print('here0', check_if_complex(energy_u, ind, nu))
        #     # print('here2')
        # energy_u = energy_u.map(remove_small_imag)
        # print('here3', check_if_complex(energy_u, ind, nu))
        # print(ind, energy_u)
        #         print(energy_u,type(energy_u))
        #         print(transitions_energy_and_fermi_energy_u(energy_u))
        transitions_energy_and_fermi_energy_u_dict = transitions_energy_and_fermi_energy_u(energy_u, nu)
        # print(transitions_energy_and_fermi_energy_u_dict)
        transitions_energy_u = transitions_energy_and_fermi_energy_u_dict['transitions_energy_u_df']
        fermi_energy.append(transitions_energy_and_fermi_energy_u_dict['fermi_energy'])
        #         print(transitions_energy_and_fermi_energy_u_dict)
        #         transitions_energy.join(transitions_energy_u)
        transitions_energy = pd.concat([transitions_energy, transitions_energy_u], axis=0)
    # transitions_energy['u'] = transitions_energy.index.values
    transitions_energy.insert(0, 'u', transitions_energy.index.values)
    transitions_energy = transitions_energy.reset_index(drop=True)
    energies = energies.applymap(np.real)
    energies['fermi_energy'] = fermi_energy
    return energies, transitions_energy


###########

markers_list = ["2", "1", "1", "2", "1", "1", "2", "2"]
color_list = ['tab:blue', 'tab:orange', 'tab:red', 'tab:purple', 'tab:brown', 'tab:gray', 'tab:olive', 'tab:cyan']

transitions_plot_dic = {}
for two_allowed_bands in allowed_transitions:
    band1, band2 = two_allowed_bands
    if 'Sdown' in band1:
        spin_sign = 'downarrow'
    elif 'Sup' in band1:
        spin_sign = 'uparrow'

    if 'Kp' in band1:
        valley_sign = '+'
        valley_shape = '-'
    elif 'Km' in band1:
        valley_sign = '-'
        valley_shape = '--'

    if ('LL1' in band1):
        LL_band1 = '1'
    elif ('LLm2' in band1):
        LL_band1 = '-2'
    elif ('LL1' in band1):
        LL_band1 = '2'

    if ('LL1' in band2):
        LL_band2 = '1'
    elif ('LLm2' in band2):
        LL_band2 = '-2'
    elif ('LL2' in band2):
        LL_band2 = '2'

    if (('LL1' in band1) and ('LL2' in band2)) or (('LLm2' in band1) and ('LL1' in band2)):
        transitions_plot_dic[band1 + '_to_' + band2] = {'color': color_list.pop(0),
                                                        'marker_shape': markers_list.pop(0),
                                                        'line_shape': valley_shape,
                                                        'label': r'$' + LL_band1 + '\mathrm{K}^{' + valley_sign + '}\\' + spin_sign + '\ \longrightarrow\ ' +
                                                                 LL_band2 + '\mathrm{K}^{' + valley_sign + '}\\' + spin_sign + '$'
                                                        }
