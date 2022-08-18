import numpy as np
import pandas as pd
from numpy import linalg as npla

from config import aux_dir_path, file_name_csv, bands
from input.parameters import occupied_bands


def frange(start, end, inc):
    return np.arange(start, end, inc).tolist()


def eigen(A):
    'returns eigenvalues and respective eigenvectors ordered by np.argsort'
    eigenValues, eigenVectors = npla.eig(A)
    idxfunc = np.argsort(eigenValues)
    eigenValues = eigenValues[idxfunc]
    eigenVectors = eigenVectors[:, idxfunc]
    # return eigenValues.real, np.transpose(eigenVectors)
    return eigenValues, np.transpose(eigenVectors)


def check_hermitian(a, tol):
    return np.all(np.abs(a - np.conjugate(a).T) < tol)


def check_real(a, tol):
    return np.all(np.abs(a - np.conjugate(a)) < tol)


def df_round(observable):
    return pd.DataFrame(observable.round(decimals=3))


def sort_dict(dict):
    sorted_dict = {k: v for k, v in sorted(dict.items())}
    return sorted_dict


def observable_to_csv(obeservables_dict, obeservable):
    obeservable_list = []
    for k, v in obeservables_dict.items():
        obeservable_list.append([k, v[obeservable]])
    obeservable_df = pd.DataFrame(obeservable_list)
    obeservable_df.to_csv(aux_dir_path + obeservable + '_' + file_name_csv, index=False, header=False)


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


def tau_func(tau):
    tau = np.array(tau)
    temp = np.zeros((8, 8), dtype='float_')
    for i in range(2):
        temp[i, i] = tau[0, 0]
        temp[i, i + 4] = tau[0, 1]
        temp[i + 4, i] = tau[1, 0]
        temp[i + 4, i + 4] = tau[1, 1]
    return np.kron(np.eye(2, dtype=int), temp)


bands_by_sector = [list(t) for t in zip(*[iter(bands)] * 4)]

allowed_transitions = []
for sector in bands_by_sector:
    for band1 in sector:
        for band2 in sector:
            if (('LL1' in band1) and ('LL2' in band2)) or (('LLm2' in band1) and ('LL1' in band2)):
                allowed_transitions.append((band1, band2))


def transtions_and_fermi_energy_u(energy_u):
    u = energy_u['u']
    energy_u.drop('u', axis=0, inplace=True)
    occupied_states = energy_u.nsmallest(occupied_bands, keep='all')
    unoccupied_states = energy_u.nlargest(len(bands) - occupied_bands, keep='all')[::-1]  # must be reversed for right fermi_energy
    # print(occupied_states)
    # print(unoccupied_states)

    fermi_energy = (unoccupied_states.iloc[0] + occupied_states.iloc[-1]) / 2

    allowed_transtions_nu = []
    for two_allowed_bands in allowed_transitions:
        band1, band2 = two_allowed_bands
        if ('LL1' in band1) and (band1 in occupied_states.index):
            allowed_transtions_nu.append(two_allowed_bands)
        elif ('LL1' in band2) and (band2 in unoccupied_states.index):
            allowed_transtions_nu.append(two_allowed_bands)

    transition_energy_u_df = pd.DataFrame([])
    for allowed_transition_nu in allowed_transtions_nu:
        from_band, to_band = allowed_transition_nu
        transition_energy_u_df[from_band + '_to_' + to_band] = [energy_u[to_band] - energy_u[from_band]]

    transition_energy_u_df.index = [u]

    transtions_energy_and_fermi_energy_u_dict = {'u': u,
                                                 'fermi_energy': fermi_energy,
                                                 'allowed_transtions_u': transition_energy_u_df.columns,
                                                 'transtions_energy_u_df': transition_energy_u_df,
                                                 }

    return transtions_energy_and_fermi_energy_u_dict


def transition_fermi_energy(energies):
    transitions_energy = pd.DataFrame([])
    fermi_energy = []
    for ind in energies.index:
        energy_u = energies.loc[ind]
        #         print(energy_u,type(energy_u))
        #         print(transtions_and_fermi_energy_u(energy_u))
        transtions_and_fermi_energy_u_dict = transtions_and_fermi_energy_u(energy_u)

        transitions_energy_u = transtions_and_fermi_energy_u_dict['transtions_energy_u_df']
        fermi_energy.append(transtions_and_fermi_energy_u_dict['fermi_energy'])
        #         print(transtions_and_fermi_energy_u_dict)
        #         transitions_energy.join(transitions_energy_u)
        transitions_energy = pd.concat([transitions_energy, transitions_energy_u], axis=0)
    transitions_energy['u'] = transitions_energy.index.values
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
