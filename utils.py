import numpy as np
import pandas as pd
from numpy import linalg as npla

from config import aux_dir_path, namecsv, bands


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
    obeservable_df.to_csv(aux_dir_path + obeservable + '_' + namecsv, index=False, header=False)


# ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
base = np.identity(16)


def idxcalc(vecs1):
    lenth = len(base)
    # idxtmp = idx.copy()
    # idx=[x for x in range(4)].copy()
    idx = []
    for i in range(lenth):
        # overlaptmp = npy.abs(npy.dot(vecs1[:, i], base[:, i]))
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


def idxcalc_v0(idx, vecs1, vecs2, u):
    lenth = len(vecs1)
    idxtmp = idx.copy()
    # idx=[x for x in range(4)].copy()
    for i in range(lenth):
        overlaptmp = np.abs(np.dot(vecs1[:, i], vecs2[:, idxtmp[i]]))
        for j in range(lenth):
            overlap = np.abs(np.dot(vecs1[:, i], vecs2[:, idxtmp[j]]))
            # print('here')
            if overlap > overlaptmp:
                idx[i] = idxtmp[j]
                overlaptmp = overlap
                print('u=%.3f' % (1000 * u), i, j)
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


def transition_energy(energies, two_allowed_bands: (str, str), extra_terms=0):
    from_band, to_band = two_allowed_bands
    if isinstance(energies, pd.DataFrame):
        pass
    else:
        print('transition_energy: invalid format for energies, pandas.DataFrame is expected')
        exit()
    transtion_label = from_band + '_to_' + to_band
    transition_energy_df = pd.DataFrame([])
    transition_energy_df[transtion_label] = energies[to_band] - energies[from_band]
    transition_energy_df.index = energies['u'].values
    # transition_tuple = transition_energy_df  # + extra_terms)
    return transition_energy_df


bands_by_sector = [list(t) for t in zip(*[iter(bands)] * 4)]

allowed_transitions = []
for sector in bands_by_sector:
    for band1 in sector:
        for band2 in sector:
            if (('LL1' in band1) and ('LL2' in band2)) or (('LLm2' in band1) and ('LL1' in band2)):
                allowed_transitions.append((band1, band2))

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
