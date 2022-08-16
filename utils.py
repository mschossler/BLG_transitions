import numpy as np
import pandas as pd
from numpy import linalg as npla
from config import aux_dir_path, namecsv

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

base = np.identity(16)

def idxcalc(vecs1):
    lenth = len(base)
    idx = []
    for i in range(lenth):
        overlap = 0
        for j in range(lenth):
            overlaptmp = np.abs(np.dot(base[i, :], vecs1[j, :]))
            if overlaptmp > overlap:
                idxtemp = j
                overlap = overlaptmp
        idx.append(idxtemp)
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