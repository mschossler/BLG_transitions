import numpy as np
from numpy import linalg as npla
import pandas as pd
from config import aux_dir_path, namecsv
# import collections

def frange(start, end, inc):
    return np.arange(start, end, inc).tolist()


def eigen(A):
    'returns eigenvalues and respective eigenvectors ordered by np.argsort'
    eigenValues, eigenVectors = npla.eig(A)
    #  eigenValuestmp=eigenValues
    #  eigenVectorstmp=eigenVectors
    idxfunc = np.argsort(eigenValues)
    eigenValues = eigenValues[idxfunc]
    eigenVectors = eigenVectors[:, idxfunc]
    # return eigenValues.real, np.transpose(eigenVectors)
    return eigenValues, np.transpose(eigenVectors)

def idxcalc(idx, vecs1, vecs2, u):
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
                print('u=%.3f' % (1000 * u),i,j)
    # if idx != idxtmp:
    #     print('crossing here', u * 10 ** 3)
    return idx

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