import math

import numpy as np
import pandas as pd

model_dir = 'model/'
if __name__ == "__main__":
    # setting path
    import sys

    sys.path.append('../')
    model_dir = '../model/'

from config import tol
from input.parameters import x

decimals_round = int(-np.log10(tol))

xsFF = pd.read_csv(model_dir + 'xsFF.csv', names=['res', 'n2', 'np', 'n', 'n1'])
xsFF = xsFF.loc[(abs(xsFF['res']) > tol)].reset_index()

xsFF['res'] = xsFF['res'].round(decimals_round)
keys_xsFF = xsFF.drop(columns=['index', 'res']).to_records(index=False)
values_xsFF = xsFF['res'].round(decimals_round).to_list()
xsFF_dict = dict([(tuple(keys_xsFF[i]), values_xsFF[i]) for i in range(len(values_xsFF))])


xdFF = pd.read_csv(model_dir + 'xdFF.csv', names=['res', 'n2', 'np', 'n', 'n1'])
xdFF = xdFF.loc[(abs(xdFF['res']) > tol)].reset_index()

xdFF['res'] = xdFF['res'].round()
keys_xdFF = xdFF.drop(columns=['index', 'res']).to_records(index=False)
values_xdFF = xdFF['res'].round(decimals_round).to_list()
xdFF_dict = dict([(tuple(keys_xdFF[i]), values_xdFF[i]) for i in range(len(values_xdFF))])


def vplf(n, eigenvectorp):
    vplf_dict = {0: 0, 1: 0, -2: eigenvectorp[2][2], 2: eigenvectorp[3][2]}
    return vplf_dict[n]


# def vplf(n, eigenvectorp):
#     if n == 0:
#         return 0
#     elif n == 1:
#         return 0
#     elif n == -2:
#         return eigenvectorp[2][2]
#     elif n == 2:
#         return eigenvectorp[3][2]

def uplf(n, eigenvectorp):
    uplf_dict = {0: 1, 1: 1, -2: eigenvectorp[2][3], 2: eigenvectorp[3][3]}
    return uplf_dict[n]


# def uplf(n, eigenvectorp):
#     if n == 0:
#         return 1
#     elif n == 1:
#         return 1
#     elif n == -2:
#         return eigenvectorp[2][3]
#     elif n == 2:
#         return eigenvectorp[3][3]

def vmnf(n, eigenvectorm):
    vmnf_dict = {0: 0, 1: 0, -2: eigenvectorm[2][3], 2: eigenvectorm[3][3]}
    return vmnf_dict[n]


# def vmnf(n, eigenvectorm):
#     if n == 0:
#         return 0
#     elif n == 1:
#         return 0
#     elif n == -2:
#         return eigenvectorm[2][3]
#     elif n == 2:
#         return eigenvectorm[3][3]

def umnf(n, eigenvectorm):
    umnf_dict = {0: 1, 1: 1, -2: eigenvectorm[2][2], 2: eigenvectorm[3][2]}
    return umnf_dict[n]


# def umnf(n, eigenvectorm):
#     if n == 0:
#         return 1
#     elif n == 1:
#         return 1
#     elif n == -2:
#         return eigenvectorm[2][2]
#     elif n == 2:
#         return eigenvectorm[3][2]
def xsFFfunc(n2, nprime, n, n1):
    restmp = xsFF_dict.get((n2, nprime, n, n1), 0)
    return restmp


# def xsFFfunc(n2, nprime, n, n1):
#     restmp = xsFF.loc[
#         (xsFF['n2'] == n2) & (xsFF['np'] == nprime) & (xsFF['n'] == n) & (xsFF['n1'] == n1)].values.tolist()
#     if len(restmp) == 0:
#         restmp = 0
#     elif len(restmp) > 1:
#         print('error with', n2, nprime, n, n1)
#         exit()
#     else:
#         restmp = restmp[0][1]
#     return restmp


def Xskp(n2, nprime, n, n1, eigenvectorp):
    absn2 = abs(n2)
    absnp = abs(nprime)
    absn = abs(n)
    absn1 = abs(n1)

    if abs((eigenvectorp[2][2]) ** 2 + (eigenvectorp[2][3]) ** 2 - 1) > tol:
        print('error with eigenvectorp', n2, nprime, n, n1)
        exit()

    def vpl(n):
        return vplf(n, eigenvectorp)

    def upl(n):
        return uplf(n, eigenvectorp)

    res = upl(n2) * upl(nprime) * upl(n) * upl(n1) * xsFFfunc(absn2, absnp, absn, absn1) + \
          upl(n2) * upl(nprime) * vpl(n) * vpl(n1) * xsFFfunc(absn2, absnp, absn - 2, absn1 - 2) + \
          vpl(n2) * vpl(nprime) * upl(n) * upl(n1) * xsFFfunc(absn2 - 2, absnp - 2, absn, absn1) + \
          vpl(n2) * vpl(nprime) * vpl(n) * vpl(n1) * xsFFfunc(absn2 - 2, absnp - 2, absn - 2, absn1 - 2)

    res = round(res, decimals_round)
    return res


def Xskm(n2, nprime, n, n1, eigenvectorm):
    absn2 = abs(n2)
    absnp = abs(nprime)
    absn = abs(n)
    absn1 = abs(n1)

    if abs((eigenvectorm[2][2]) ** 2 + (eigenvectorm[2][3]) ** 2 - 1) > tol:
        print('error with eigenvectorm', n2, nprime, n, n1)
        exit()

    def vmn(n):
        return vmnf(n, eigenvectorm)

    def umn(n):
        return umnf(n, eigenvectorm)

    res = umn(n2) * umn(nprime) * umn(n) * umn(n1) * xsFFfunc(absn2, absnp, absn, absn1) + \
          umn(n2) * umn(nprime) * vmn(n) * vmn(n1) * xsFFfunc(absn2, absnp, absn - 2, absn1 - 2) + \
          vmn(n2) * vmn(nprime) * umn(n) * umn(n1) * xsFFfunc(absn2 - 2, absnp - 2, absn, absn1) + \
          vmn(n2) * vmn(nprime) * vmn(n) * vmn(n1) * xsFFfunc(absn2 - 2, absnp - 2, absn - 2, absn1 - 2)

    res = round(res, decimals_round)
    return res


def xdFFfunc(n2, nprime, n, n1):
    restmp = xdFF_dict.get((n2, nprime, n, n1), 0)
    return restmp


# def xdFFfunc(n2, nprime, n, n1):
#     restmp = xdFF.loc[
#         (xdFF['n2'] == n2) & (xdFF['np'] == nprime) & (xdFF['n'] == n) & (xdFF['n1'] == n1)].values.tolist()
#     if len(restmp) == 0:
#         restmp = 0
#     elif len(restmp) > 1:
#         print('error with', n2, nprime, n, n1)
#         exit()
#     else:
#         restmp = restmp[0][1]
#     return restmp


def Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm):
    absn2 = abs(n2)
    absnp = abs(nprime)
    absn = abs(n)
    absn1 = abs(n1)

    def vpl(n):
        return vplf(n, eigenvectorp)

    def upl(n):
        return uplf(n, eigenvectorp)

    def vmn(n):
        return vmnf(n, eigenvectorm)

    def umn(n):
        return umnf(n, eigenvectorm)

    res = upl(n2) * upl(nprime) * umn(n) * umn(n1) * xdFFfunc(absn2, absnp, absn, absn1) + \
          upl(n2) * upl(nprime) * vmn(n) * vmn(n1) * xdFFfunc(absn2, absnp, absn - 2, absn1 - 2) + \
          vpl(n2) * vpl(nprime) * umn(n) * umn(n1) * xdFFfunc(absn2 - 2, absnp - 2, absn, absn1) + \
          vpl(n2) * vpl(nprime) * vmn(n) * vmn(n1) * xdFFfunc(absn2 - 2, absnp - 2, absn - 2, absn1 - 2)

    res = round(res, decimals_round)
    return res


def Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp):
    absn2 = abs(n2)
    absnp = abs(nprime)
    absn = abs(n)
    absn1 = abs(n1)

    def vpl(n):
        return vplf(n, eigenvectorp)

    def upl(n):
        return uplf(n, eigenvectorp)

    def vmn(n):
        return vmnf(n, eigenvectorm)

    def umn(n):
        return umnf(n, eigenvectorm)

    res = umn(n2) * umn(nprime) * upl(n) * upl(n1) * xdFFfunc(absn2, absnp, absn, absn1) + \
          umn(n2) * umn(nprime) * vpl(n) * vpl(n1) * xdFFfunc(absn2, absnp, absn - 2, absn1 - 2) + \
          vmn(n2) * vmn(nprime) * upl(n) * upl(n1) * xdFFfunc(absn2 - 2, absnp - 2, absn, absn1) + \
          vmn(n2) * vmn(nprime) * vpl(n) * vpl(n1) * xdFFfunc(absn2 - 2, absnp - 2, absn - 2, absn1 - 2)

    res = round(res, decimals_round)
    return res


Xzs = 1;
Xzd = np.exp((x ** 2) / 2) * math.erfc(x / (np.sqrt(2)));
Xos = 3 / 4;
Xod = (1 / 8) * np.sqrt(2 / np.pi) * (-2 * x * (1 + x ** 2) + np.exp((x ** 2) / 2) * np.sqrt(2 * np.pi) * (3 + 2 * x ** 2 + x ** 4) * math.erfc(x / (np.sqrt(2))));
Xfs = 1 / 2;
Xfd = (1 / 4) * np.sqrt(2 / np.pi) * (-2 * x + np.exp((x ** 2) / 2) * np.sqrt(2 * np.pi) * (1 + x ** 2) * math.erfc(x / (np.sqrt(2))));
Xsts = 1 / 2;
Xstd = (1 / 4) * np.sqrt(2 / np.pi) * (2 * x - np.exp((x ** 2) / 2) * np.sqrt(2 * np.pi) * (-1 + x ** 2) * math.erfc(x / (np.sqrt(2))));
