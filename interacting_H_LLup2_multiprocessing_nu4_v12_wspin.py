import multiprocessing
import os
import time

import pandas as pd
from numpy import linalg as npla

from config import aux_dir_path, namecsv, itmax, nprocesses, t0, title
from input.parameters import *
from model.densities import rho0
from model.hamiltonians import full_hp, full_hm, full_hpm, full_hmp, h0p2, h0m2, h0p, h0m, idp, idps, idm, idms
from utils import eigen, frange, idxcalc, df_round, sort_dict, observable_to_csv


# from model.exchange_int import *

# now = datetime.now()
# current_time = now.strftime("%Y-%m-%d %H:%M:%S")
# current_time_file = now.strftime("%d%m%Y%H%M%S")

# dir_path = os.path.dirname(os.path.realpath(__file__))
# cwd = os.getcwd()
#
# aux_dir_path = cwd + '/aux2/'
# input_dir_path = cwd + '/input/'

# if os.path.isfile('screenlog.0'):
#     os.remove(aux_dir_path + 'screenlog.0')
#
# t0 = time.time()
#
# pd.set_option('display.max_columns', 300)
# pd.set_option('display.width', 1000)
# pd.set_option('display.max_rows', 500)
# # np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
# # np.set_printoptions( threshold=20, edgeitems=10, linewidth=140, formatter = dict( float = lambda x: "%.3g" % x ))  # float arrays %.3g
# np.set_printoptions( precision=5, suppress=True, threshold=20, edgeitems=10, linewidth=140, formatter={'float': '{: 0.3f}'.format})
#
# title = 'nu4_v12_wspin_random_rho0_hermitian_rho0phbroken_ep'+str(ep)+'_Zm' + str(Zm)  # copy v6
# print(title)
# # name=title+'.csv'
# namecsv = title + '.csv'
# # name = 'test.csv'
#
# machine=platform.node()
# infos='\n'+' Starting this script ('+title+'.py) at date/time: ' + current_time+'. \n'+ ' This script is running at: '+machine+', directory: ' + cwd +'\n'
#
# print(infos)
# with open(aux_dir_path + 'progress.txt', 'a') as f:
#     print(infos, file=f)

# idp_dic = {0:0,1:1,-2:2,2:3}
# def idp(n):
#     return idp_dic[n]
#
# idm_dic = {0:4,1:5,-2:6,2:7}
# def idm(n):
#     return idm_dic[n]
#
# idps_dic = {0:8,1:9,-2:10,2:11}
# def idps(n):
#     return idps_dic[n]
#
# idms_dic = {0:12,1:13,-2:14,2:15}
# def idms(n):
#     return idms_dic[n]

# def frange(start, end, inc):
#     return np.arange(start, end, inc).tolist()
#
#
# def eigen(A):
#     'returns eigenvalues and respective eigenvectors ordered by np.argsort'
#     eigenValues, eigenVectors = npla.eig(A)
#     #  eigenValuestmp=eigenValues
#     #  eigenVectorstmp=eigenVectors
#     idxfunc = np.argsort(eigenValues)
#     eigenValues = eigenValues[idxfunc]
#     eigenVectors = eigenVectors[:, idxfunc]
#     # return eigenValues.real, np.transpose(eigenVectors)
#     return eigenValues, np.transpose(eigenVectors)
#
# def idxcalc(idx, vecs1, vecs2, u):
#     lenth = len(vecs1)
#     idxtmp = idx.copy()
#     # idx=[x for x in range(4)].copy()
#     for i in range(lenth):
#         overlaptmp = np.abs(np.dot(vecs1[:, i], vecs2[:, idxtmp[i]]))
#         for j in range(lenth):
#             overlap = np.abs(np.dot(vecs1[:, i], vecs2[:, idxtmp[j]]))
#             # print('here')
#             if overlap > overlaptmp:
#                 idx[i] = idxtmp[j]
#                 overlaptmp = overlap
#                 print('u=%.3f' % (1000 * u),i,j)
#     # if idx != idxtmp:
#     #     print('crossing here', u * 10 ** 3)
#     return idx
#
# def check_hermitian(a, tol):
#     return np.all(np.abs(a - np.conjugate(a).T) < tol)
#
# def check_real(a, tol):
#     return np.all(np.abs(a - np.conjugate(a)) < tol)

# rho0rand = pd.read_csv(  aux_dir_path + 'rhoconstphbroken.csv', header=None).values.tolist()
# # print(rho0rand)
# rho = rho0rand

# mat=[[0.1 for i in range(16)] for j in range(16)]
# 0p-, 1p-, -2p-, 2p-, 0m-, 1m-, -2m-, 2m-,   0p+, 1p+, -2p+, 2p+, 0m+, 1m+, -2m+, 2m+
# rho0phbroken = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
#                 [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
#
# # 0p-, 1p-, -2p-, 2p-, 0m-, 1m-, -2m-, 2m-,   0p+, 1p+, -2p+, 2p+, 0m+, 1m+, -2m+, 2m+
# rho0ph = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
#           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
#
# # rho0rand = (1 - ep) * np.array(rho0phbroken) + ep * np.array(rho0ph)
#
#
# # rho0rand = [[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
# #        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
#
# # randvec = [np.random.uniform(0, 1) for i in range(16)]
# randvec = pd.read_csv(  input_dir_path + 'randvec.csv', header=None).values.tolist()[0]
# # print(randvec)
# randmatrix = np.outer(randvec, randvec)
# randeigenvalue, randeigenvector = eigen(randmatrix)
# rho0rand = sum(np.outer(randeigenvector[i, :], randeigenvector[i, :]) for i in range(nu))
#
# rho0 = (1 - ep) * np.array(rho0rand) + ep * np.array(rho0phbroken)
#
# res=check_hermitian(rho0, tol)
# resreal=check_real(rho0, tol)
#
# if res:
#     print('rho0 is hermitian')
#     if resreal:
#         print('rho0 is real')
#         rho0=rho0.real
#     else:
#         print('rho0 is not real')
#         exit()
#
# else:
#     print('rho0rand is not hermitian')
#     pd.DataFrame(rho0).to_csv( aux_dir_path + 'rho0_nsymmetric_' + namecsv, mode='w', index=False, header=False)
#     exit()
#
# pd.DataFrame(rho0).to_csv( aux_dir_path +  'rho0' + namecsv , mode='w', index=False, header=False)
# # # rho = sum(np.outer(np.absolute(np.array(randeigenvector[i])), np.absolute(np.array(randeigenvector[i]))) for i in
# # #            range(len(randeigenvector)))
# # # rho0rand = sum(np.outer(randeigenvector[i, :], randeigenvector[i, :]) for i in range(len(randeigenvector)))
# # rho0rand = sum( np.outer(randeigenvector[i, :], randeigenvector[i, :]) for i in range( nu ) )
# # # rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(nu));
# # rho = rho0rand
# # # pd.DataFrame(rho0rand).to_csv( aux_dir_path + 'rho0rand.csv', mode='w', index=False, header=False)
#
# eigenU = []
# # rhoU=[]
# rhoUdf = pd.DataFrame([])
# rhoUdf.to_csv( aux_dir_path + 'rhoU_tmp_' + namecsv, mode='w', index=False, header=False)
#
#
# def hp(n, nprime, s1, s2):
#     if (n == nprime) and (s1 == s2):
#         tmp = Eh * deltatb
#     else:
#         tmp = 0
#
#     def id1(n):
#         idxfsu = round((1 + s1) / 2)
#         idxfsd = round((1 - s1) / 2)
#         return idxfsu * idp(n) + idxfsd * idps(n)
#
#     def id2(n):
#         idxfsu = round((1 + s2) / 2)
#         idxfsd = round((1 - s2) / 2)
#         return idxfsu * idp(n) + idxfsd * idps(n)
#
#     res = -sum(
#         [Xskp(n2, nprime, n, n1, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
#             Xskp(n2, nprime, n, n1, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in
#             [-2, 2] for n2 in [-2, 2]]) - tmp
#     return res
#
#
# def hm(n, nprime, s1, s2):
#     if (n == nprime) and (s1 == s2):
#         tmp = Eh * deltatb
#     else:
#         tmp = 0
#
#     def id1(n):
#         idxfsu = round((1 + s1) / 2)
#         idxfsd = round((1 - s1) / 2)
#         return idxfsu * idm(n) + idxfsd * idms(n)
#
#     def id2(n):
#         idxfsu = round((1 + s2) / 2)
#         idxfsd = round((1 - s2) / 2)
#         return idxfsu * idm(n) + idxfsd * idms(n)
#
#     res = -sum(
#         [Xskm(n2, nprime, n, n1, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
#             Xskm(n2, nprime, n, n1, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in
#             [-2, 2] for n2 in [-2, 2]]) + tmp
#     return res
#
#
# def hpm(n, nprime, s1, s2):
#     def id1(n):
#         idxfsu = round((1 + s1) / 2)
#         idxfsd = round((1 - s1) / 2)
#         return idxfsu * idp(n) + idxfsd * idps(n)
#
#     def id2(n):
#         idxfsu = round((1 + s2) / 2)
#         idxfsd = round((1 - s2) / 2)
#         return idxfsu * idm(n) + idxfsd * idms(n)
#
#     res = -sum(
#         [Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in
#          [0, 1]] + [
#             Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[id1(n2)][id2(n1)] for n1 in [-2, 2] for n2 in
#             [-2, 2]])
#     return res
#
#
# def hmp(n, nprime, s1, s2):
#     def id1(n):
#         idxfsu = round((1 + s1) / 2)
#         idxfsd = round((1 - s1) / 2)
#         return idxfsu * idm(n) + idxfsd * idms(n)
#
#     def id2(n):
#         idxfsu = round((1 + s2) / 2)
#         idxfsd = round((1 - s2) / 2)
#         return idxfsu * idp(n) + idxfsd * idps(n)
#
#     res = -sum(
#         [Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [0, 1] for n2 in
#          [0, 1]] + [
#             Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[id1(n2)][id2(n1)] for n1 in [-2, 2] for n2 in
#             [-2, 2]])
#     return res
def hp(n, nprime, s1, s2):
    return full_hp(n, nprime, s1, s2, Eh, deltatb, eigenvectorp, rho)


def hm(n, nprime, s1, s2):
    return full_hm(n, nprime, s1, s2, Eh, deltatb, eigenvectorm, rho)


def hpm(n, nprime, s1, s2):
    return full_hpm(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


def hmp(n, nprime, s1, s2):
    return full_hmp(n, nprime, s1, s2, eigenvectorp, eigenvectorm, rho)


# tmp_obeservables = {}


def loopU(u):
    global eigenvectorp, eigenvectorm, deltatb, rho
    print('u=%.2f' % (u * 1000))
    rho = rho0

    eigenvaluep2, eigenvectorp2 = eigen(h0p2(u))
    eigenvaluem2, eigenvectorm2 = eigen(h0m2(u))
    eigenvaluep = [h0p(u)[0][0]] + [h0p(u)[1][1]] + eigenvaluep2.tolist()
    eigenvaluem = [h0m(u)[0][0]] + [h0m(u)[1][1]] + eigenvaluem2.tolist()

    h0 = np.diag(eigenvaluep + eigenvaluem + eigenvaluep + eigenvaluem)
    # print('here')
    mZm = k * np.diag([-Zm for i in range(8)] + [Zm for i in range(8)])
    # print(u * 1000, (1000 * np.diag(h0)).tolist())
    # pd.DataFrame([[u * 1000, pd.DataFrame((1000 * h0).round(decimals=3))]]).to_csv(aux_dir_path + 'h0_tmp_' + namecsv, mode='a', index=False, header=False)

    # eigenvectorp = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()]))
    # eigenvectorm = np.transpose(np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()]))

    eigenvectorp = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorp2.tolist()])
    eigenvectorm = np.array([[1, 0, 0, 0], [0, 1, 0, 0]] + [[0, 0] + x for x in eigenvectorm2.tolist()])

    # diff = 1;

    # def hps(n, nprime):
    #     if n == np:
    #         tmp = Eh * deltatb
    #     else:
    #         tmp = 0
    #     res = -sum(
    #         [Xskp(n2, nprime, n, n1, eigenvectorp) * rho[idps(n2)][idps(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [Xskp(n2, nprime, n, n1, eigenvectorp) * rho[idps(n2)][idps(n1)] for n1 in
    #                                                                                                            [-2, 2] for n2 in [-2, 2]]) - tmp
    #     return res
    #
    # def hms(n, nprime):
    #     if n == np:
    #         tmp = Eh * deltatb
    #     else:
    #         tmp = 0
    #     res = -sum(
    #         [Xskm(n2, nprime, n, n1, eigenvectorm) * rho[idms(n2)][idms(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [Xskm(n2, nprime, n, n1, eigenvectorm) * rho[idms(n2)][idms(n1)] for n1 in
    #                                                                                                            [-2, 2] for n2 in [-2, 2]]) + tmp
    #     return res
    #
    # def hpms(n, nprime):
    #     res = -sum(
    #         [Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[idps(n2)][idms(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
    #             Xdkpm(n2, nprime, n, n1, eigenvectorp, eigenvectorm) * rho[idps(n2)][idms(n1)] for n1 in [-2, 2] for n2 in [-2, 2]])
    #     return res
    #
    # def hmps(n, nprime):
    #     res = -sum(
    #         [Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[idms(n2)][idps(n1)] for n1 in [0, 1] for n2 in [0, 1]] + [
    #             Xdkmp(n2, nprime, n, n1, eigenvectorm, eigenvectorp) * rho[idms(n2)][idps(n1)] for n1 in [-2, 2] for n2 in [-2, 2]])
    #     return res
    it = 1;
    # while it < itmax:
    while it < itmax:
        deltatb = sum([rho[idp(n)][idp(n)] + rho[idps(n)][idps(n)] - rho[idm(n)][idm(n)] - rho[idms(n)][idms(n)] for n in setH])

        # hphpm = [[hp(n, nprime) for n in setH] + [hpm(n, nprime) for n in setH] for nprime in setH]
        # hmphm = [[hmp(n, nprime) for n in setH] + [hm(n, nprime) for n in setH] for nprime in setH]
        hphpm = [[hp(n, nprime, 1, 1) for nprime in setH] + [hpm(n, nprime, 1, 1) for nprime in setH] for n in setH]
        hmphm = [[hmp(n, nprime, 1, 1) for nprime in setH] + [hm(n, nprime, 1, 1) for nprime in setH] for n in setH]

        # for el in hmphm:
        #     hphpm.append(el)
        Hintup = np.vstack((hphpm, hmphm))

        hphpms = [[hp(n, nprime, -1, -1) for nprime in setH] + [hpm(n, nprime, -1, -1) for nprime in setH] for n in setH]
        hmphms = [[hmp(n, nprime, -1, -1) for nprime in setH] + [hm(n, nprime, -1, -1) for nprime in setH] for n in setH]

        Hintdown = np.vstack((hphpms, hmphms))

        hphpmsud = [[hp(n, nprime, 1, -1) for nprime in setH] + [hpm(n, nprime, 1, -1) for nprime in setH] for n in setH]
        hmphmsud = [[hmp(n, nprime, 1, -1) for nprime in setH] + [hm(n, nprime, 1, -1) for nprime in setH] for n in setH]

        Hintupdown = np.vstack((hphpmsud, hmphmsud))

        hphpmsdu = [[hp(n, nprime, -1, 1) for nprime in setH] + [hpm(n, nprime, -1, 1) for nprime in setH] for n in setH]
        hmphmsdu = [[hmp(n, nprime, -1, 1) for nprime in setH] + [hm(n, nprime, -1, 1) for nprime in setH] for n in setH]

        Hintdownup = np.vstack((hphpmsdu, hmphmsdu))

        Hint = k * np.vstack((np.hstack((Hintup, Hintupdown)), np.hstack((Hintdownup, Hintdown))))
        # Hint = k * np.vstack((np.hstack((Hintup, Hintupdown)), np.hstack((Hintupdown, Hintdown))))
        H = Hint + h0 + mZm  # np.add(Hint, h0)
        eigenvalue_loop, eigenvector_loop = eigen(H);
        rho = sum(np.outer(eigenvector_loop[i, :], eigenvector_loop[i, :]) for i in range(nu));

        it += 1

    # df = pd.DataFrame(rho.round(decimals=3))

    # .set_axis(['0', '1', '-2','2','0', '1', '-2','2'], axis=0).set_axis(['0', '1', '-2','2','0', '1', '-2','2'], axis=1)
    # pd.DataFrame([[u * 1000, df]]).to_csv(aux_dir_path + 'rhoU_tmp_' + namecsv, mode='a', index=False, header=False)
    # pd.DataFrame([[u, pd.DataFrame(rho )]]).to_csv( aux_dir_path + 'rhoU_tmp_'+namecsv, mode='a', index=False, header=False)

    # pd.DataFrame([[u * 1000, 1000 * k * Eh * deltatb]]).to_csv(aux_dir_path + 'Eh-delta_tmp_' + namecsv, mode='a',
    #                                                            index=False, header=False)

    # pd.DataFrame([[u, pd.DataFrame(1000 * Hint ) ]]).to_csv( aux_dir_path + 'Hint_tmp_'+namecsv, mode='a', index=False, header=False)
    # df = pd.DataFrame((1000 * Hint).round(decimals=3))  # .set_axis(['0', '1', '-2','2','0', '1', '-2','2'], axis=0).set_axis(['0', '1', '-2','2','0', '1', '-2','2'], axis=1)
    # pd.DataFrame([[u * 1000, df]]).to_csv(aux_dir_path + 'Hint_tmp_' + namecsv, mode='a', index=False, header=False)

    eigenvalue, eigenvector = npla.eig(H)
    #     #  eigenValuestmp=eigenValues
    #     #  eigenVectorstmp=eigenVectors
    idxfunc = np.argsort(eigenvalue)
    #
    eigenvalue = eigenvalue[idxfunc]
    eigenvector = eigenvector[:, idxfunc]

    ehf = - sum(
        [Hint[idp(n)][idp(nprime)] * rho[idp(nprime)][idp(n)] for n in setH for nprime in setH] + [Hint[idm(n)][idm(nprime)] * rho[idm(nprime)][idm(n)] for n in setH for nprime in
                                                                                                   setH] + [
            Hint[idps(n)][idps(nprime)] * rho[idps(nprime)][idps(n)] for n in setH for nprime in setH] + [Hint[idms(n)][idms(nprime)] * rho[idms(nprime)][idms(n)] for n in setH for
                                                                                                          nprime in setH] + [
            2 * Hint[idp(n)][idps(nprime)] * rho[idp(nprime)][idps(n)] for n in setH for nprime in setH] + [2 * Hint[idm(n)][idms(nprime)] * rho[idm(nprime)][idms(n)] for n in setH
                                                                                                            for nprime in
                                                                                                            setH])
    Et = sum([eigenvalue[i] for i in range(nu)]) + ehf

    dictionary = {'u': u * 1e3,
                  'eigenvalue': 1e3 * eigenvalue,
                  'eigenvector': eigenvector,
                  'Et': 1e3 * Et,
                  'h0': df_round(1e3 * h0),
                  'rhoU': df_round(rho),
                  'Eh_deltaU': 1e3 * k * Eh * deltatb,
                  'Hint': df_round(1e3 * Hint),
                  }
    # eigensystemU.append([u, (1000 * eigenvalue), eigenvector, 1000 * Et])
    # return u, (1000 * eigenvalue.real), eigenvector, 1000 * Et
    return u, (1000 * eigenvalue), eigenvector, 1000 * Et, dictionary


a_pool = multiprocessing.Pool(processes=nprocesses)

# calculates the function loopU over the range (U0minD, U0maxD, dU0D) using multiprocessing

eigensystemU = a_pool.map(loopU, frange(U0minD, U0maxD, dU0D))
# eigensystemU = a_pool.map(loopU, [-0.001,0])

# print(eigensystemU[0][1])
print(time.time() - t0)

# pd.read_csv(aux_dir_path + 'h0_tmp_' + namecsv, header=None).sort_values([0]).to_csv(aux_dir_path + 'h0_' + namecsv, index=False, header=False)
# os.remove(aux_dir_path + 'h0_tmp_' + namecsv)

# pd.read_csv(aux_dir_path + 'rhoU_tmp_' + namecsv, header=None).sort_values([0]).to_csv(aux_dir_path + 'rhoU_' + namecsv, index=False, header=False)
# os.remove(aux_dir_path + 'rhoU_tmp_' + namecsv)

# pd.read_csv(aux_dir_path + 'Eh-delta_tmp_' + namecsv, header=None).sort_values([0]).to_csv(aux_dir_path + 'Eh-deltaU_' + namecsv, index=False, header=False)
# os.remove(aux_dir_path + 'Eh-delta_tmp_' + namecsv)

# pd.read_csv(aux_dir_path + 'Hint_tmp_' + namecsv, header=None).sort_values([0]).to_csv(aux_dir_path + 'HintU_' + namecsv, index=False, header=False)
# os.remove(aux_dir_path + 'Hint_tmp_' + namecsv)

# print(tmp_obeservables)
# tmp_obeservables = sort_dict(tmp_obeservables)
# h0_df = pd.DataFrame([])
# rhoU_df = pd.DataFrame([])
# Eh_deltaU_df = pd.DataFrame([])
# HintU_df = pd.DataFrame([])
#
#
# h0_tmp = []
# rhoU_tmp = []
# Eh_deltaU_tmp = []
# HintU_tmp = []
# print(tmp_obeservables)


# for k,v in tmp_obeservables.items():
#     h0_tmp.append([k, v['h0_tmp']])
#     rhoU_tmp.append([k, v['rhoU_tmp']])
#     Eh_deltaU_tmp.append([k, v['Eh_deltaU_tmp']])
#     HintU_tmp.append([k, v['HintU_tmp']])

full_dict = {}

for el in eigensystemU:
    u,_,_,_,dict = el
    full_dict[dict['u']] = dict

tmp_obeservables = sort_dict(full_dict)

print(len(tmp_obeservables))

# def observable_to_csv(obeservables_dict, obeservable):
#     obeservable_list = []
#     for k, v in obeservables_dict.items():
#         obeservable_list.append([k, v[obeservable]])
#     obeservable_df = pd.DataFrame(obeservable_list)
#     obeservable_df.to_csv(aux_dir_path + obeservable + '_' + namecsv, index=False, header=False)


observable_to_csv(tmp_obeservables, 'h0')
observable_to_csv(tmp_obeservables, 'rhoU')
observable_to_csv(tmp_obeservables, 'Eh_deltaU')
observable_to_csv(tmp_obeservables, 'Hint')
observable_to_csv(tmp_obeservables, 'Et')

print(time.time() - t0)
idx = np.argsort(eigensystemU[0][1]).tolist()
eigenU = [[eigensystemU[0][0]] + eigensystemU[0][1][idx].tolist()]
eigenvector = eigensystemU[0][2][:, idx]

for i in range(len(eigensystemU) - 1):
    idx = idxcalc(idx, eigenvector, eigensystemU[i + 1][2], eigensystemU[i + 1][0])
    # idx = np.argsort(eigensystemU[i + 1][1]).tolist()
    # print(idx)
    eigenU.append([eigensystemU[i + 1][0]] + eigensystemU[i + 1][1][idx].tolist())
    eigenvector = eigensystemU[i + 1][2][:, idx]

df = pd.DataFrame(eigenU)
# print(df)
df.to_csv(aux_dir_path + namecsv, index=False, header=False)

eigenU = [[el.real for el in elU] for elU in eigenU]
# eigenU=np.array(eigenU).real.tolist

df = pd.DataFrame(eigenU)
# print(df)
df.to_csv(aux_dir_path + namecsv, index=False, header=False)

eigenUimag = [[elU[0]] + [el.imag for el in elU[1:len(elU)]] for elU in eigenU]

df = pd.DataFrame(eigenUimag)
# print(df)
df.to_csv(aux_dir_path + title + 'imag.csv', index=False, header=False)

EtU = [[eigensystemU[i][0], eigensystemU[i][3].real] for i in range(len(eigensystemU))]

pd.DataFrame(EtU).to_csv(aux_dir_path + 'EtU' + namecsv, index=False, header=False)
# print(df)

EtU = [[eigensystemU[i][0], eigensystemU[i][3].imag] for i in range(len(eigensystemU))]

pd.DataFrame(EtU).to_csv(aux_dir_path + 'EtUimag' + namecsv, index=False, header=False)
# print(df)


print('file ' + namecsv + ' saved')
# print('done')
print(" \n done in ", time.time() - t0)
