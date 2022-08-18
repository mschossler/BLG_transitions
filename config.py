import os
import platform
import time
from datetime import datetime

import numpy as np
import pandas as pd

from input.parameters import Zm, asym, x, alpha_H_oct_int, uz, uperp, alpha_state, alpha_rand, dens, occupied_bands

pd.set_option('display.max_columns', 300)
pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 500)
np.set_printoptions(precision=5, suppress=True, threshold=20, edgeitems=10, linewidth=140, formatter={'float': '{: 0.3f}'.format})

nprocesses = 12
itmax = 5
tol = 1e-8
setH = [0, 1, -2, 2]
now = datetime.now()
current_time_formated = now.strftime("%Y-%m-%d %H:%M:%S")
current_time = now.strftime("%d%m%Y%H%M%S")
current_date = now.strftime("%d%m%Y")
t0 = time.time()
cwd = os.getcwd()  # working directory
path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)
print(dir_path)

aux_dir_path = dir_path + '/results/' + current_date + '/occupation_' + str(occupied_bands) + '/'
input_dir_path = dir_path + '/input/'
# bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
bands = ['LL0_Kp_Sdown', 'LL1_Kp_Sdown', 'LLm2_Kp_Sdown', 'LL2_Kp_Sdown',
         'LL0_Km_Sdown', 'LL1_Km_Sdown', 'LLm2_Km_Sdown', 'LL2_Km_Sdown',
         'LL0_Kp_Sup', 'LL1_Kp_Sup', 'LLm2_Kp_Sup', 'LL2_Kp_Sup',
         'LL0_Km_Sup', 'LL1_Km_Sup', 'LLm2_Km_Sup', 'LL2_Km_Sup']

file_name = 'nu_' + str(occupied_bands)
file_name_csv = file_name + '.csv'
script_name = __file__
print(script_name)
# file_name = file_name + '.csv'
machine = platform.node()
infos = '\n' + ' Starting this script at date/time: ' + current_time + '. \n' + ' This script is running at: ' + machine + ', directory: ' + cwd + '\n'
folder_name = 'files_' + 'asym_' + str(round(asym, 2)) + '__itmax_' + str(round(itmax, 2)) + '__Zm_' + str(round(Zm * 1e3, 3)) + \
              '__alpha_H_oct_int_' + str(round(alpha_H_oct_int, 2)) + '__uz_' + str(round(uz * 1e3, 3)) + '__uperp_' + str(round(uperp * 1e3, 3)) + \
              '__x_' + str(round(x, 3)) + '__alpha_state_' + str(round(alpha_state, 3)) + '__alpha_rand_' + str(round(alpha_rand, 3)) + '__dens_' + str(round(dens, 1))

if not os.path.isdir(aux_dir_path):
    os.makedirs(aux_dir_path)

if os.path.isfile('screenlog.0'):
    os.remove(aux_dir_path + 'screenlog.0')

print(infos)
with open(aux_dir_path + 'progress.txt', 'a') as f:
    print(infos, file=f)
