import os
import platform
from datetime import datetime

import numpy as np
import pandas as pd

from input.parameters import Zm, asym, x, alpha_H_oct_int, uz, uperp, nu, number_occupied_bands, model_regime, tests_mode  # , alpha_state, dens,

pd.set_option('display.max_columns', 300)
pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 500)
np.set_printoptions(precision=5, suppress=True, threshold=20, edgeitems=10, linewidth=140, formatter={'float': '{: 0.3f}'.format})

if abs(nu) > 4:
    model_regime = 'full_range'

if model_regime == 'full_range':
    print('model: %s' % model_regime)
elif model_regime == 'near_zero_dielectric_field':
    print('model: %s' % model_regime)

nprocesses = 10
# print(np.finfo(float).eps)
tol = 1e-12
setH = [0, 1, -2, 2]
now = datetime.now()
current_time_formated = now.strftime("%Y-%m-%d %H:%M:%S")
current_time = now.strftime("%d%m%Y%H%M%S")
current_date = now.strftime("%d%m%Y")

# current_date = 'server'

cwd = os.getcwd()  # working directory
path = os.path.abspath(__file__)
dir_path = os.path.dirname(path)
# print(dir_path)

if tests_mode == 'on':
    tests_mode = '/' + current_time + '__uz_' + str(round(uz * 1e3, 3)) + '__uperp_' + str(round(uperp * 1e3, 3)) + '/'
elif tests_mode == 'off':
    tests_mode = '/'

results_dir_path = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands) + tests_mode
results_dir_path_plot_vs_nu = dir_path + '/results/results_' + current_date + '/vs_nu/'
# results_dir_path = dir_path + '/results_old/occupation_' + str(number_occupied_bands) + '/'
# print(results_dir_path_plot_vs_nu)

input_dir_path = dir_path + '/input/'
# bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
bands = ['LL0_Kp_Sdown', 'LL1_Kp_Sdown', 'LLm2_Kp_Sdown', 'LL2_Kp_Sdown',
         'LL0_Km_Sdown', 'LL1_Km_Sdown', 'LLm2_Km_Sdown', 'LL2_Km_Sdown',
         'LL0_Kp_Sup', 'LL1_Kp_Sup', 'LLm2_Kp_Sup', 'LL2_Kp_Sup',
         'LL0_Km_Sup', 'LL1_Km_Sup', 'LLm2_Km_Sup', 'LL2_Km_Sup']
index_octet_on_bands_oct = [bands.index(band) for band in bands if '2' not in band]
# print(index_octet_on_bands_oct)
base_octet = ['LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL1_Km_Sdown', 'LL1_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sup', 'LL1_Kp_Sup']  #
indexes_octet_on_bands = [bands.index(band) for band in base_octet]
# print(indexes_octet_on_bands)
# print(base_octet)
bands_LLm2_LL2 = [band for band in bands if '2' in band]
bands_LL2 = [band for band in bands if 'LL2' in band]
bands_LLm2 = [band for band in bands if 'm2' in band]

file_name = 'nu_' + str(nu)
file_name_csv = file_name + '.csv'
script_name = __file__
# print(script_name)
# file_name = file_name + '.csv'
machine = platform.node()
infos = '\n' + ' Starting this script at date/time: ' + current_time_formated + '. \n' + ' This script is running at: ' + machine + ', directory: ' + cwd + '\n'
folder_name = 'files_' + 'asym_' + str(round(asym, 2)) + '__Zm_' + str(round(Zm * 1e3, 3)) + \
              '__alpha_H_oct_int_' + str(round(alpha_H_oct_int, 2)) + '__uz_' + str(round(uz * 1e3, 3)) + '__uperp_' + str(round(uperp * 1e3, 3)) + \
              '__x_' + str(round(x, 3)) + '__alpha_state_'  # + str(round(alpha_state, 3)) + '__dens_' + str(round(dens, 1))

if not os.path.isdir(results_dir_path):
    os.makedirs(results_dir_path)

if not os.path.isdir(results_dir_path_plot_vs_nu):
    os.makedirs(results_dir_path_plot_vs_nu)

if os.path.isfile('screenlog.0'):
    os.remove(results_dir_path + 'screenlog.0')

print(infos)
with open(results_dir_path + 'progress.txt', 'a') as f:
    print(infos, file=f)
