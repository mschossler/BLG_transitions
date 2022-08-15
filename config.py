import os
from datetime import datetime
import time
import platform
import pandas as pd
import numpy as np
from input.parameters import ep, Zm, nu

now = datetime.now()
current_time = now.strftime("%Y-%m-%d %H:%M:%S")
current_time_file = now.strftime("%d%m%Y%H%M%S")

cwd = os.getcwd()

aux_dir_path = cwd + '/aux2/'
input_dir_path = cwd + '/input/'

nprocesses=40
itmax = 3
tol = 1e-8


if os.path.isfile('screenlog.0'):
    os.remove(aux_dir_path + 'screenlog.0')

t0 = time.time()

pd.set_option('display.max_columns', 300)
pd.set_option('display.width', 1000)
pd.set_option('display.max_rows', 500)
# np.set_printoptions(formatter={'float': '{: 0.3f}'.format})
# np.set_printoptions( threshold=20, edgeitems=10, linewidth=140, formatter = dict( float = lambda x: "%.3g" % x ))  # float arrays %.3g
np.set_printoptions( precision=5, suppress=True, threshold=20, edgeitems=10, linewidth=140, formatter={'float': '{: 0.3f}'.format})

title = 'nu4_v12_wspin_random_rho0_hermitian_rho0phbroken_ep'+str(ep)+'_Zm' + str(Zm)  # copy v6
print(title)
# name=title+'.csv'
namecsv = title + '.csv'
# name = 'test.csv'

machine=platform.node()
infos='\n'+' Starting this script ('+title+'.py) at date/time: ' + current_time+'. \n'+ ' This script is running at: '+machine+', directory: ' + cwd +'\n'

print(infos)
with open(aux_dir_path + 'progress.txt', 'a') as f:
    print(infos, file=f)