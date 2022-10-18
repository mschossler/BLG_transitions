import multiprocessing
import subprocess
import time

t0_for_total = time.time()


def execute_nu(nu):
    subprocess.call(' python main.py ' + str(nu), shell=True)


pool = multiprocessing.Pool(processes=13)
pool.map(execute_nu, range(-6, 7))
# for nu in range(-6, 7):
# #     # for nu in (-6, -5, 5, 6):
# #     # for nu in range(-4, 4 + 1):
# #     # for nu in [0]:
# #     subprocess.call(' python main.py ' + str(nu), shell=True)
#     execute_nu(nu)
subprocess.call(' cd visualization \n python plots.py ', shell=True)
print('total working duration of execute_main: %.1fs' % (time.time() - t0_for_total))
subprocess.call('history -a ./results/history.txt', shell=True)
