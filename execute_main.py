import multiprocessing
import subprocess
import time

t0_for_total = time.time()
a_pool = multiprocessing.Pool(processes=4)


def execute_nu(nu):
    subprocess.call(' python main.py ' + str(nu), shell=True)


a_pool.map(execute_nu, range(-6, 7))
# for nu in range(-6, 7):
#     # for nu in (-6, -5, 5, 6):
#     # for nu in range(-4, 4 + 1):
#     # for nu in [0]:
#     subprocess.call(' python main.py ' + str(nu), shell=True)

print('total working duration of execute_main: %.1fs' % (time.time() - t0_for_total))
