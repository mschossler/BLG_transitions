import subprocess

for nu in range(-4, 5):
    # for nu in (-6, -5, 5, 6):
    subprocess.call(' python main.py ' + str(nu), shell=True)
