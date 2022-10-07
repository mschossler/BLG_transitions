import numpy as np
import pandas as pd

# print('here_density')
# import time
# time.sleep(.1)

input_dir = 'input/'

if __name__ == "__main__":
    # setting path
    import sys

    sys.path.append('../')
    input_dir = '../input/'

from input.parameters import nu, alpha_rand_full_range, alpha_rand_asymmetric_calcs, file_seed  # , seed_large_u  # , model_regime
from config import bands, tol
from utils import eigen, remove_small_imag

class Density_Seed:
    def __init__(self, model_regime, nu):
        self.rho0constUp = None
        self.rho0constUm = None
        self.model_regime = model_regime
        self.nu = nu
        self.number_occupied_bands = nu + 8

    filling_order_Upositive = ['LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sdown',
                               'LL1_Kp_Sdown', 'LL1_Km_Sup', 'LL1_Kp_Sup', 'LL2_Kp_Sdown', 'LL2_Kp_Sup']

    filling_order_Unegative = ['LLm2_Km_Sdown', 'LLm2_Km_Sup', 'LLm2_Kp_Sdown', 'LLm2_Kp_Sup', 'LL0_Kp_Sdown', 'LL0_Km_Sdown', 'LL0_Kp_Sup', 'LL0_Km_Sup', 'LL1_Kp_Sdown',
                               'LL1_Km_Sdown', 'LL1_Kp_Sup', 'LL1_Km_Sup', 'LL2_Km_Sdown', 'LL2_Km_Sup']

    filling_order_Upositive_small_u = ['LLm2_Kp_Sdown',  # -7
                                       'LLm2_Kp_Sup',  # -6
                                       'LLm2_Km_Sdown',  # -5
                                       'LLm2_Km_Sup',  # -4
                                       'LL0_Km_Sdown',  # -3
                                       'LL0_Kp_Sdown',  # -2
                                       'LL1_Km_Sdown',  # -1
                                       'LL1_Kp_Sdown',  # 0
                                       'LL0_Kp_Sup',  # 1
                                       'LL0_Km_Sup',  # 2
                                       'LL1_Km_Sup',  # 3
                                       'LL1_Kp_Sup',  # 4
                                       'LL2_Kp_Sdown',  # 5
                                       'LL2_Kp_Sup']  # 6

    # filling_order_Upositive_small_u = ['LLm2_Kp_Sdown',  # -7
    #                                    'LLm2_Kp_Sup',  # -6
    #                                    'LLm2_Km_Sdown',  # -5
    #                                    'LLm2_Km_Sup',  # -4
    #                                    'LL0_Km_Sdown',  # -3
    #                                    'LL0_Kp_Sdown',  # -2
    #                                    'LL0_Kp_Sup',  # -1
    #                                    'LL0_Km_Sup',  # 0
    #                                    'LL1_Km_Sdown',  # 1
    #                                    'LL1_Kp_Sdown',  # 2
    #                                    'LL1_Km_Sup',  # 3
    #                                    'LL1_Kp_Sup',  # 4
    #                                    'LL2_Kp_Sdown',  # 5
    #                                    'LL2_Kp_Sup']  # 6

    # base_octet = ['LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL1_Km_Sdown', 'LL1_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sup', 'LL1_Kp_Sup']
    # seed_oct_dict = {-4: (0, 0, 0, 0, 0, 0, 0, 0),
    #                  -3: (1, 0, 0, 0, 0, 0, 0, 0),
    #                  -2: (1, 1, 0, 0, 0, 0, 0, 0),
    #                  -1: (1, 1, 1, 0, 0, 0, 0, 0),
    #                  0: (1, 1, 1, 1, 0, 0, 0, 0),
    #                  1: (1, 1, 1, 1, 1, 0, 0, 0),
    #                  2: (1, 1, 1, 1, 1, 1, 0, 0),
    #                  3: (1, 1, 1, 1, 1, 1, 1, 0),
    #                  4: (1, 1, 1, 1, 1, 1, 1, 1)
    #                  }

    def ramdom_16x16_density(self):
        if (self.model_regime == 'full_range') and (self.nu == 0) and file_seed:
            rhorand8 = pd.read_csv(input_dir + 'rho0phbroken.csv', header=None).values.tolist()

            rho0file8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
            zero_rho0file8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

            rhorand16 = rho0file8_zero + zero_rho0file8
        else:
            rhorand8 = np.random.rand(8, 8)
            eigenvalue, eigenvector = eigen(rhorand8)
            rhorand8 = np.real(sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(self.number_occupied_bands // 2)))

            rhorand8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
            zero_rhorand8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

            rhorand16 = rhorand8_zero + zero_rhorand8

        return rhorand16

    def diag_full_regime_high_u(self, u_signal):

        if u_signal >= 0:
            filling_order = self.filling_order_Upositive
        if u_signal < 0:
            filling_order = self.filling_order_Unegative
        diag = [(1 if band in filling_order[0:self.number_occupied_bands] else 0) for band in bands]
        if abs(sum(diag) - self.number_occupied_bands) > tol:
            print('Wrong filling factor for nu %i' % self.nu)
            exit()

        rhorand16 = self.ramdom_16x16_density()
        return (1 - alpha_rand_full_range) * np.diag(diag) + alpha_rand_full_range * rhorand16

    def diag_full_regime_small_u(self):

        diag = [(1 if band in self.filling_order_Upositive_small_u[0:self.number_occupied_bands] else 0) for band in bands]
        if abs(sum(diag) - self.number_occupied_bands) > tol:
            print('Wrong filling factor for nu %i' % self.nu)
            exit()

        rhorand16 = self.ramdom_16x16_density()
        if self.model_regime == 'full_range':
            alpha_rand = alpha_rand_full_range
        elif self.model_regime == 'no_LL2_mixing_and_asym':
            alpha_rand = alpha_rand_asymmetric_calcs
        return (1 - alpha_rand) * np.diag(diag) + alpha_rand * rhorand16

    # def seed_asymmetric_calcs(self):
    #
    #     diag = [(1 if band in self.filling_order_Upositive_small_u[0:self.number_occupied_bands] else 0) for band in bands]
    #     if abs(sum(diag) - number_occupied_bands) > tol:
    #         print('Wrong filling factor for nu %i' % self.nu)
    #         exit()
    #     rhorand16 = self.ramdom_16x16_density()
    #
    #     rho0const = (1 - alpha_rand_asymmetric_calcs) * np.diag(diag) + alpha_rand_asymmetric_calcs * rhorand16
    #
    #     return rho0const

    def assign_densities(self):
        self.rho0const_small_u = remove_small_imag(self.diag_full_regime_small_u())
        self.rhoRandom = self.ramdom_16x16_density()
        if self.model_regime == 'full_range':
            self.rho0constUp = remove_small_imag(self.diag_full_regime_high_u(+1))
            self.rho0constUm = remove_small_imag(self.diag_full_regime_high_u(-1))
        # elif self.model_regime == 'no_LL2_mixing_and_asym':
        #     self.rho0constUp = remove_small_imag(self.seed_asymmetric_calcs())
        #     self.rho0constUm = remove_small_imag(self.seed_asymmetric_calcs())


def density_by_model_regime(model_regime):
    densities = Density_Seed(model_regime, nu)
    densities.assign_densities()
    rho0constUp = densities.rho0constUp
    rho0constUm = densities.rho0constUm
    rhoRandom = densities.rhoRandom
    rho0const_small_u = densities.rho0const_small_u

    return {'rho0constUp': rho0constUp, 'rho0constUm': rho0constUm, 'rho0const_small_u': rho0const_small_u, 'rhoRandom':rhoRandom}
