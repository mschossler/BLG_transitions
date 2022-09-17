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

from input.parameters import nu, number_occupied_bands, alpha_rand_full_range, alpha_rand_asymmetric_calcs, file_seed  # , seed_large_u  # , model_regime
from config import bands, base_octet, tol
from utils import eigen, remove_small_imag


# remove_small_imag = np.vectorize(remove_small_imag)
# def array_map(x):
#     return np.vectorize(remove_small_imag)

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
    # base_octet = ['LL0_Km_Sdown', 'LL0_Kp_Sdown', 'LL1_Km_Sdown', 'LL1_Kp_Sdown', 'LL0_Km_Sup', 'LL0_Kp_Sup', 'LL1_Km_Sup', 'LL1_Kp_Sup']
    seed_oct_dict = {-4: (0, 0, 0, 0, 0, 0, 0, 0),
                     -3: (1, 0, 0, 0, 0, 0, 0, 0),
                     -2: (1, 1, 0, 0, 0, 0, 0, 0),
                     -1: (1, 1, 1, 0, 0, 0, 0, 0),
                     0: (1, 1, 1, 1, 0, 0, 0, 0),
                     1: (1, 1, 1, 1, 1, 0, 0, 0),
                     2: (1, 1, 1, 1, 1, 1, 0, 0),
                     3: (1, 1, 1, 1, 1, 1, 1, 0),
                     4: (1, 1, 1, 1, 1, 1, 1, 1)
                     }

    def ramdom_16x16_density(self):
        # if (self.model_regime == 'full_range_error') and (self.nu == 0):
        if (self.model_regime == 'full_range') and (self.nu == 0) and file_seed:
            # print('inside condition (model_regime == full_range) and (nu == 0)')
            rhorand8 = pd.read_csv(input_dir + 'rho0phbroken.csv', header=None).values.tolist()

            rho0file8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
            zero_rho0file8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

            rhorand16 = rho0file8_zero + zero_rho0file8
        else:
            # print('here_random_seed')
            # import time
            # time.sleep(.1)
            rhorand8 = np.random.rand(8, 8)
            eigenvalue, eigenvector = eigen(rhorand8)
            # rhorand8 = sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(int(number_occupied_bands / 2)))
            rhorand8 = np.real(sum(np.outer(eigenvector[i, :], eigenvector[i, :]) for i in range(self.number_occupied_bands // 2)))

            rhorand8_zero = np.pad(rhorand8, ((0, 8), (0, 8)), mode='constant')
            zero_rhorand8 = np.pad(rhorand8, ((8, 0), (8, 0)), mode='constant')

            rhorand16 = rhorand8_zero + zero_rhorand8
        # print('______________________________________________________________________________________________________')

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

        # if u_signal >= 0:
        #     rho0const = rho0constUp
        # if u_signal < 0:
        #     filling_order = filling_order_Unegative
        # seed_dict = {'rho0constUp': rho0constUp, 'rho0constUm': rho0constUm}
        rhorand16 = self.ramdom_16x16_density()
        # if nu==0:
        #     return np.diag([1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0])
        # return np.diag([0.618, 0.618, 1, 0, 0.618, 0.618, 1, 0, 0.382, 0.382, 1, 0, 0.382, 0.382, 1, 0])
        return (1 - alpha_rand_full_range) * np.diag(diag) + alpha_rand_full_range * rhorand16

    def diag_full_regime_small_u(self):
        if self.nu <= -4:
            seed_oct = self.seed_oct_dict[-4]
            number_occupied_bands_local = self.nu + 8
        elif abs(self.nu) < 4:
            seed_oct = self.seed_oct_dict[nu]
            number_occupied_bands_local = sum(seed_oct) + 4
        elif self.nu >= 4:
            seed_oct = self.seed_oct_dict[4]
            number_occupied_bands_local = self.nu + 8
        # number_occupied_bands_local = sum(seed_oct) + 4
        # number_occupied_bands = nu + 8
        occupied_octet_states = [base_octet[i] for i in range(8) if seed_oct[i]]

        filling_order = self.filling_order_Upositive[0:4] + occupied_octet_states + self.filling_order_Upositive[-2::1]

        diag = [(1 if band in filling_order[0:number_occupied_bands_local] else 0) for band in bands]

        diag = [(1 if band in self.filling_order_Upositive_small_u[0:self.number_occupied_bands] else 0) for band in bands]
        if abs(sum(diag) - self.number_occupied_bands) > tol:
            print('Wrong filling factor for nu %i' % self.nu)
            exit()

        # if u_signal >= 0:
        #     rho0const = rho0constUp
        # if u_signal < 0:
        #     filling_order = filling_order_Unegative
        # seed_dict = {'rho0constUp': rho0constUp, 'rho0constUm': rho0constUm}
        rhorand16 = self.ramdom_16x16_density()
        # if nu==0:
        #     return np.diag([1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0])
        # return np.diag([0.618, 0.618, 1, 0, 0.618, 0.618, 1, 0, 0.382, 0.382, 1, 0, 0.382, 0.382, 1, 0])
        return (1 - alpha_rand_full_range) * np.diag(diag) + alpha_rand_full_range * rhorand16

    #     # return np.diag(diag)

    # def diag_full_regime_v2(self):
    #     seed_oct = self.seed_oct_dict[self.nu]
    #     number_occupied_bands_local = sum(seed_oct) + 4
    #     # number_occupied_bands = nu + 8
    #     occupied_octet_states = [base_octet[i] for i in range(8) if seed_oct[i]]
    #
    #     filling_order = self.filling_order_Upositive[0:4] + occupied_octet_states + self.filling_order_Upositive[-2::1]
    #
    #     diag = [(1 if band in filling_order[0:number_occupied_bands_local] else 0) for band in bands]
    #     # diag= [0.618, 0.618, 1, 0, 0.618, 0.618, 1, 0, 0.382, 0.382, 1, 0, 0.382, 0.382, 1, 0]
    #     if abs(sum(diag) - number_occupied_bands) > tol:
    #         print('Wrong filling factor for nu %i' % self.nu)
    #         exit()
    #     # print('here_seed_asymmetric')
    #     # import time
    #     # time.sleep(.1)
    #     rhorand16 = self.ramdom_16x16_density()
    #
    #     rho0const = (1 - alpha_rand_full_range) * np.diag(diag) + alpha_rand_full_range * rhorand16
    #
    #     return rho0const

    def seed_asymmetric_calcs(self):
        seed_oct = self.seed_oct_dict[self.nu]
        number_occupied_bands_local = sum(seed_oct) + 4
        # number_occupied_bands = nu + 8
        occupied_octet_states = [base_octet[i] for i in range(8) if seed_oct[i]]

        filling_order = self.filling_order_Upositive[0:4] + occupied_octet_states + self.filling_order_Upositive[-2::1]

        diag = [(1 if band in filling_order[0:number_occupied_bands_local] else 0) for band in bands]
        # diag= [0.618, 0.618, 1, 0, 0.618, 0.618, 1, 0, 0.382, 0.382, 1, 0, 0.382, 0.382, 1, 0]
        if abs(sum(diag) - number_occupied_bands) > tol:
            print('Wrong filling factor for nu %i' % self.nu)
            exit()
        # print('here_seed_asymmetric')
        # import time
        # time.sleep(.1)
        rhorand16 = self.ramdom_16x16_density()

        rho0const = (1 - alpha_rand_asymmetric_calcs) * np.diag(diag) + alpha_rand_asymmetric_calcs * rhorand16

        return rho0const

    def assign_densities(self):
        if self.model_regime == 'full_range':
            self.rho0constUp = remove_small_imag(self.diag_full_regime_high_u(+1))
            self.rho0constUm = remove_small_imag(self.diag_full_regime_high_u(-1))
            self.rho0const_small_u = remove_small_imag(self.diag_full_regime_small_u())
            self.rhoRandom = self.ramdom_16x16_density()
            # self.rho0constUp = remove_small_imag(self.diag_full_regime_v2())
            # self.rho0constUm = remove_small_imag(self.diag_full_regime_v2())
        # elif self.model_regime == 'full_range':
        #     self.rho0constUp = remove_small_imag(self.diag_full_regime_high_u(+1))
        #     self.rho0constUm = remove_small_imag(self.diag_full_regime_high_u(-1))
        #     self.rho0const_small_u = remove_small_imag(self.diag_full_regime_small_u(-1))
        elif self.model_regime == 'no_LL2_mixing_and_asym':
            self.rho0constUp = remove_small_imag(self.seed_asymmetric_calcs())
            self.rho0constUm = remove_small_imag(self.seed_asymmetric_calcs())
            self.rho0const_small_u = remove_small_imag(self.diag_full_regime_small_u())
            self.rhoRandom = self.ramdom_16x16_density()
        # if self.model_regime == 'full_range':
        #     self.rho0constUp = np.real(self.diag_full_regime(+1))
        #     self.rho0constUm = np.real(self.diag_full_regime(-1))
        # elif self.model_regime == 'no_LL2_mixing_and_asym':
        #     self.rho0constUp = np.real(self.seed_asymmetric_calcs())
        #     self.rho0constUm = np.real(self.seed_asymmetric_calcs())


def density_by_model_regime(model_regime):
    # print('here')
    densities = Density_Seed(model_regime, nu)
    # densities = Density_Seed('no_LL2_mixing_and_asym', nu)
    densities.assign_densities()
    rho0constUp = densities.rho0constUp
    rho0constUm = densities.rho0constUm
    rhoRandom = densities.rhoRandom
    rho0const_small_u = densities.rho0const_small_u

    return {'rho0constUp': rho0constUp, 'rho0constUm': rho0constUm, 'rho0const_small_u': rho0const_small_u, 'rhoRandom':rhoRandom}

# print(density_by_model_regime('no_LL2_mixing_and_asym'))
