import sys

import numpy as np

nu_default = 0

if len(sys.argv) == 2:
    nu = int(sys.argv[1])
else:
    nu = nu_default

if abs(nu) > 6:
    print('Error: nu (filling factor) out of range')
    exit()
else:
    print('executing script with nu = %i' % nu)

u_critical = 12e-3  # value of u in eV for phase transition
asym = 1
alpha_Zm = 1  # 0.04227165829987071 # k / alpha_k
alpha_H_oct_int = 1
alpha_int_H = 1  # 0 for none int calculations on the full_range model
apha_H_asym = 0
alpha_reg = 1
alpha_reg_asym_calcs = 1
alpha_x = 1
uz = 7e0 * 1e-3
uperp = -1.6e0 * 1e-3

U0minD = -10e-3
U0maxD = 25e-3
dU0D = 1e-3

u_zero = 1
u_zero = round(u_zero, 4)

tests_mode = 'on'  # change this to off when done with tests
# tests_mode = 'off'

itmax_full_range = int(1e5)
itmax_asymmetric_calcs = int(500)
if nu == 0:
    alpha_rand_full_range = 0.6  # 0.6
else:
    alpha_rand_full_range = 0
alpha_rand_asymmetric_calcs = 0.005
alpha_rho = 0.005  # controls numerical regularization for rho (memory of rho from previews loop)
# variables_dict['asym']=asym
# variables_dict['alpha_Zm']=alpha_Zm
# variables_dict['alpha_H_oct_int']=alpha_H_oct_int
# variables_dict['alpha_x']=alpha_x

alpha_k_nu4 = 0.244
alpha_k_nu0 = 0.260
alpha_k_dic = {i: alpha_k_nu4 for i in range(-6, 7)}
alpha_k_dic[0] = alpha_k_nu0

# variables_dict['alpha_k_dic']=alpha_k_dic
# print(alpha_k_dic)
# alpha_k_dic = {-6: alpha_k_nu4,
#                -4: alpha_k_nu4,
#                -3: alpha_k_nu4,
#                -2: alpha_k_nu4,
#                -1: alpha_k_nu4,
#                0: alpha_k_nu0,
#                1: alpha_k_nu4,
#                2: alpha_k_nu4,
#                3: alpha_k_nu4,
#                4: alpha_k_nu4,
#                5: alpha_k_nu4,
#                6: alpha_k_nu4
#                }

# asymmetry parameters
# model_regime = 'near_zero_dielectric_field'
model_regime = 'full_range'

#### appoximation mode for LL2 and LLm2 ###
mode = 'hartree_fock_and_regularization_calcs'
add_int_to_bands_LLm2_LL2_low_u = True  # if false this is effectivelly equivalent to fast_none_interact mode for low u regime
# mode = 'fast_none_interact'
# mode = 'fast_from_file' # L2 and LLm2 are set by fully interacting hamiltonian (no asymmetric int)
# mode = 'fast_from_constant' # LL2 and LLm2 are set by none-interacting hamiltonian + constant (from avereage LL2 and LLm2 set by interacting hamiltonian)
#########
# if tests_mode=='off':
#     mode = 'hartree_fock_and_regularization_calcs'
# alpha_state = 1

# dens = 3
print('uz=%(uz).1fmeV, uperp=%(uperp).1fmeV' % {'uz': uz * 1e3, 'uperp': uperp * 1e3})
# variables_dict['alpha_rho']=alpha_rho
# variables_dict['alpha_rand']=alpha_rand
# variables_dict['alpha_reg']=alpha_reg


# variables_dict['uz']=uz
# variables_dict['uperp']=uperp
# variables_dict['uz']=uz

B = 13
omega = 35 * np.sqrt(B) / 1000
gamma0 = 3
gamma1 = 0.41
gamma4 = 0.15 * asym
eta4 = gamma4 / gamma0
gamma3 = -0.3 * asym  # gamma3 = 0.3;
eta3 = gamma3 / gamma0
Delta = 0.018 * asym
Delta_ab = 0.001 * 0
Delta_td = (2 * gamma1 * gamma4 / gamma0 + Delta)
beta = (omega / gamma1) ** 2
Zm = 57.9e-6 * B * alpha_Zm  # in eV,  57.9e-6 is from Zhang2011PRB below equation 17 #temperature: 4 K=4 * 0.0862=0.3448meV # k * 0.0178 #
# print(Zm)
# ep = 0.5
number_occupied_bands = nu + 8  # nu = 8 is the charge neutrality point in this schema
alpha_k = alpha_k_dic[nu]

clight = 299792458
hbar = (6.62607 * 10 ** (-34)) / (2 * np.pi)
el = 1.602176634 * 10 ** (-19)
dlayer = 3.35 * 10 ** (-10)
ep0 = 8.8541878128 * 10 ** (-12)

alpha_tilda = alpha_k * (dlayer * 10 ** 9)

Lb = np.sqrt(hbar / (el * B))
x = dlayer / Lb * alpha_x
epr = 6
Eh = x / np.sqrt(2 * np.pi)
k = (np.sqrt(np.pi / 2) * el) / (4 * np.pi * ep0 * epr * Lb) * alpha_k

# Zm = k * 0.0178 #use def above, this is slightly off due to alpha_k
parameters_to_plot = {'nu': nu,
                      # 'number_occupied_bands': number_occupied_bands,
                      # 'u_critical_meV': u_critical * 1e3,
                      'asym': asym,
                      'alpha_H_oct_int': alpha_H_oct_int,
                      'apha_H_asym': apha_H_asym,
                      'alpha_reg_asym_calcs': alpha_reg_asym_calcs,
                      # 'Zm_meV': Zm * 1e3,
                      # 'x': x,
                      # 'itmax_full_range': itmax_full_range,
                      'itmax_asymmetric_calcs': itmax_asymmetric_calcs,
                      # 'alpha_rand_full_range': alpha_rand_full_range,
                      # 'alpha_rand_asymmetric_calcs': alpha_rand_asymmetric_calcs,
                      # 'alpha_rho': alpha_rho,
                      # 'alpha_k': alpha_k_dic[nu],
                      'model_regime': model_regime,
                      'mode': mode,
                      'alpha_reg': alpha_reg,
                      'replace_LLm2_LL2_low_u': add_int_to_bands_LLm2_LL2_low_u,
                      'alpha_int_H': alpha_int_H,
                      'uz_meV': uz * 1e3,
                      'uperp_meV': uperp * 1e3,
                      # 'u_zero_meV': u_zero,
                      # 'range_meV': (U0minD * 1e3, U0maxD * 1e3, dU0D * 1e3),
                      # 'tests_mode': tests_mode,
                      # 'add_int_to_bands_LLm2_LL2_low_u':add_int_to_bands_LLm2_LL2_low_u
                      }

parameters_to_plot_text = []
for key in sorted(list(parameters_to_plot.keys()), key=str.lower):
    parameters_to_plot_text.append(str('%s: %s' % (key, parameters_to_plot[key])))

parameters_to_folder_text = []
for key in sorted(list(parameters_to_plot.keys()), key=str.lower):
    parameters_to_folder_text.append(str('%s' % (parameters_to_plot[key])))

parameters_to_save = {'nu': nu,
                      'number_occupied_bands': number_occupied_bands,
                      'u_critical_meV': u_critical * 1e3,
                      'asym': asym,
                      'alpha_H_oct_int': alpha_H_oct_int,
                      'apha_H_asym': apha_H_asym,
                      'alpha_reg_asym_calcs': alpha_reg_asym_calcs,
                      'Zm_meV': Zm * 1e3,
                      'x': x,
                      'itmax_full_range': itmax_full_range,
                      'itmax_asymmetric_calcs': itmax_asymmetric_calcs,
                      'alpha_rand_full_range': alpha_rand_full_range,
                      'alpha_rand_asymmetric_calcs': alpha_rand_asymmetric_calcs,
                      'alpha_rho': alpha_rho,
                      'alpha_k': alpha_k_dic[nu],
                      'model_regime': model_regime,
                      'mode': mode,
                      'alpha_reg': alpha_reg,
                      'replace_LLm2_LL2_low_u': add_int_to_bands_LLm2_LL2_low_u,
                      'alpha_int_H': alpha_int_H,
                      'uz_meV': uz * 1e3,
                      'uperp_meV': uperp * 1e3,
                      # 'u_zero_meV': u_zero,
                      'range_meV': (U0minD * 1e3, U0maxD * 1e3, dU0D * 1e3),
                      'tests_mode': tests_mode,
                      # 'add_int_to_bands_LLm2_LL2_low_u':add_int_to_bands_LLm2_LL2_low_u
                      }

# if __name__ == "__main__":
#     # setting path
#     import sys
#     sys.path.append('../')
#     from config import results_dir_path
# print(sorted(parameters_to_save, key=parameters_to_save.get))
# print(sorted(parameters_to_save, key=lambda key: parameters_to_save[key]))
# with open("myfile.txt", 'w') as f:
#     for key, value in details.items():
#         f.write('%s:%s\n' % (key, value))
# print((  ))

# from datetime import datetime
# with open('parameters.txt', 'w') as parameters_file:
#     # for key in sorted(parameters_to_save.keys()):
#     parameters_file.write('created on (Y-m-d): %s \n' % datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
#     for key in sorted(list(parameters_to_save.keys()),key=str.lower):
#         parameters_file.write('%s = %s\n' % (key, parameters_to_save[key]))

# print(x)
# print(Zm, alpha_H_oct_int, uz, uperp, x)
# variables_dict = {}
# variables_dict.update({k:v for k,v in locals().copy().iteritems() if k[:2] != '__' and k != 'variables_dict'})
# variables_dict = {k:v for k,v in locals().copy().iteritems() if k[:2] != '__' and k != 'variables_dict'}
# variables_dict = {str(i): locals()[i] for i in locals()}
# print(variables_dict)
