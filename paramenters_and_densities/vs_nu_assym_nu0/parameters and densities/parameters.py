import sys

import numpy as np

nu_default = 4

if len(sys.argv) == 2:
    nu = int(sys.argv[1])
else:
    nu = nu_default

if abs(nu) > 6:
    print('Error: nu (filling factor) out of range')
    exit()
else:
    print('executing script with nu = %i' % nu)
number_occupied_bands = nu + 8  # number_occupied_bands = 8 is the charge neutrality point

u_critical = 9e-3  # value of u in eV for phase transition based on experiment data with screening factor of 0.26 for nu=0
asym = 1
alpha_Zm = 1  # 0.04227165829987071 # k / alpha_k
alpha_int_H = 1 # 0 for none int calculations on the full_range model
apha_H_asym = 1
potential_asym_layers = 0

alpha_reg = 1
alpha_x = 1
uz = 7e-3
uperp = -2e-3
Orbital_Polarized_seed = False
if Orbital_Polarized_seed:
    uz = 6e-3
    uperp = -3e-3

if (not apha_H_asym) or (abs(nu) > 0):
    uz = 0
    uperp = 0

U0minD = -0.5e-3
U0maxD = 1e-3
dU0D = 0.5e-3

fraction_part = round(U0minD * 1e3 % 1, 2)

add_legend_curve = 0  # if activated k values must be changed accordingly to possible values of u inside first for loop of plot.plot_energies

u_zero = 0.5
u_zero = round(u_zero, 4)

tests_mode = 'on'
save_folder_name = 1  # change this to off when done with tests
# tests_mode = 'off'

itmax_full_range = int(2e4)
alpha_rand_full_range_small_u = 0
alpha_rand_full_range_high_u = 0
same_rhoRandom = 1
alpha_rho = 0.05 * 0  # controls numerical regularization for rho (small memory of rho from previews loop)

alpha_k_nu4 = 0.244
alpha_k_nu0 = 0.26


# alpha_k_dic = {i: alpha_k_nu4 for i in range(-6, 7)}
# alpha_k_dic[0] = alpha_k_nu0
# alpha_k = alpha_k_dic[nu]
# def alpha_k_linear(nu):
#    y = (alpha_k_nu4 - alpha_k_nu0) * nu / 4 + alpha_k_nu0
#    return round(y, 4)
def alpha_k_linear(nu):
    if nu == 0:
        return alpha_k_nu0
    elif nu < 0:
        return alpha_k_nu4
    elif nu > 0:
        return alpha_k_nu4


alpha_k = alpha_k_linear(nu)

file_seed = 0
if file_seed:
    alpha_rand_full_range_small_u = 0.6
    # alpha_rand_full_range_small_u = 0.83

print('uz=%(uz).1fmeV, uperp=%(uperp).1fmeV' % {'uz': uz * 1e3, 'uperp': uperp * 1e3})

B = 13
omega = 35 * np.sqrt(B) / 1e3
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
clight = 299792458
hbar = 6.62607e-34 / (2 * np.pi)
el = 1.602176634e-19
dlayer = 3.35e-10
ep0 = 8.8541878128e-12
alpha_tilda = alpha_k * (dlayer * 1e9)
Lb = np.sqrt(hbar / (el * B))
x = dlayer / Lb * alpha_x
epr = 6
Eh = x / np.sqrt(2 * np.pi)
k = (np.sqrt(np.pi / 2) * el) / (4 * np.pi * ep0 * epr * Lb) * alpha_k
couplings_dict = {'gamma0': gamma0, 'gamma1': gamma1, 'gamma4': gamma4, 'gamma3': gamma3, 'Delta': Delta, 'Delta_td': Delta_td, 'k': k}
model_regime = 'full_range'
effective_H0 = False
projected_four_band_H0 = not effective_H0
######################################################################################################################
######################################### hartree_fock_with_asymmetric_interactions.py #############################
# model_regime = 'no_LL2_mixing_and_asym'
alpha_H_oct_int = 1
itmax_asymmetric_calcs = int(1e4)
alpha_reg_asym_calcs = 1
alpha_rand_asymmetric_calcs = 0  # 0 for Ferro, 0.1 for CAF phase and uperp_meV: -3.2 uz_meV: 14.0
alpha_rho_asymmetric_calcs = 0.0001
replace_LLm2_LL2_low_u = 0  # if false this is effectivelly equivalent to fast_none_interact mode for low u regime
######################################################################################################################
######################################################################################################################
parameters_to_save = {'nu': nu,
                      'number_occupied_bands': number_occupied_bands,
                      'u_critical_meV': u_critical * 1e3,
                      'asym': asym,
                      'alpha_H_oct_int': alpha_H_oct_int,
                      'alpha_int_H': alpha_int_H,
                      'apha_H_asym': apha_H_asym,
                      'alpha_reg_asym_calcs': alpha_reg_asym_calcs,
                      'file_seed': file_seed,
                      'screening': alpha_k,
                      # 'seed_large_u': seed_large_u,
                      'Zm_meV': Zm * 1e3,
                      'x': x,
                      'itmax_full_range': itmax_full_range,
                      'itmax_asymmetric_calcs': itmax_asymmetric_calcs,
                      'alpha_rand_full_range_small_u': alpha_rand_full_range_small_u,
                      'alpha_rand_full_range_high_u': alpha_rand_full_range_high_u,
                      'alpha_rand_asymmetric_calcs': alpha_rand_asymmetric_calcs,
                      'alpha_rho': alpha_rho,
                      'alpha_rho_asymmetric_calcs': alpha_rho_asymmetric_calcs,
                      'model_regime': model_regime,
                      'alpha_reg': alpha_reg,
                      'replace_LLm2_LL2_low_u': replace_LLm2_LL2_low_u,
                      'alpha_int_H': alpha_int_H,
                      'uz_meV': uz * 1e3,
                      'uperp_meV': uperp * 1e3,
                      'u_zero_meV': u_zero,
                      'range_meV': (U0minD * 1e3, U0maxD * 1e3, dU0D * 1e3),
                      'tests_mode': tests_mode,
                      'save_folder_name': save_folder_name,
                      'potential_asym_layers': potential_asym_layers,
                      'couplings_dict': couplings_dict,
                      'effective_H0': effective_H0,
                      'Orbital_Polarized_seed': Orbital_Polarized_seed
                      # 'add_int_to_bands_LLm2_LL2_low_u':add_int_to_bands_LLm2_LL2_low_u
                      }
# Zm = k * 0.0178 #use def above, this is slightly off due to alpha_k
if model_regime == 'full_range':
    parameters_in_plot_list = ['nu', 'asym', 'screening', 'alpha_int_H', 'apha_H_asym',
                               'itmax_full_range', 'alpha_rand_full_range_small_u', 'alpha_rand_full_range_high_u', 'alpha_rho', 'alpha_reg', 'uz_meV', 'uperp_meV', 'effective_H0']
elif model_regime == 'no_LL2_mixing_and_asym':
    parameters_in_plot_list = ['nu', 'asym', 'screening', 'alpha_int_H', 'apha_H_asym', 'alpha_H_oct_int', 'replace_LLm2_LL2_low_u',
                               'itmax_asymmetric_calcs', 'alpha_rand_asymmetric_calcs', 'alpha_rho_asymmetric_calcs', 'alpha_reg_asym_calcs', 'uz_meV', 'uperp_meV']
parameters_to_plot = {key: parameters_to_save[key] for key in parameters_in_plot_list}

parameters_to_plot_text = []
for key in sorted(list(parameters_to_plot.keys()), key=str.lower):
    parameters_to_plot_text.append(str('%s: %s' % (key, parameters_to_plot[key])))

parameters_to_folder_text = []
for key in sorted(list(parameters_to_plot.keys()), key=str.lower):
    parameters_to_folder_text.append(str('%s' % (parameters_to_plot[key])))
