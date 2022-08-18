import sys

import pandas as pd
from matplotlib import pyplot as plt

if __name__ == "__main__":
    # setting path
    sys.path.append('../')

from config import aux_dir_path, file_name_csv, bands, input_dir_path
from input.parameters import nu, alpha_tilda

style_dict = {'LL0_Kp_Sdown': {'color': 'lightblue', 'line_shape': '-', 'marker_shape': 'v', 'label': '$\\ \\ \\,0\\mathrm{K}^{+}\\downarrow$'},
              'LL1_Kp_Sdown': {'color': 'salmon', 'line_shape': '-', 'marker_shape': 'v', 'label': '$\\ \\ \\,1\\mathrm{K}^{+}\\downarrow$'},
              'LLm2_Kp_Sdown': {'color': 'gray', 'line_shape': '-', 'marker_shape': 'v', 'label': '$-2\\mathrm{K}^{+}\\downarrow$'},
              'LL2_Kp_Sdown': {'color': 'gray', 'line_shape': '-', 'marker_shape': 'v', 'label': '$\\ \\ \\,2\\mathrm{K}^{+}\\downarrow$'},
              'LL0_Km_Sdown': {'color': 'lightblue', 'line_shape': '--', 'marker_shape': 'v', 'label': '$\\ \\ \\,0\\mathrm{K}^{-}\\downarrow$'},
              'LL1_Km_Sdown': {'color': 'salmon', 'line_shape': '--', 'marker_shape': 'v', 'label': '$\\ \\ \\,1\\mathrm{K}^{-}\\downarrow$'},
              'LLm2_Km_Sdown': {'color': 'gray', 'line_shape': '--', 'marker_shape': 'v', 'label': '$-2\\mathrm{K}^{-}\\downarrow$'},
              'LL2_Km_Sdown': {'color': 'gray', 'line_shape': '--', 'marker_shape': 'v', 'label': '$\\ \\ \\,2\\mathrm{K}^{-}\\downarrow$'},
              'LL0_Kp_Sup': {'color': 'blue', 'line_shape': '-', 'marker_shape': '^', 'label': '$\\ \\ \\,0\\mathrm{K}^{+}\\uparrow$'},
              'LL1_Kp_Sup': {'color': 'red', 'line_shape': '-', 'marker_shape': '^', 'label': '$\\ \\ \\,1\\mathrm{K}^{+}\\uparrow$'},
              'LLm2_Kp_Sup': {'color': 'black', 'line_shape': '-', 'marker_shape': '^', 'label': '$-2\\mathrm{K}^{+}\\uparrow$'},
              'LL2_Kp_Sup': {'color': 'black', 'line_shape': '-', 'marker_shape': '^', 'label': '$\\ \\ \\,2\\mathrm{K}^{+}\\uparrow$'},
              'LL0_Km_Sup': {'color': 'blue', 'line_shape': '--', 'marker_shape': '^', 'label': '$\\ \\ \\,0\\mathrm{K}^{-}\\uparrow$'},
              'LL1_Km_Sup': {'color': 'red', 'line_shape': '--', 'marker_shape': '^', 'label': '$\\ \\ \\,1\\mathrm{K}^{-}\\uparrow$'},
              'LLm2_Km_Sup': {'color': 'black', 'line_shape': '--', 'marker_shape': '^', 'label': '$-2\\mathrm{K}^{-}\\uparrow$'},
              'LL2_Km_Sup': {'color': 'black', 'line_shape': '--', 'marker_shape': '^', 'label': '$\\ \\ \\,2\\mathrm{K}^{-}\\uparrow$'}}

plt.rcParams['figure.dpi'] = 150


def plot_energies(energies):
    f = plt.figure()
    ax = plt.gca()

    for band in bands:
        styleX = style_dict.get(band)
        energies.plot(x='u', y=band, color=styleX['color'], style=styleX['line_shape'], markersize=3, linewidth=0.7, label=styleX['label'], ax=ax)  # , marker='o')
    energies.plot(x='u', y='fermi_energy', color='purple', style='-', markersize=3, linewidth=1, label=r'Fermi energy', ax=ax)
    plt.title('Energy bands  nu=' + str(nu) + ' as function of U with self-energy alpha 1')
    plt.xlabel('U(meV)')
    plt.ylabel('Energy bands(meV)')

    # plt.legend(bbox_to_anchor=(1, 0.55))
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    f.savefig(aux_dir_path + "LL(U)_HF_interactions_w_SE_warping_alpha1_nu_" + str(nu) + ".pdf", bbox_inches='tight')




transitions_style_dic = {'LL1_Kp_Sdown_to_LL2_Kp_Sdown': {'color': 'tab:blue', 'marker_shape': '2', 'line_shape': '-',
                                                          'label': '$1\\mathrm{K}^{+}\\downarrow\\ \\longrightarrow\\ 2\\mathrm{K}^{+}\\downarrow$'},
                         'LLm2_Kp_Sdown_to_LL1_Kp_Sdown': {'color': 'tab:orange', 'marker_shape': '1', 'line_shape': '-',
                                                           'label': '$-2\\mathrm{K}^{+}\\downarrow\\ \\longrightarrow\\ 1\\mathrm{K}^{+}\\downarrow$'},
                         'LL1_Km_Sdown_to_LL2_Km_Sdown': {'color': 'tab:red', 'marker_shape': '1', 'line_shape': '--',
                                                          'label': '$1\\mathrm{K}^{-}\\downarrow\\ \\longrightarrow\\ 2\\mathrm{K}^{-}\\downarrow$'},
                         'LLm2_Km_Sdown_to_LL1_Km_Sdown': {'color': 'tab:purple', 'marker_shape': '2', 'line_shape': '--',
                                                           'label': '$-2\\mathrm{K}^{-}\\downarrow\\ \\longrightarrow\\ 1\\mathrm{K}^{-}\\downarrow$'},
                         'LL1_Kp_Sup_to_LL2_Kp_Sup': {'color': 'tab:brown', 'marker_shape': '1', 'line_shape': '-',
                                                      'label': '$1\\mathrm{K}^{+}\\uparrow\\ \\longrightarrow\\ 2\\mathrm{K}^{+}\\uparrow$'},
                         'LLm2_Kp_Sup_to_LL1_Kp_Sup': {'color': 'tab:gray', 'marker_shape': '1', 'line_shape': '-',
                                                       'label': '$-2\\mathrm{K}^{+}\\uparrow\\ \\longrightarrow\\ 1\\mathrm{K}^{+}\\uparrow$'},
                         'LL1_Km_Sup_to_LL2_Km_Sup': {'color': 'tab:olive', 'marker_shape': '2', 'line_shape': '--',
                                                      'label': '$1\\mathrm{K}^{-}\\uparrow\\ \\longrightarrow\\ 2\\mathrm{K}^{-}\\uparrow$'},
                         'LLm2_Km_Sup_to_LL1_Km_Sup': {'color': 'tab:cyan', 'marker_shape': '2', 'line_shape': '--',
                                                       'label': '$-2\\mathrm{K}^{-}\\uparrow\\ \\longrightarrow\\ 1\\mathrm{K}^{-}\\uparrow$'}}


def plot_transitions(transitions_df):
    f = plt.figure()
    ax = plt.gca()

    for transition in transitions_df.drop('u', axis=1).columns:
        style_transition = transitions_style_dic.get(transition)
        transitions_df.plot(x='u', y=transition, color=style_transition['color'], style=style_transition['line_shape'], markersize=3, linewidth=0.7,
                            label=style_transition['label'], ax=ax)  # , marker='o')
    plt.title('Transition nu=' + str(nu) + ' as function of U with self-energy')
    plt.xlabel('U(meV)')
    plt.ylabel('transitions(meV)')
    if nu == 0:
        nu0_exp_transition_energies = pd.read_csv(input_dir_path + 'DGBLG_nu0_data.dat', sep='\s\s+|\t', engine='python')
        nu0_exp_transition_energies['U'] = nu0_exp_transition_energies['D_mV_per_nm'] * alpha_tilda
        nu0_exp_transition_energies.plot(x='U', y=['peak_A_meV', 'peak_B_meV', 'peak_C_meV'], linestyle='None', label=['experiment', '_nolegend_', '_nolegend_'], color='green',
                                         marker='o', fillstyle='none', ax=ax)
        # print(nu0_exp_transition_energies)

    if nu == 4:
        nu0_exp_transition_energies = pd.read_csv(input_dir_path + 'DGBLG122118_nu4_peak_energies.csv')
        # print(nu0_exp_transition_energies)
        nu0_exp_transition_energies['U'] = nu0_exp_transition_energies['e_nu4'] * alpha_tilda
        # nu0_exp_transition_energies.plot(x='U', y=['peak_A_meV', 'peak_B_meV', 'peak_C_meV'], linestyle='None', label=['experiment', '_nolegend_', '_nolegend_'], color='green',
        #                                  marker='o', fillstyle='none', ax=ax)
        nu0_exp_transition_energies.plot(x='U', y=['peak_nu4_HIGH', 'peak_nu4_LOW'], linestyle='None', label=['experiment', '_nolegend_'], color='green', marker='o',
                                         fillstyle='none', ax=ax)
        # print(nu0_exp_transition_energies)

    # plt.legend(bbox_to_anchor=(1, 0.55))
    plt.legend(loc='upper right', bbox_to_anchor=(.55, 1.020))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    f.savefig(aux_dir_path + "Transition_nu_" + str(nu) + ".pdf", bbox_inches='tight')


if __name__ == "__main__":
    energies_df = pd.read_csv(aux_dir_path + 'energies_' + file_name_csv)
    transitions_df = pd.read_csv(aux_dir_path + 'trans_' + file_name_csv)
    # print(energies_df)
    plot_energies(energies_df)
    plot_transitions(transitions_df)
