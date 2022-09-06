import pandas as pd
from matplotlib import pyplot as plt

if __name__ == "__main__":
    # setting path
    import sys

    sys.path.append('../')

from config import bands, input_dir_path, dir_path, current_date, tests_mode, results_dir_path_plot_vs_nu
from input.parameters import alpha_tilda, u_zero

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
              'LL2_Km_Sup': {'color': 'black', 'line_shape': '--', 'marker_shape': '^', 'label': '$\\ \\ \\,2\\mathrm{K}^{-}\\uparrow$'},
              'fermi_energy': {'color': 'purple', 'line_shape': '-', 'label': r'Fermi energy'}
              }

plt.rcParams['figure.dpi'] = 150


def plot_energies(energies, nu):
    f = plt.figure()
    ax = plt.gca()

    for band in bands:
        styleX = style_dict.get(band)
        energies.plot(x='u', y=band, color=styleX['color'], style=styleX['line_shape'], markersize=3, linewidth=0.7, label=styleX['label'], ax=ax)  # , marker='o')
    fermi_energy_style = style_dict['fermi_energy']
    energies.plot(x='u', y='fermi_energy', color=fermi_energy_style['color'], style=fermi_energy_style['line_shape'], markersize=3, linewidth=1, label=fermi_energy_style['label'],
                  ax=ax)
    plt.title('Energy bands  nu=' + str(nu) + ' as function of U with self-energy alpha 1')
    plt.xlabel('U(meV)')
    plt.ylabel('Energy bands(meV)')

    # plt.legend(bbox_to_anchor=(1, 0.55))
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    number_occupied_bands_local = nu + 8
    results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + tests_mode

    f.savefig(results_dir_path_local + "LL(U)_HF_interactions_w_SE_warping_alpha1_nu_" + str(nu) + ".pdf", bbox_inches='tight')


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
                                                       'label': '$-2\\mathrm{K}^{-}\\uparrow\\ \\longrightarrow\\ 1\\mathrm{K}^{-}\\uparrow$'},
                         }


def plot_transitions(transitions_df, nu):
    f = plt.figure()
    ax = plt.gca()

    for transition in transitions_df.drop('u', axis=1).columns:
        # print(transition)
        style_transition = transitions_style_dic.get(transition)
        transitions_df.plot(x='u', y=transition, color=style_transition['color'], style=style_transition['line_shape'], markersize=3, linewidth=0.7,
                            label=style_transition['label'], ax=ax)  # , marker='o')
    # print(transitions_df)
    plt.title('Transition nu=' + str(nu) + ' as function of U with self-energy')
    plt.xlabel('U(meV)')
    plt.ylabel('transitions(meV)')
    if nu == 0:
        nu0_exp_transition_energies = pd.read_csv(input_dir_path + 'DGBLG_nu0_data.csv')
        nu0_exp_transition_energies['u'] = nu0_exp_transition_energies['D_mV_per_nm'] * alpha_tilda
        nu0_exp_transition_energies.plot(x='u', y=['peak_A_meV', 'peak_B_meV', 'peak_C_meV'], linestyle='None', label=['experiment', '_nolegend_', '_nolegend_'], color='green',
                                         marker='o', fillstyle='none', ax=ax)  # print(nu0_exp_transition_energies)
        # print('here')
        # print(nu0_exp_transition_energies)

    if nu == 4:
        nu0_exp_transition_energies = pd.read_csv(input_dir_path + 'DGBLG122118_nu4_peak_energies.csv')
        # print(nu0_exp_transition_energies)
        nu0_exp_transition_energies['u'] = nu0_exp_transition_energies['e_nu4'] * alpha_tilda
        # nu0_exp_transition_energies.plot(x='u', y=['peak_A_meV', 'peak_B_meV', 'peak_C_meV'], linestyle='None', label=['experiment', '_nolegend_', '_nolegend_'], color='green',
        #                                  marker='o', fillstyle='none', ax=ax)
        nu0_exp_transition_energies.plot(x='u', y=['peak_nu4_HIGH', 'peak_nu4_LOW'], linestyle='None', label=['experiment', '_nolegend_'], color='green', marker='o',
                                         fillstyle='none', ax=ax)  # print(nu0_exp_transition_energies)

    # plt.legend(bbox_to_anchor=(1, 0.55))
    plt.legend(loc='upper right', bbox_to_anchor=(.99, 0.5))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    number_occupied_bands_local = nu + 8
    results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + tests_mode

    f.savefig(results_dir_path_local + "Transition_nu_" + str(nu) + ".pdf", bbox_inches='tight')


def plot_energies_vs_nu():
    # pass
    energies_df = pd.DataFrame([])
    for nu in range(-6, 7):
        number_occupied_bands_local = nu + 8
        results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + '/'
        tmp = pd.read_csv(results_dir_path_local + 'energies_' + 'nu_' + str(nu) + '.csv')  # ,header=None
        tmp.insert(0, 'nu', nu)
        tmp = tmp[round(tmp['u'], 4) == u_zero]
        if tmp.empty:
            print('u_zero=' + str(u_zero) + 'meV not found for energies at nu=' + str(nu))
        # tmp.drop('u',axis=1,inplace=True)
        energies_df = pd.concat([energies_df, tmp], ignore_index=True)
    energies_df.to_csv(results_dir_path_plot_vs_nu + 'energies_vs_nu.csv', index=False)
    print(energies_df)

    f = plt.figure()
    font = {'size': 15}
    plt.rc('font', **font)
    ax = plt.gca()
    plt.rcParams['figure.dpi'] = 150
    # transitions.plot(x='nu', linestyle=':', linewidth=1, markersize=12, ax=ax)
    # for i, line in enumerate(ax.get_lines()):
    #     line.set_marker(markers_list[i])
    #     line.set_label(label_list[i])
    #     line.set_color(color_list[i])
    for band in energies_df.drop(['nu', 'u', 'fermi_energy'], axis=1).columns:
        # print(band)
        style_band = style_dict.get(band)
        energies_df.plot(x='nu', y=band, color=style_band['color'], linestyle=style_band['line_shape'], linewidth=1, label=style_band['label'], ax=ax)  # , marker='o')
    fermi_energy_style = style_dict['fermi_energy']
    energies_df.plot(x='nu', y='fermi_energy', color=fermi_energy_style['color'], style=fermi_energy_style['line_shape'], markersize=3, linewidth=1,
                     label=fermi_energy_style['label'], ax=ax)
    # transitions_exp = pd.read_csv(input_dir_path + 'DGBLG122118_D0_peak_energies_02232022.csv')
    # transitions_exp.plot(x='nu_D0', y=['peak_High_meV_D0', 'peak_Mid_meV_D0', 'peak_Low_meV_D0', 'peak_LowPrime_meV_D0'], linestyle='--', linewidth=1,
    #                      label=['experiment', '_nolegend_', '_nolegend_', '_nolegend_'], color='tab:green', marker='o', fillstyle='none', ax=ax)

    plt.xticks(range(-6, 5, 1))
    plt.legend(bbox_to_anchor=(0.8, 0.7))
    plt.rcParams["figure.figsize"] = (10, 5)
    plt.title('Energies as function of filling factor for U ~ 0meV.', fontsize=19, y=-0.24, x=0.55)
    plt.xlabel(r'$\nu$')
    plt.ylabel('Energies(meV)')
    f.savefig(results_dir_path_plot_vs_nu + 'energies_vs_nu.pdf', bbox_inches='tight')


def plot_transitions_vs_nu():
    # pass

    transitions_df = pd.DataFrame([])
    for nu in range(-6, 7):
        number_occupied_bands_local = nu + 8
        results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + '/'
        tmp = pd.read_csv(results_dir_path_local + 'transitions_' + 'nu_' + str(nu) + '.csv')  # ,header=None
        tmp.insert(0, 'nu', nu)
        tmp = tmp[round(tmp['u'], 4) == u_zero]
        if tmp.empty:
            print('u_zero=' + str(u_zero) + 'meV not found for transitions at nu=' + str(nu))
        # tmp.drop('u',axis=1,inplace=True)
        transitions_df = pd.concat([transitions_df, tmp], ignore_index=True)
    transitions_df.to_csv(results_dir_path_plot_vs_nu + 'transitions_vs_nu.csv', index=False)
    # print(transitions_df)

    f = plt.figure()
    font = {'size': 15}
    plt.rc('font', **font)
    ax = plt.gca()
    plt.rcParams['figure.dpi'] = 150
    # transitions.plot(x='nu', linestyle=':', linewidth=1, markersize=12, ax=ax)
    # for i, line in enumerate(ax.get_lines()):
    #     line.set_marker(markers_list[i])
    #     line.set_label(label_list[i])
    #     line.set_color(color_list[i])
    for transition in transitions_df.drop(['nu', 'u'], axis=1).columns:
        style_transition = transitions_style_dic.get(transition)
        transitions_df.plot(x='nu', y=transition, color=style_transition['color'], linestyle=':', linewidth=1, markersize=12, marker=style_transition['marker_shape'],
                            # style=style_transition['line_shape'], # markersize=3, linewidth=0.7,
                            label=style_transition['label'], ax=ax)  # , marker='o')

    transitions_exp = pd.read_csv(input_dir_path + 'DGBLG122118_D0_peak_energies_02232022.csv')
    transitions_exp.plot(x='nu_D0', y=['peak_High_meV_D0', 'peak_Mid_meV_D0', 'peak_Low_meV_D0', 'peak_LowPrime_meV_D0'], linestyle='--', linewidth=1,
                         label=['experiment', '_nolegend_', '_nolegend_', '_nolegend_'], color='tab:green', marker='o', fillstyle='none', ax=ax)

    plt.xticks(range(-6, 5, 1))
    plt.legend(bbox_to_anchor=(0.8, 0.7))
    plt.rcParams["figure.figsize"] = (10, 5)
    plt.title('Transition energies as function of filling factor for U ~ 0meV.', fontsize=19, y=-0.24, x=0.55)
    plt.xlabel(r'$\nu$')
    plt.ylabel('Transitions(meV)')
    f.savefig(results_dir_path_plot_vs_nu + 'transitions_vs_nu.pdf', bbox_inches='tight')


#
# def plot_energies_with_asymmetry(nu):
#     folder_name_CAF7 = 'files_asym_1__itmax_10000__Zm_0.753__alpha_H_oct_int_1__uz_7.0__uperp_-1.6__x_0.047__alpha_state_1__alpha_rand_0.01__dens_3'
#     # folder_name_CAF6 = 'files_asym_1__itmax_1000.0__Zm_0.753__alpha_H_oct_int_1__uz_6.0__uperp_-1.6__x_0.047__alpha_state_1__alpha_rand_0.01__dens_3'#CAF for U<12meV ocupation 0.62:0.38 $$$$$$$$$$
#     folder_name_CAF_uz6_uperpm1 = 'files_asym_1__itmax_5000__Zm_0.753__alpha_H_oct_int_1__uz_6.0__uperp_-1.0__x_0.047__alpha_state_1__alpha_rand_0.01__dens_2'  # CAF for U<12meV ocupation 0.69:0.31 $$$$$$$$$$
#     folder_name_CAF5 = 'files_asym_1__itmax_1000__Zm_0.753__alpha_H_oct_int_1__uz_5.0__uperp_-1.6__x_0.047__alpha_state_1__alpha_rand_0.01__dens_3'
#     folder_name_CAF6 = 'files_asym_1__itmax_1000.0__Zm_0.753__alpha_H_oct_int_1__uz_6.0__uperp_-1.6__x_0.047__alpha_state_1__alpha_rand_0.01__dens_3'  # CAF for U<12meV ocupation 0.62:0.38 $$$$$$$$$$
#     folder_name_CAF8 = 'files_asym_1__itmax_80000__Zm_0.753__alpha_H_oct_int_1__uz_8.0__uperp_-1.6__x_0.047__alpha_state_1__alpha_rand_0.01__dens_3'
#     folder_name_CAF9 = 'files_asym_1__itmax_10000__Zm_0.753__alpha_H_oct_int_1__uz_9.0__uperp_-1.6__x_0.047__alpha_state_1__alpha_rand_0.01__dens_3'
#     folder_name_with_asymmetry = folder_name_CAF7
#     # folder_name_with_asymmetry = folder_name_CAF_uz6_uperpm1
#     folder_name_file_name = dir_path + '/canted_antiferromagnetic/CAF_oct_interaction/' + folder_name_with_asymmetry + '/eigenU0_fullH0_Delta_CAF_tests.csv'
#     energies_df = pd.read_csv(folder_name_file_name)
#     energies_df.columns = ['u'] + bands
#     # print(energies_df[energies_df['u']==2.0])
#
#     print(energies_df.iloc[[2, 4]])
#
#     energies_df, transition_energy_df = transitions_energy_fermi_energy(energies_df, nu)
#     # if nu==0:
#     #     for column_name in energies_df.columns:
#     #         if 'LLm2' in column_name:
#     #             energies_df[column_name] = energies_df[column_name]
#
#     number_occupied_bands_local = nu + 8
#     results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + '/'
#     print(results_dir_path_local)
#     energies_df.to_csv(results_dir_path_local + 'energies_' + 'nu_' + str(nu) + '.csv', index=False)
#     transition_energy_df.to_csv(results_dir_path_local + 'transitions_' + 'nu_' + str(nu) + '.csv', index=False)
#     plot_energies(energies_df, nu)
#     plot_transitions(transition_energy_df, nu)
#     # # eigenU = pd.read_csv(folder_name+namecsv,header=None).applymap(complex).applymap( lambda x: round(x.real,decimals) )
#     # # eigenU.columns = ['u'] + bands
#     #
#     # # display(eigenU)
#     #
#     # style_dict = {'0p-': ('lightblue', '-', 'v', r'$\ \ \,0\mathrm{K}^{+}\downarrow$'),
#     #               '1p-': ('salmon', '-', 'v', r'$\ \ \,1\mathrm{K}^{+}\downarrow$'),
#     #               '-2p-': ('gray', '-', 'v', r'$-2\mathrm{K}^{+}\downarrow$'),
#     #               '2p-': ('gray', '-', 'v', r'$\ \ \,2\mathrm{K}^{+}\downarrow$'),
#     #               '0m-': ('lightblue', '--', 'v', r'$\ \ \,0\mathrm{K}^{-}\downarrow$'),
#     #               '1m-': ('salmon', '--', 'v', r'$\ \ \,1\mathrm{K}^{-}\downarrow$'),
#     #               '-2m-': ('gray', '--', 'v', r'$-2\mathrm{K}^{-}\downarrow$'),
#     #               '2m-': ('gray', '--', 'v', r'$\ \ \,2\mathrm{K}^{-}\downarrow$'),
#     #               '0p+': ('blue', '-', '^', r'$\ \ \,0\mathrm{K}^{+}\uparrow$'),
#     #               '1p+': ('red', '-', '^', r'$\ \ \,1\mathrm{K}^{+}\uparrow$'),
#     #               '-2p+': ('black', '-', '^', r'$-2\mathrm{K}^{+}\uparrow$'),
#     #               '2p+': ('black', '-', '^', r'$\ \ \,2\mathrm{K}^{+}\uparrow$'),
#     #               '0m+': ('blue', '--', '^', r'$\ \ \,0\mathrm{K}^{-}\uparrow$'),
#     #               '1m+': ('red', '--', '^', r'$\ \ \,1\mathrm{K}^{-}\uparrow$'),
#     #               '-2m+': ('black', '--', '^', r'$-2\mathrm{K}^{-}\uparrow$'),
#     #               '2m+': ('black', '--', '^', r'$\ \ \,2\mathrm{K}^{-}\uparrow$')}
#     # #             ,'fermi_energy':('purple','-','.')}
#     # bands = ['0p-', '1p-', '-2p-', '2p-', '0m-', '1m-', '-2m-', '2m-', '0p+', '1p+', '-2p+', '2p+', '0m+', '1m+', '-2m+', '2m+']
#     #
#     # f = plt.figure()
#     # ax = plt.gca()
#     #
#     # plt.rcParams['figure.dpi'] = 150
#     #
#     # for x in bands:
#     #     styleX = style_dict.get(x)
#     #     eigenU.plot(x='u', y=x, color=styleX[0], style=styleX[1], markersize=3, linewidth=0.7, label=styleX[3], ax=ax)  # , marker='o')
#     # #     eigenU.plot(x='u', y=x, color=styleX[0], style=styleX[1],marker='.', markersize=.5,markerfacecolor='Black', linestyle = 'None', label=styleX[3], ax=ax)#, marker='o')
#     # # fermi_energy.plot(x='u', y='fermi_energy', color='purple', style='-', markersize=3, linewidth=1, label=r'Fermi energy', ax=ax)


if __name__ == "__main__":
    # energies_df = pd.read_csv(results_dir_path + 'energies_' + file_name_csv)
    # transitions_df = pd.read_csv(results_dir_path + 'transitions_' + file_name_csv)
    # print(energies_df)
    # plot_energies(energies_df)
    # # plot_transitions(transitions_df)

    # plot_energies_with_asymmetry(nu=0)
    plot_energies_vs_nu()
    plot_transitions_vs_nu()
