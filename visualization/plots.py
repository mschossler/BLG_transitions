import pandas as pd
from matplotlib import pyplot as plt

# import warnings
# warnings.filterwarnings( "ignore", module = "pandas/plotting/_matplotlib\..*" )
if __name__ == "__main__":
    # setting path
    # import numpy as np
    import sys

    sys.path.append('../')

from config import bands, base_octet, input_dir_path, dir_path, current_date, tests_mode, results_dir_path_plot_vs_nu, current_time, results_dir_path
from input.parameters import alpha_tilda, u_zero, parameters_to_plot_text, add_legend_curve, nu

if __name__ == '__main__':
    import os

    os.remove(results_dir_path + '/progress.txt')
    os.rmdir(results_dir_path)
    # current_date = '10102022'

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

    i = 0
    for band in base_octet:
        styleX = style_dict.get(band)
        energies.plot(x='u', y=band, color=styleX['color'], style=styleX['line_shape'], markersize=3, linewidth=0.7, label=styleX['label'], ax=ax)  # , marker='o')
        if add_legend_curve:
            if '1_Kp' in band: k, rotation = 20, 15
            if '1_Km' in band: k, rotation = 50, -15
            if '0_Km' in band: k, rotation = 30, -15
            if '0_Kp' in band: k, rotation = 30, 15
            if 'Sdown' in band: xytext, ha, va = (energies[energies['u'] == k]['u'], energies[energies['u'] == k][band] + 1), 'center', 'top'
            if 'Sup' in band: xytext, ha, va = (energies[energies['u'] == k]['u'], energies[energies['u'] == k][band] - 2), 'center', 'bottom'
            # ax.annotate(styleX['label'], color=styleX['color'], xy=(energies[energies['u'] == k]['u'], energies[energies['u'] == k][band]), xytext=xytext, arrowprops=dict(
            #     arrowstyle="->"))
            plt.annotate(styleX['label'], color=styleX['color'], xy=xytext, fontsize=7, ha=ha, va=va, rotation=rotation)
            i += 1

    for band in [band for band in bands if '2' in band]:
        styleX = style_dict.get(band)
        energies.plot(x='u', y=band, color=styleX['color'], style=styleX['line_shape'], markersize=3, linewidth=0.7, label=styleX['label'], ax=ax)  # , marker='o')
        # xy = (energies[energies['u']==15]['u'], energies[energies['u']==15][band])
        # ax.annotate("test", xy=xy, xytext=(15, 0), arrowprops=dict(arrowstyle="->"))
    fermi_energy_style = style_dict['fermi_energy']
    energies.plot(x='u', y='fermi_energy', color=fermi_energy_style['color'], style=fermi_energy_style['line_shape'], markersize=3, linewidth=1, label=fermi_energy_style['label'],
                  ax=ax)
    plt.title('Energy bands  nu=' + str(nu) + ' as function of U')
    textstr = '\n'.join(parameters_to_plot_text)
    # print(textstr)
    plt.text(4, -40, textstr, fontsize=4, verticalalignment='top')

    plt.xlabel('U(meV)')
    plt.ylabel('Energy bands(meV)')

    # plt.legend(bbox_to_anchor=(1, 0.55))
    # anchored_text=AnchoredText(textstr, loc='best')
    # ax.add_artist(anchored_text)
    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    number_occupied_bands_local = nu + 8
    results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + tests_mode

    f.savefig(results_dir_path_local + "LL(U)_nu_" + str(nu) + '_' + current_time + ".pdf", bbox_inches='tight')


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
    plt.title('Transition nu=' + str(nu) + ' as function of U')

    textstr = '\n'.join(parameters_to_plot_text)
    # print(textstr)
    plt.text(22, 68, textstr, fontsize=5, verticalalignment='top')
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
        nu4_exp_transition_energies = pd.read_csv(input_dir_path + 'DGBLG122118_nu4_peak_energies.csv')
        # print(nu0_exp_transition_energies)
        nu4_exp_transition_energies['u'] = nu4_exp_transition_energies['e_nu4'] * alpha_tilda
        # nu0_exp_transition_energies.plot(x='u', y=['peak_A_meV', 'peak_B_meV', 'peak_C_meV'], linestyle='None', label=['experiment', '_nolegend_', '_nolegend_'], color='green',
        #                                  marker='o', fillstyle='none', ax=ax)

        nu4_exp_transition_energies.plot(x='u', y=['peak_nu4_HIGH', 'peak_nu4_LOW'], linestyle='None', label=['experiment', '_nolegend_'], color='green', marker='o',
                                         fillstyle='none', ax=ax)

    # print(nu0_exp_transition_energies)

    # plt.legend(bbox_to_anchor=(1, 0.55))

    # plt.show()
    plt.legend(loc='upper right', bbox_to_anchor=(.99, 0.5))
    plt.rcParams["figure.figsize"] = (10, 5)
    number_occupied_bands_local = nu + 8
    results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local) + tests_mode

    f.savefig(results_dir_path_local + "Transition_nu_" + str(nu) + '_' + current_time + ".pdf", bbox_inches='tight')


def plot_energies_vs_nu():
    # pass
    energies_df = pd.DataFrame([])
    folder_names = []
    for nu in range(-6, 7):
        number_occupied_bands_local = nu + 8
        results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local)

        try:
            folder_names_file = open(results_dir_path_local + '/folder_list.txt', 'r').read()
        except FileNotFoundError:
            print('add file folders to folder_list')
        else:
            folder_name = folder_names_file.splitlines()[-1]
            folder_names.append(folder_name)
        tmp = pd.read_csv(results_dir_path_local + folder_name + 'energies_' + 'nu_' + str(nu) + '.csv')  # ,header=None
        tmp.insert(0, 'nu', nu)
        tmp = tmp[round(tmp['u'], 4) == u_zero]
        if tmp.empty:
            print('u_zero=' + str(u_zero) + 'meV not found for energies at nu=' + str(nu))
        # tmp.drop('u',axis=1,inplace=True)
        energies_df = pd.concat([energies_df, tmp], ignore_index=True)
    if not os.path.isdir(results_dir_path_plot_vs_nu):
        os.makedirs(results_dir_path_plot_vs_nu)
    with open(results_dir_path_plot_vs_nu + 'folder_names.txt', 'w') as folder_names_file:
        for folder in folder_names:
            folder_names_file.write(folder + '\n')
    energies_df.to_csv(results_dir_path_plot_vs_nu + 'energies_vs_nu.csv', index=False)
    # print(energies_df.dtypes)

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
    plt.legend(loc='upper right', bbox_to_anchor=(1.45, 1.1))
    # plt.legend(bbox_to_anchor=(0.8, 0.7))
    plt.rcParams["figure.figsize"] = (10, 5)
    plt.title('Energies as function of filling factor for U=' + str(u_zero) + 'meV.', fontsize=19, y=-0.24, x=0.55)
    plt.xlabel(r'$\nu$')
    plt.ylabel('Energies(meV)')
    f.savefig(results_dir_path_plot_vs_nu + 'energies_vs_nu.pdf', bbox_inches='tight')


def plot_transitions_vs_nu():
    # pass

    transitions_df = pd.DataFrame([])
    for nu in range(-6, 7):
        number_occupied_bands_local = nu + 8
        results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local)

        try:
            folder_names_file = open(results_dir_path_local + '/folder_list.txt', 'r').read()
        except FileNotFoundError:
            print('add file folders to folder_list')
        else:
            folder_name = folder_names_file.splitlines()[-1]
        tmp = pd.read_csv(results_dir_path_local + folder_name + 'transitions_' + 'nu_' + str(nu) + '.csv')  # ,header=None

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
    plt.legend(loc='upper right', bbox_to_anchor=(1.22, 0.6))
    # plt.legend(bbox_to_anchor=(0.7, 0.75))
    plt.rcParams["figure.figsize"] = (10, 5)
    plt.title('Transition energies as function of filling factor for U=' + str(u_zero) + 'meV.', fontsize=19, y=-0.24, x=0.55)
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

def plot_total_hf_energy(nu):
    number_occupied_bands_local = nu + 8
    results_dir_path_local = dir_path + '/results/results_' + current_date + '/occupation_' + str(number_occupied_bands_local)
    try:
        folder_names_file = open(results_dir_path_local + '/folder_list.txt', 'r').read()
    except FileNotFoundError:
        print('add file folders to folder_list')
    else:
        folder_names_list = folder_names_file.splitlines()
        # print(folder_names_list)
        total_energy_df = pd.DataFrame([], columns=['u'])
        for index, folder, in enumerate(folder_names_list):
            total_energy_file = results_dir_path_local + folder + 'Et_' + 'nu_' + str(nu) + '.csv'
            tmp = pd.read_csv(total_energy_file, names=['u', 'Et_' + str(index)])  # ,header=None
            tmp['Et_' + folder[1:20] + '_real'] = tmp['Et_' + str(index)].apply(lambda x: complex(x).real)
            tmp['Et_' + folder[1:20] + '_imag'] = tmp['Et_' + str(index)].apply(lambda x: complex(x).imag)
            # tmp.index=tmp['u']
            tmp.drop(columns=['Et_' + str(index)], inplace=True)
            # tmp.drop(columns=['Et_'+str(index),'u'],inplace=True)
            # tmp.insert(0, 'nu', nu)
            # tmp = tmp[round(tmp['u'], 4) == u_zero]
            # if tmp.empty:
            #     print('u_zero=' + str(u_zero) + 'meV not found for transitions at nu=' + str(nu))
            # tmp.drop('u',axis=1,inplace=True)
            # total_energy_df = pd.concat([total_energy_df, tmp])
            # total_energy_df = total_energy_df.merge(tmp,how='outer',on=['u'])
            # total_energy_df = pd.merge_ordered(total_energy_df,tmp,how='outer',on=['u'])
            try:
                total_energy_df = pd.merge_ordered(total_energy_df, tmp, how='outer', on='u')
            except ValueError:
                total_energy_df = pd.concat([total_energy_df, tmp])
        # tmp = pd.read_csv(results_dir_path_local + 'energies_' + 'nu_' + str(nu) + '.csv')  # ,header=None
        print(total_energy_df)
        total_energy_df.to_csv(results_dir_path_local + '/total_hf_energy.csv', index=False)

        f = plt.figure()
        # font = {'size': 15}
        # plt.rc('font', **font)
        ax = plt.gca()
        plt.rcParams['figure.dpi'] = 150
        # transitions.plot(x='nu', linestyle=':', linewidth=1, markersize=12, ax=ax)
        # for i, line in enumerate(ax.get_lines()):
        #     line.set_marker(markers_list[i])
        #     line.set_label(label_list[i])
        #     line.set_color(color_list[i])
        for total_energy in total_energy_df.filter(like='real', axis=1).columns:
            total_energy
            # style_transition = transitions_style_dic.get(transition)
            total_energy_df[total_energy_df[total_energy].notna()].plot(x='u', y=total_energy, label=total_energy, style='-', ax=ax)  # , marker='o')

        plt.legend(bbox_to_anchor=(0.62, 1))
        plt.rcParams["figure.figsize"] = (10, 5)
        plt.title('total energies ', fontsize=19, y=-0.24, x=0.55)
        plt.xlabel(r'$u$(meV)')
        plt.ylabel('Total Hartree-Fock energy (meV)')
        f.savefig(results_dir_path_local + '/total_hf_energy.pdf', bbox_inches='tight')


if __name__ == "__main__":
    # energies_df = pd.read_csv(results_dir_path + 'energies_' + file_name_csv)
    # # transitions_df = pd.read_csv(results_dir_path + 'transitions_' + file_name_csv)
    # # print(energies_df)
    # plot_energies(energies_df)
    # # # plot_transitions(transitions_df)
    #
    # # plot_energies_with_asymmetry(nu=0)
    # plot_energies_vs_nu()
    # plot_transitions_vs_nu()
    plot_total_hf_energy(nu)
