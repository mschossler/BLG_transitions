import sys

import pandas as pd
from matplotlib import pyplot as plt

# from config import aux_dir_path, namecsv, title, t0, bands
# from input.parameters import nu

if __name__ == "__main__":
    # setting path
    sys.path.append('../')

from config import aux_dir_path, namecsv, bands
from input.parameters import nu

# style_dict = {'LL0_Kp_Sdown': ('lightblue', '-', 'v', r'$\ \ \,0\mathrm{K}^{+}\downarrow$'),
#               'LL1_Kp_Sdown': ('salmon', '-', 'v', r'$\ \ \,1\mathrm{K}^{+}\downarrow$'),
#               'LLm2_Kp_Sdown': ('gray', '-', 'v', r'$-2\mathrm{K}^{+}\downarrow$'),
#               'LL2_Kp_Sdown': ('gray', '-', 'v', r'$\ \ \,2\mathrm{K}^{+}\downarrow$'),
#               'LL0_Km_Sdown': ('lightblue', '--', 'v', r'$\ \ \,0\mathrm{K}^{-}\downarrow$'),
#               'LL1_Km_Sdown': ('salmon', '--', 'v', r'$\ \ \,1\mathrm{K}^{-}\downarrow$'),
#               'LLm2_Km_Sdown': ('gray', '--', 'v', r'$-2\mathrm{K}^{-}\downarrow$'),
#               'LL2_Km_Sdown': ('gray', '--', 'v', r'$\ \ \,2\mathrm{K}^{-}\downarrow$'),
#               'LL0_Kp_Sup': ('blue', '-', '^', r'$\ \ \,0\mathrm{K}^{+}\uparrow$'),
#               'LL1_Kp_Sup': ('red', '-', '^', r'$\ \ \,1\mathrm{K}^{+}\uparrow$'),
#               'LLm2_Kp_Sup': ('black', '-', '^', r'$-2\mathrm{K}^{+}\uparrow$'),
#               'LL2_Kp_Sup': ('black', '-', '^', r'$\ \ \,2\mathrm{K}^{+}\uparrow$'),
#               'LL0_Km_Sup': ('blue', '--', '^', r'$\ \ \,0\mathrm{K}^{-}\uparrow$'),
#               'LL1_Km_Sup': ('red', '--', '^', r'$\ \ \,1\mathrm{K}^{-}\uparrow$'),
#               'LLm2_Km_Sup': ('black', '--', '^', r'$-2\mathrm{K}^{-}\uparrow$'),
#               'LL2_Km_Sup': ('black', '--', '^', r'$\ \ \,2\mathrm{K}^{-}\uparrow$')}

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
        # energies.plot(x='u', y=band, color=styleX[0], style=styleX[1], markersize=3, linewidth=0.7, label=styleX[3], ax=ax)  # , marker='o')
        energies.plot(x='u', y=band, color=styleX['color'], style=styleX['line_shape'], markersize=3, linewidth=0.7, label=styleX['label'], ax=ax)  # , marker='o')
    plt.title('Energy bands  nu=' + str(nu) + ' as function of U with self-energy alpha 1')
    plt.xlabel('U(meV)')
    plt.ylabel('Energy bands(meV)')

    # plt.legend(bbox_to_anchor=(1, 0.55))
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    f.savefig(aux_dir_path + "LL(U)_HF_interactions_w_SE_warping_alpha1_nu_" + str(nu) + ".pdf", bbox_inches='tight')


if __name__ == "__main__":
    energies_df = pd.read_csv(aux_dir_path + namecsv)
    print(energies_df)
    plot_energies(energies_df)

transitions_plot_dic = {'LL1_Kp_Sdown_to_LL2_Kp_Sdown': {'color': 'tab:blue', 'marker_shape': '2', 'line_shape': '-',
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


def plot_transitions(energies):
    f = plt.figure()
    ax = plt.gca()

    for band in bands:
        styleX = style_dict.get(band)
        # energies.plot(x='u', y=band, color=styleX[0], style=styleX[1], markersize=3, linewidth=0.7, label=styleX[3], ax=ax)  # , marker='o')
        energies.plot(x='u', y=band, color=styleX['color'], style=styleX['line_shape'], markersize=3, linewidth=0.7, label=styleX['label'], ax=ax)  # , marker='o')
    plt.title('Energy bands  nu=' + str(nu) + ' as function of U with self-energy alpha 1')
    plt.xlabel('U(meV)')
    plt.ylabel('Energy bands(meV)')

    # plt.legend(bbox_to_anchor=(1, 0.55))
    plt.legend(loc='upper right', bbox_to_anchor=(1, 1))
    plt.rcParams["figure.figsize"] = (10, 5)
    # plt.show()
    f.savefig(aux_dir_path + "LL(U)_HF_interactions_w_SE_warping_alpha1_nu_" + str(nu) + ".pdf", bbox_inches='tight')
