import pandas as pd
import time
from matplotlib import pyplot as plt
import sys
# from config import aux_dir_path, namecsv, title, t0, bands
# from input.parameters import nu

if __name__ == "__main__":
    # setting path
    sys.path.append('../')
    from config import aux_dir_path, namecsv, bands
    from input.parameters import  nu

else:
    from config import aux_dir_path, namecsv, title, t0, bands
    from input.parameters import nu

f = plt.figure()
ax = plt.gca()

plt.rcParams['figure.dpi'] = 150


style_dict={'LL0_Kp_Sdown':('lightblue','-','v',r'$\ \ \,0\mathrm{K}^{+}\downarrow$'),
            'LL1_Kp_Sdown':('salmon','-','v',r'$\ \ \,1\mathrm{K}^{+}\downarrow$'),
            'LLm2_Kp_Sdown':('gray','-','v',r'$-2\mathrm{K}^{+}\downarrow$'),
            'LL2_Kp_Sdown':('gray','-','v',r'$\ \ \,2\mathrm{K}^{+}\downarrow$'),
            'LL0_Km_Sdown':('lightblue','--','v',r'$\ \ \,0\mathrm{K}^{-}\downarrow$'),
            'LL1_Km_Sdown':('salmon','--','v',r'$\ \ \,1\mathrm{K}^{-}\downarrow$'),
            'LLm2_Km_Sdown':('gray','--','v',r'$-2\mathrm{K}^{-}\downarrow$'),
            'LL2_Km_Sdown':('gray','--','v',r'$\ \ \,2\mathrm{K}^{-}\downarrow$'),
            'LL0_Kp_Sup':('blue','-','^',r'$\ \ \,0\mathrm{K}^{+}\uparrow$'),
            'LL1_Kp_Sup':('red','-','^',r'$\ \ \,1\mathrm{K}^{+}\uparrow$'),
            'LLm2_Kp_Sup':('black','-','^',r'$-2\mathrm{K}^{+}\uparrow$'),
            'LL2_Kp_Sup':('black','-','^',r'$\ \ \,2\mathrm{K}^{+}\uparrow$'),
            'LL0_Km_Sup':('blue','--','^',r'$\ \ \,0\mathrm{K}^{-}\uparrow$'),
            'LL1_Km_Sup':('red','--','^',r'$\ \ \,1\mathrm{K}^{-}\uparrow$'),
            'LLm2_Km_Sup':('black','--','^',r'$-2\mathrm{K}^{-}\uparrow$'),
            'LL2_Km_Sup':('black','--','^',r'$\ \ \,2\mathrm{K}^{-}\uparrow$')}



def plot_energies(energies):
    for band in bands:
        styleX = style_dict.get(band)
        energies.plot(x='u', y=band, color=styleX[0], style=styleX[1], markersize=3, linewidth=0.7, label=styleX[3], ax=ax)#, marker='o')
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
# dfv0 = pd.read_csv(  aux_dir_path + 'nu4_v8_fuLL_range_rho_v3_25itmax.csv', names=['u', 'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7', 'band8'])

#
# dfv0.plot(kind='line',x='u',y='band1', color='black', ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band2', color='black', ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band3', color='black',ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band4', color='black', ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band5', color='black', ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band6', color='black', ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band7', color='black', ax=ax, linestyle='--',legend=None)
# dfv0.plot(kind='line',x='u',y='band8', color='black', ax=ax, linestyle='--',legend=None)

# plt.legend('')
# plt.savefig('plot2.pdf')
#
# df = pd.read_csv(  aux_dir_path + namecsv, names=['u', 'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7', 'band8'])
#
# df.plot(kind='line', x='u', y='band1', ax=ax)
# df.plot(kind='line', x='u', y='band2', color='red', ax=ax)
# df.plot(kind='line', x='u', y='band3', color='blue', ax=ax)
# df.plot(kind='line', x='u', y='band4', color='orange', ax=ax)
# df.plot(kind='line', x='u', y='band5', color='green', ax=ax)
# df.plot(kind='line', x='u', y='band6', color='purple', ax=ax)
# df.plot(kind='line', x='u', y='band7', color='brown', ax=ax)
# df.plot(kind='line', x='u', y='band8', color='coral', ax=ax)
# # df.plot(kind='line',x='u',y='band2', color='gray', ax=ax)
# # df.plot(kind='line',x='u',y='band2', color='olive', ax=ax)
# plt.legend('', frameon=False)
# # plt.show()
#
# # ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
# #
# # ts = ts.cumsum()
# #
# # ts.plot()
# plt.xlim([-0.04, 0.05])
# plt.ylim([-125, 75])
#
# plt.grid(True)
#
# plt.savefig('plot_' + title + '.png')
# print(" \n done in ", time.time() - t0)