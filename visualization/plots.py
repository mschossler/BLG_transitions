import pandas as pd
import time
from matplotlib import pyplot as plt
from config import aux_dir_path, namecsv, title, t0


ax = plt.gca()

# dfv0 = pd.read_csv(  aux_dir_path + 'nu4_v8_full_range_rho_v3_25itmax.csv', names=['u', 'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7', 'band8'])

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

df = pd.read_csv(  aux_dir_path + namecsv, names=['u', 'band1', 'band2', 'band3', 'band4', 'band5', 'band6', 'band7', 'band8'])

df.plot(kind='line', x='u', y='band1', ax=ax)
df.plot(kind='line', x='u', y='band2', color='red', ax=ax)
df.plot(kind='line', x='u', y='band3', color='blue', ax=ax)
df.plot(kind='line', x='u', y='band4', color='orange', ax=ax)
df.plot(kind='line', x='u', y='band5', color='green', ax=ax)
df.plot(kind='line', x='u', y='band6', color='purple', ax=ax)
df.plot(kind='line', x='u', y='band7', color='brown', ax=ax)
df.plot(kind='line', x='u', y='band8', color='coral', ax=ax)
# df.plot(kind='line',x='u',y='band2', color='gray', ax=ax)
# df.plot(kind='line',x='u',y='band2', color='olive', ax=ax)
plt.legend('', frameon=False)
# plt.show()

# ts = pd.Series(np.random.randn(1000), index=pd.date_range('1/1/2000', periods=1000))
#
# ts = ts.cumsum()
#
# ts.plot()
plt.xlim([-0.04, 0.05])
plt.ylim([-125, 75])

plt.grid(True)

plt.savefig('plot_' + title + '.png')
print(" \n done in ", time.time() - t0)