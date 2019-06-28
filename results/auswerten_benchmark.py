# Copyright (C) 2018-2019 Bho Matthiesen
# 
# This program is used in the article:
# 
# Bho Matthiesen and Eduard Jorswieck, "Efficient Global Optimal Resource
# Allocation in Non-Orthogonal Interference Networks," submitted to IEEE
# Transactions on Signal Processing.
# 
# 
# License:
# This program is licensed under the GPLv2 license. If you in any way use this
# code for research that results in publications, please cite our original
# article listed above.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data_cutoff = .028 # 2.5%

f = h5py.File('benchmark.h5','r')

def getData(g):
    useless = np.sum(np.isnan(g), axis=-1) >= int(g.shape[-1] * data_cutoff)
    m = np.nanmean(g, axis=-1)
    m[useless] = np.nan

    return m

def fit(d, start = None, end = None):
    if start is None:
        start = min(np.where(~np.isnan(d))[0])

    if end is None:
        end = start + 1

    y = d[start:end+1]
    x = d.index[start:end+1]

    a,b=np.polyfit(x, np.log(y), 1)
    b = np.exp(b)

    return lambda x: b * np.exp(a*x)

g = f['runtime_numUE']

gd = g['data'][...]
#gdnan = np.isnan(gd[1,-3:])
#assert(np.sum(gdnan) == 17)
#gd[1,-3:][gdnan] = 7*24*3600

mUE = getData(gd).T
rUE = pd.DataFrame(mUE, index=g['dim2_numUE'][:], columns=['NC{}'.format(d) for d in g['dim1_numNC']])

idx = np.asarray(rUE.index)
rUE['NC2_fit']=pd.Series(fit(rUE['NC2'], start=9, end=12)(idx), index=idx)
rUE['NC3_fit']=pd.Series(fit(rUE['NC3'], start=15, end=18)(idx), index=idx)

g = f['runtime_numNC']
mNC = getData(g['data'])
rNC = pd.DataFrame(mNC, index = g['dim1_numNC'][:], columns=['UE7'])

rUE.to_csv('benchUE.dat', na_rep='NaN')
rNC.to_csv('benchNC.dat', na_rep='NaN')

#rUE.plot(logy=True)
#plt.show()
