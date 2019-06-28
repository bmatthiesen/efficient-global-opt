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

import pathlib
import h5py
import numpy as np
import itertools as it

ofn = 'benchmark.h5'

pn = ['2_benchmark', '3_benchmark', '7_benchmark']
pbase = pathlib.Path('/scratch/p_mwrc/')
p = [pbase / pp / 'results' for pp in pn]
glob = 'res_sndBench*h5'

numChan = 1000

# select files
bench2 = p[0].glob('res_sndBench*_2_*h5')
bench3 = p[1].glob('res_sndBench*_3_*h5')
bench7 = p[2].glob('res_sndBench*h5') # add 2 and 3 from other results!

# raw results
of = h5py.File(pbase / ofn, 'w')
g = of.create_group('raw results')
g.create_group('NC_2')
g.create_group('NC_3')
g.create_group('C_7')

for fn in it.chain(bench2, bench3):
    #print(fn)
    with h5py.File(fn, 'r') as f:
        ue = max(f['results']['xopt_C'].shape)
        nc = max(f['results']['xopt_NC'].shape)

        if nc not in [2, 3]:
            print('error: {}'.format(fn))
            continue

        g = of['raw results/NC_{}'.format(nc)]

        uek = 'C_{}'.format(ue)

        if uek not in g:
            g.create_dataset(uek, dtype = f['results'].dtype, shape = (numChan, ))

        g[uek][f['results']['WP index'][0]] = f['results'][0]

for fn in bench7:
    #print(fn)
    with h5py.File(fn, 'r') as f:
        ue = max(f['results']['xopt_C'].shape)
        nc = max(f['results']['xopt_NC'].shape)

        if ue != 7:
            print('error: {}'.format(fn))
            continue

        g = of['raw results/C_{}'.format(ue)]

        uek = 'NC_{}'.format(nc)

        if uek not in g:
            g.create_dataset(uek, dtype = f['results'].dtype, shape = (numChan, ))

        g[uek][f['results']['WP index'][0]] = f['results'][0]

# runtime_numUE
idx2 = [int(s.split('_')[1]) for s in of['raw results']['NC_2'].keys()]
idx3 = [int(s.split('_')[1]) for s in of['raw results']['NC_3'].keys()]

assert(len(range(min(idx2), max(idx2)+1)) == len(idx2))
assert(len(range(min(idx3), max(idx3)+1)) == len(idx3))

numNC = (2, 3)
numUE = range(min(min(idx2), min(idx3)), max(max(idx2), max(idx3))+1)

g = of.create_group('runtime_numUE')
g.create_dataset('dim2_numUE', data = numUE)
g.create_dataset('dim1_numNC', data = numNC)
g.create_dataset('data', dtype = 'f8', shape = (2, len(numUE), numChan), fillvalue = np.nan)

for ncIdx in numNC:
    nc = 'NC_{}'.format(ncIdx)

    for ue in of['raw results'][nc].keys():
        ueIdx = int(ue.split('_')[1])
        
        g['data'][numNC.index(ncIdx), numUE.index(ueIdx)] = of['raw results'][nc][ue]['Runtime']

zeros = g['data'][...]==0
g['data'][zeros] = np.nan

# runtime_numNC
ueIdx = 7
ue = 'C_{}'.format(ueIdx)

idx7 = [int(s.split('_')[1]) for s in of['raw results'][ue].keys()] + [2,3]
numNC = range(min(idx7), max(idx7)+1)

g = of.create_group('runtime_numNC')
g.create_dataset('dim1_numNC', data = numNC)
g.create_dataset('data', dtype = 'f8', shape = (len(numNC), numChan), fillvalue = np.nan)

for ncIdx in numNC:
    nc = 'NC_{}'.format(ncIdx)

    try:
        rt = of['raw results'][ue][nc]['Runtime']
    except KeyError:
        rt = of['raw results'][nc][ue]['Runtime']

    g['data'][numNC.index(ncIdx)] = rt

zeros = g['data'][...]==0
g['data'][zeros] = np.nan

of.close()
