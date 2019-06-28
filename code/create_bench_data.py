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
import sys
import platform

numChan = 1000
maxUE = 50

def crandn(*args, **kwargs):
    return 1/np.sqrt(2) * (np.random.randn(*args, **kwargs) + 1j * np.random.randn(*args, **kwargs))

f = h5py.File('../data/wp-bench.h5','w')

# store some info
g = f.create_group('channel_generation')

# save source code
g.create_dataset('source code', (1,), dtype = h5py.special_dtype(vlen=str))
with open(__file__, "r", encoding='utf') as src:
    g['source code'][:] = "".join(src.readlines())

# save relevant versions
g.create_dataset('python version', (1,), dtype = h5py.special_dtype(vlen=str))
g['python version'][:] = sys.version

g.create_dataset('numpy version', (1,), dtype = h5py.special_dtype(vlen=str))
g['numpy version'][:] = np.version.version

g.create_dataset('platform', (1,), dtype = h5py.special_dtype(vlen=str))
g['platform'][:] = platform.platform()

# save random state
state_dt = np.dtype([('before channel idx',np.uint64), ('state', [('keys', np.uint32, 624), ('pos', np.int), ('has_gauss', np.int), ('cached_gaussian', np.float)])])
g.create_dataset('MT19937 state', dtype = state_dt, shape = (1,), maxshape = (None, ), chunks = (1, ))

state = np.random.get_state()
assert(state[0] == 'MT19937')
g['MT19937 state'][0] = np.asarray((0,state[1:]),dtype=state_dt)

# make channel
f.create_dataset("channel", dtype = 'c16', shape = (numChan,maxUE,maxUE), fillvalue = np.nan + 1j *np.nan)

f['channel'][:100] = crandn(100,maxUE,maxUE)
f['channel'][100:] = crandn(numChan-100,maxUE,maxUE)

f.close()
