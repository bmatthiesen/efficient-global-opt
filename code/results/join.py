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

dset_names = ('afsndEE_wp2_gurobi', 'afsndEE_wp_gurobi', 'afsndEE_wp3_gurobi')
joined_name = 'afsndEE_gurobi'

f = h5py.File('results.h5','r+')

P = np.concatenate([f[d]['input']['P'] for d in dset_names])

if not [x<y for x,y in zip(P,P[1:])]: raise RuntimeError('Powers not increasing')
if not len(set([y-x for x,y in zip(P,P[1:])]))==1: raise RuntimeError('Powers not evenly spaced')

dset = f.create_group(joined_name)
f[dset_names[0]].copy('input', dset)
del dset['input']['P']
dset['input']['P'] = P
del P

dset['raw_results'] = np.concatenate([f[d]['raw_results'] for d in dset_names], axis=1)
dset['joint_results'] = np.concatenate([f[d]['joint_results'] for d in dset_names], axis=1)

if 'individual_wp_results' not in f:
    f.create_group('individual_wp_results')
old = f['individual_wp_results']

for d in dset_names:
    old[d] = f[d]
    del f[d]

f.close()
