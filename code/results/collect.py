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
from pathlib import Path
import sys

def collect(wpfile, outfile, dsetname, respath, glob_pattern):
    wp = h5py.File(wpfile, 'r')
    out = h5py.File(outfile)

    if dsetname in out:
        del out[dsetname]

    grp = out.create_group(dsetname)
    wp.copy('input', grp)

    p = Path(respath)
    first = True
    process = False
    for fn in p.glob(glob_pattern):
        try:
            with h5py.File(str(fn),'r') as res:
                if first:
                    process = res['results'].shape[0] > 1

                    init = np.empty((wp['input']['h'].shape[0], wp['input']['P'].shape[0], res['results'].shape[0]), dtype=res['results'].dtype)
                    init[:] = np.nan
                    init[:]['Status'] = h5py.check_dtype(enum = res['results']['Status'].dtype)['Unsolved']
                    grp.create_dataset('raw_results',data=init)

                    if process:
                        init = np.empty((wp['input']['h'].shape[0], wp['input']['P'].shape[0]), dtype=res['results'].dtype)

                        init[:] = np.nan
                        init[:]['Status'] = h5py.check_dtype(enum = res['results']['Status'].dtype)['Unsolved']
                        grp.create_dataset('joint_results',data=init)

                    first = False

                idx = res['results'][0]['WP index']
                if not all(np.logical_or(res['results'][:]['WP index'] == idx, res['results'][:]['WP index'] == 0)):
                    raise RuntimeError("Can't handle different wp indices (" + str(fn) + ")")

                grp['raw_results'][wp['WP'][idx][0],wp['WP'][idx][1],:] = res['results'][:]

                if process:
                    idx2 = np.argmax(res['results']['Objective Value'])
                    tmp = res['results'][idx2]
                    tmp['Runtime'] = sum(res['results']['Runtime'])
                    tmp['Peak RSS'] = max(res['results']['Peak RSS'])
                    tmp['WP index'] = idx2
                    grp['joint_results'][wp['WP'][idx][0],wp['WP'][idx][1]] = tmp
        except OSError as e:
            print('{}: OSError: {}'.format(fn,e), file=sys.stderr)


if __name__=="__main__":
    import argparse

    def hdf_file(f):
        try:
            with open(f,'r+'):
                pass
        except Exception as e:
            raise argparse.ArgumentTypeError(str(e))

        if not h5py.is_hdf5(f):
            raise argparse.ArgumentTypeError("Not an HDF5 file: '{}'".format(f))

        return str(f)

    def hdf_file2(f):
        try:
            with open(f,'r+'):
                pass
        except FileNotFoundError:
            return str(f)
        except Exception as e:
            raise argparse.ArgumentTypeError(str(e))

        if not h5py.is_hdf5(f):
            raise argparse.ArgumentTypeError("Not an HDF5 file: '{}'".format(f))

        return str(f)

    def folder(f):
        p = Path(f)

        if not p.is_dir():
            raise argparse.ArgumentTypeError("Not a directory: '{}'".format(f))

        return str(f)

    
    p = argparse.ArgumentParser()
    p.add_argument("wpfile", help="the WP file", type=hdf_file)
    p.add_argument("outfile", help="the output file", type=hdf_file2)
    p.add_argument("dsetname", help="The dataset name")
    p.add_argument("respath", help="The path to the results folder", type=folder)
    p.add_argument("glob_pattern", help="The glob pattern for results files")

    args = p.parse_args()
    
    collect(**(vars(args)))
