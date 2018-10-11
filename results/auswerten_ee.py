import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

TINIDX = 0
SNDIDX = 7

try:
    h5py.enable_ipython_completer()
except RuntimeError:
    pass

hdf = h5py.File('results.h5')

objval = pd.DataFrame(index=hdf['afsndEE_gurobi']['input']['P'][:])
runtime = pd.DataFrame(index=hdf['afsndEE_gurobi']['input']['P'][:])

def addDset(name, dset):
    if np.any(np.isnan(dset['Objective Value'])):
        print("Dataset '{}' contains {} NaNs".format(name, np.sum(np.isnan(dset['Objective Value']))))

    global objval, runtime
    objval = objval.assign(**{name: np.nanmean(np.log2(np.e)*dset['Objective Value'], axis=0)})
    runtime = runtime.assign(**{name + '_mean': np.nanmean(dset['Runtime'], axis=0), name + '_median': np.nanmedian(dset['Runtime'], axis=0)})

addDset("SND", hdf['afsndEE_1e-3']['joint_results'])
addDset("TIN", hdf['afsndEE_1e-3']['raw_results'][:,:,TINIDX])
addDset("Joint", hdf['afsndEE_1e-3']['raw_results'][:,:,SNDIDX])

addDset("SND_gurobi", hdf['afsndEE_gurobi']['joint_results'])
addDset("TIN_gurobi", hdf['afsndEE_gurobi']['raw_results'][:,:,TINIDX])
addDset("Joint_gurobi", hdf['afsndEE_gurobi']['raw_results'][:,:,SNDIDX])

addDset("SND_mosek", hdf['afsndEE_mosek']['joint_results'])
addDset("TIN_mosek", hdf['afsndEE_mosek']['raw_results'][:,:,TINIDX])
addDset("Joint_mosek", hdf['afsndEE_mosek']['raw_results'][:,:,SNDIDX])

addDset("Dinkelbach", hdf['afsnd_dinkelbach']['raw_results'])

avg_gain = pd.concat((
        pd.concat((objval['SND_gurobi']/objval['TIN_gurobi'], objval['SND_gurobi']/objval['Joint_gurobi']), axis=1, keys=['SND vs TIN [%]', 'SND vs Joint [%]'])*100-100,
        pd.concat((objval['SND_gurobi']-objval['TIN_gurobi'], objval['SND_gurobi']-objval['Joint_gurobi']), axis=1, keys=['SND vs TIN [bpcu]', 'SND vs Joint [bpcu]'])
        ), axis=1)

objval.to_csv('ee.dat', index_label='snr')
runtime.to_csv('ee_runtime.dat', index_label='snr')
