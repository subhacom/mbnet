# plot_shifting_pn_lognorm_syn.py --- 
# Author: Subhasis Ray
# Created: Thu Apr 11 16:23:57 2019 (-0400)
# Last-Updated: Thu Apr 11 16:43:50 2019 (-0400)
#           By: Subhasis Ray
# Version: $Id$

# Code:


"""This script plots shifting PN and GGN response for lognormally
distributed synaptic strengths"""
#* imports
from __future__ import print_function
import sys
if sys.platform == 'win32':
    sys.path.append('D:/subhasis_ggn/model/analysis')
else:
    sys.path += ['/home/rays3/projects/ggn/analysis',
                 '/home/rays3/projects/ggn/morphutils',
                 '/home/rays3/projects/ggn/nrn']
import os
import h5py as h5
import numpy as np
import random
from matplotlib import pyplot as plt
import yaml
import pint
from collections import defaultdict
import pandas as pd
import network_data_analysis as nda
import matplotlib.gridspec as gridspec
import pn_kc_ggn_plot_mpl as myplot

plt.rc('font', size=8)

_ur = pint.UnitRegistry()
Q_ = _ur.Quantity
# datadir = 'Z:/Subhasis/ggn_model_data/olfactory_network'
# datadir = '/data/rays3/ggn/olfactory_network/'
datadir = 'D:/biowulf_stage/olfactory_network'

#** jid_sc for only shifting PN, constant GGN->KC inhibition - but this did not
# have PN-KC clustered conn

jid_sc = '22087969'  
fname_sc = nda.find_h5_file(jid_sc, datadir)
fd_sc = h5.File(fname_sc, 'r')
print('shifting PN, jid: {}, spiking KCs {}'.format(jid_sc, len(nda.get_spiking_kcs(fd_sc))))
print(yaml.dump(nda.load_config(fd_sc), default_style=''))
print('-' * 20)

stiminfo = nda.get_stimtime(fd_sc)
pn_st = []
pn_id = []
for pn in fd_sc[nda.pn_st_path].values():
    pn_st.append(pn[:])
    pn_id.append([int(pn.name.rpartition('_')[-1])] * len(pn))

kc_st = []
kc_id = []
for kc, st in fd_sc[nda.kc_st_path].items():
    kc_st.append(st[:])
    kc_id.append([int(kc)] * len(st))

fig, ax = plt.subplots(nrows=3, sharex='all')

ax[0].plot(np.concatenate(pn_st[::10]), np.concatenate(pn_id[::10]),
         color='#fdb863', ls='', marker=',')
for a in ax[:2]:
    a.xaxis.set_visible(False)
    a.yaxis.set_visible(False)

ax[1].plot(np.concatenate(kc_st[::10]), np.concatenate(kc_id[::10]),
         color='#e66101', ls='', marker=',')

#* plot ggn vm from both simulations
myplot.plot_ggn_vm(ax[2], fd_sc, fd_sc['/data/uniform/ggn_basal/GGN_basal_Vm'],
                   'dend_b', 1, color='k', alpha=1.0)
ax[2].set_ylim((-53, -40))
ax[2].set_yticks([-50, -45])
ax[2].set_xlim(200, 2500)
ax[2].set_xticks([200, 1000, 2000])
ax[2].set_xlabel('Time (ms)')
ax[2].hlines(y=-53,
           xmin=stiminfo['onset'],
           xmax=stiminfo['onset'] + stiminfo['duration'],
           color='gray', lw=10)
ax[2].hlines(y=-53,
           xmin=stiminfo['onset'] + stiminfo['duration'],
           xmax=stiminfo['onset'] + stiminfo['duration'] + stiminfo['offdur'],
           color='lightgray', lw=10)
for a in ax.flat:
    a.tick_params(top=False, right=False)
    a.xaxis.set_visible(False)
    [sp.set_visible(False) for sp in a.spines.values()]
fig.tight_layout()
fig.set_size_inches(6/2.54, 8.0/2.54)
fig.savefig('shifting_pn_lognorm_syn.svg', transparent=True)
plt.show()


# 
# plot_shifting_pn_lognorm_syn.py ends here
