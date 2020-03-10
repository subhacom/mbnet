# fig_supp_steady_pn_lognorm_syn.py --- 
# Author: Subhasis  Ray
# Created: Fri Mar 15 16:04:34 2019 (-0400)
# Last-Updated: Fri Mar 15 16:13:02 2019 (-0400)
#           By: Subhasis  Ray
# Version: $Id$

# Code:

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
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import pn_kc_ggn_plot_mpl as myplot
import neurograph as ng
import timeit
from sklearn import cluster, preprocessing

plt.rc('font', size=11)

_ur = pint.UnitRegistry()
Q_ = _ur.Quantity
datadir = 'D:\\biowulf_stage\\olfactory_network\\'
datadir = '/data/rays3/ggn/olfactory_network/'


jid = '22251511'

fname = nda.find_h5_file(jid, datadir)

fig, ax = plt.subplots()
with h5.File(fname, 'r') as fd:
    print(jid, len(nda.get_spiking_kcs(fd)))
    config = nda.load_config(fd)
    print(yaml.dump(config, default_style=''))
    # pn_st = []
    # pn_id = []
    # for pn in fd[nda.pn_st_path].values():
    #     pn_st.append(pn.value)
    #     pn_id.append([int(pn.name.rpartition('_')[-1])] * len(pn))    
    # ax0.plot(np.concatenate(pn_st[::10]), np.concatenate(pn_id[::10]), marker='s', ms=1, color='#fdb863', ls='none')
    # kc_x, kc_y = nda.get_event_times(fd[nda.kc_st_path])
    # ax1.plot(np.concatenate(kc_x[::10]), np.concatenate(kc_y[::10]), color='#e66101', marker='s', ms=1, ls='none')
    myplot.plot_ggn_vm(ax, fd, fd['/data/uniform/ggn_basal/GGN_basal_Vm'], 'dend_b', 1, color='#009292', alpha=1.0)
    ax.set_ylim(-53, -45)
    ax.set_xlim(200, 2000)
    xticks = np.array([200.0, 1000.0, 2000.0])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks/1000.0)
    ax.set_yticks([-50.0, -45.0])
    ax.hlines(y=-53.0, xmin=Q_(config['stimulus']['onset']).to('ms').m, xmax=Q_(config['stimulus']['onset']).to('ms').m + Q_(config['stimulus']['duration']).to('ms').m, lw=3, color='lightgray')
    [sp.set_visible(False) for sp in ax.spines.values()]
    ax.tick_params(right=False, top=False)
    ax.xaxis.set_visible(False)
    # ax.yaxis.set_visible(False)

fig.subplots_adjust(hspace=0.01)    
fig.set_frameon(False)
fig.set_size_inches(5.5/2.54, 3.0/2.54)
fig.savefig('steady_pn_lognorm_syn_sustained_ggn_dep.svg')
plt.show()


# 
# fig_supp_steady_pn_lognorm_syn.py ends here
