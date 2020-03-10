# plot_ggn_vm_ts.py --- 
# Author: Subhasis Ray
# Created: Fri Sep  7 16:57:34 2018 (-0400)
# Last-Updated: Mon Jan  7 14:54:43 2019 (-0500)
#           By: Subhasis Ray
# Version: $Id$

# Code:
from __future__ import print_function
import sys
sys.path.append('D:/subhasis_ggn/model/nrn')
sys.path.append('D:/subhasis_ggn/model/morphutils')
sys.path.append('D:/subhasis_ggn/model/mb')
import argparse
import random
from itertools import chain
from datetime import datetime
import numpy as np
import h5py as h5
import matplotlib.lines as mlines
from matplotlib import pyplot as plt
from matplotlib import gridspec
from collections import defaultdict
import neurograph as ng
from morph3d_vtk import neuron3d
import nrnutils as nu
from nhpp import nhpp_thinning
from neuron import h

plt.rc('font', size=8)

if __name__ == '__main__':
    figsize = (6/2.54, 5.5/2.54)
    outfilename = sys.argv[1]
    colormap = sys.argv[2]
    with h5.File(outfilename, 'r') as fd:        
        fig = plt.figure(num=1, figsize=figsize, dpi=300)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3,1])
        ax_vm = fig.add_subplot(gs[0])   # Plot the Vm
        syninfo = [s for s in fd['syninfo']]
        sid_data_dict = defaultdict(list)
        for nodename in fd:
            if nodename.startswith('v_'):
                sid = fd[nodename].attrs['sid']
                sid_data_dict[sid].append(nodename)
        # Here I plot a random sampling of 5 or fewer datasets from each structure id
        cmap = ng.colormaps[colormap]
        syn_plotted = []
        for sid, nodelist in sid_data_dict.items():
            print('sid', sid)
            nodes = random.sample(nodelist, min(5, len(nodelist)))
            color = np.array(cmap[sid % len(cmap)]) / 255.0
            for nodename in nodes:
                sec = fd[nodename].attrs['section']
                ls = '-'
                # if sec in syninfo:
                #     ls = '--'
                #     syn_plotted.append(sec)
                # else:
                #     ls = '-'                
                ax_vm.plot(fd['time'], fd[nodename], color=color, ls=ls, lw=2.0, label=sec)
        handles = []
        sid_label_map = {1: 'soma', 3: 'basal', 5: 'LCA', 6: 'MCA', 7: 'LH', 8: r'$\alpha$ lobe'}
        for sid in sid_data_dict.keys():
            if sid == 1:
                continue # skip soma
            handles.append(mlines.Line2D([], [], color=np.array(cmap[sid % len(cmap)]) / 255.0,
                                         label=sid_label_map[sid], lw=2))
        ax_vm.legend(handles=handles, loc='lower left', bbox_to_anchor=(1.05, 0), borderaxespad=0.)
        # ax_vm.legend(handles=handles, loc='best')
        # ax_vm.set_xlim(onset-10, t_stop)
        # Plot the input spike raster, put the ones delivered to the plotted sections on top and in red
        ax_stim = fig.add_subplot(gs[1], sharex=ax_vm)
        syn_color = 'yellow'
        other_color = 'blue'
        syn_stim = []
        other_stim = []
        for ii, stimnode in enumerate(fd['stimulus']):
            ds = fd['stimulus'][stimnode]
            # if ds.attrs['section'] in syn_plotted:
            #     syn_stim.append((ii, ds))
            # else:
            #     other_stim.append((ii, ds))
            ax_stim.plot(ds, np.ones_like(ds)*ii, color=other_color, ls='', marker='|', mew=2, alpha=0.3)
        # for ii, ds in other_stim:
        #     ax_stim.plot(ds, np.ones_like(ds)*ii, color=other_color, ls='', marker='|', mew=2, alpha=0.3)
        # for ii, ds in syn_stim:
        #     ax_stim.plot(ds, np.ones_like(ds)*ii, color=syn_color, ls='', marker='|', mew=2, alpha=1.0)
    
        handles = [mlines.Line2D([], [], color=syn_color,
                                 label='plotted', lw=2),
                   mlines.Line2D([], [], color=other_color,
                                 label='not plotted', lw=2)]
        # ax_stim.legend(handles=handles)
        # ax_vm.set_ylim(-55, -35.0)
        ax_vm.set_yticks([-50, -45, -40])
        # ax_stim.set_ylim(0, syncount)
        # plt.tight_layout()
        ax_vm.xaxis.set_visible(False)
        for ax in [ax_vm, ax_stim]:
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
        # fig.set_size_inches(*figsize)
        # fig.tight_layout()
        fig.subplots_adjust(left=0.15, right=0.6, top=0.99, bottom=0.15, wspace=0.1, hspace=0.1)
        figfile = '{}.pdf'.format(outfilename)
        fig.savefig(figfile, transparent=True)
        print('Figure saved in', figfile)
        # ax_stim.set_xlim((onset-10, h.tstop))
        # figfile = '{}.zoomed.png'.format(outfilename)
        # fig.savefig(figfile, transparent=True)                
        # print('Figure saved in', figfile)
        plt.show()



# 
# plot_ggn_vm_ts.py ends here
