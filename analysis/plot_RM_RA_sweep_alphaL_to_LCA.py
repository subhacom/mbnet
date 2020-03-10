# -*- coding: utf-8 -*-
"""
Graphical summary of the RA-RM sweep.

@author: Subhasis Ray / ray dor subhasis at gmail dot com
"""
from __future__ import print_function
import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import neurograph as ng

plt.rc('font', size=32)
#summary_file = 'D:\subhasis_ggn\model\summary_data\single_syn_summary_20170403.csv'
#summary_file = 'D:\subhasis_ggn\model\summary_data\single_syn_summary_20170405.csv'
#summary_file = 'D:\\subhasis_ggn\\model\\summary_data\\single_syn_summary_20170406.csv'

colors = ['r', 'g', 'b', 'c', 'y', 'k']  # colors for recording region
markers = ['s', 'o', '^', 'v', 'x', '+']   # Markers for synaptic region
regions = ['5', '6', '7', '8', 'dend', 'soma']
# 5: 'LCA', 6:'MCA', 7:'LH', 8:'alphaL'

summary_stats_file = 'D:\\subhasis_ggn\\model\\summary_data\\single_syn_summary_20170407.csv_normalized_vm_deflection_stats.csv'

#############################
# 2D plot the distribution of relative deflection at LCA for input at alphaL
##############################################
from matplotlib import ticker
# plt.style.use('light_background')
df = pd.read_csv(summary_stats_file)
alphaL_to_LCA = df[(df['synregion'] == 'alphaL') & (df['recregion'] == 'LCA')]
RA_vals = sorted(list(set(alphaL_to_LCA['RA'])))
RM_vals = sorted(list(set(alphaL_to_LCA['RM'])))
print(RM_vals)
RA_grp = alphaL_to_LCA.groupby('RA')    
fig, ax = plt.subplots(nrows=1, ncols=1)
for ii, RA in enumerate(RA_vals[:6]):
    df_RA = RA_grp.get_group(RA)
    deflections = [df_RA[df_RA['RM'] == RM]['mean'].item() for RM in RM_vals]
    yerr = [df_RA[df_RA['RM'] == RM]['std'].item() for RM in RM_vals]
    ax.plot([1e-3*RM for RM in RM_vals], deflections, lw=5)
    # ax.errorbar([1e-3*RM for RM in RM_vals], deflections, yerr=yerr, lw=5)
    ax.text(RM_vals[int(len(RM_vals)*2/3)]*1e-3 , deflections[-1], '{:.2f}'.format(RA))
    # ax.text(RM_vals[int(len(RM_vals)*2/3)]*1e-3 , deflections[-1], 'RA={:.2f} ohm-cm'.format(RA))
ax.set_frame_on(False)    
#     ax.tick_params(top='off', bottom='off', left='off', right='off', labelbottom='off')
#     ax.set_title('RA={:.2f} ohm-cm'.format(RA))
#     ax.set_yticks([0, 0.5, 1.0])        
#     ax.spines['top'].set_visible(False)
#     ax.spines['right'].set_visible(False)
#     #ax.spines['bottom'].set_visible(False)
#     ax.spines['left'].set_visible(False)
# ax.tick_params(labelleft='on', labelbottom='on')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
ax.set_xlabel(r'Membrane resistivity ($\mathrm{k\Omega\/cm^2}$)')
ax.set_ylabel('Mean Vm deflection (relative)')
#ax.set_ylabel('Normalized Vm defelction')
# fig.tight_layout()
fig.set_size_inches((16, 8))
fig.savefig('20170407_summary_stats.svg')
plt.show()
    
