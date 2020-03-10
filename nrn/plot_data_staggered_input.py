# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 15:23:31 2016

@author: Subhasis
"""

import neurograph as ng
from nrnutils import plot_syn_Vm

outfilename = 'data/Vm_inhpoisson_stim_series_20161122_152140.h5'
with h5.File(outfilename) as fd:
    
fig = plot_syn_Vm(outfilename, ng.colormaps['7q'], alpha=0.5)
for ax in fig.axes:
    ax.set_xlim(onset - 10, t_stop)
    ax.set_ylim(-55, -10)
plt.show()
fig.savefig('{}.png'.format(outfilename))
