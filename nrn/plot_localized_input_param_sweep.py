# -*- coding: utf-8 -*-
"""
Plot parameter sweep data for synchronized synaptic input to GGN

Created on Fri May 12 12:40:50 2017

@author: Subhasis
"""
from __future__ import print_function
import numpy as np
import h5py as h5
from matplotlib import pyplot as plt

        
with h5.File('data/B_Vm_multi_syn_20170510_171812.h5', 'r') as fd:
    for ii in range(10):
        fig = plt.figure()        
        for jj in range(10):            
            simno = ii * 10 + jj
            simgrp = fd['simulation_{}'.format(simno)]
            print('simgroup', simgrp)
            syninfo = simgrp['syninfo']
            RA = simgrp.attrs['RA']
            RM = 1.0/simgrp.attrs['g_pas']
            print(RA, RM)
            ax = fig.add_subplot(5,2, jj+1)
            ax.set_title('RA {:0.2f} ohm-cm, RM {:.2f} kohm-cm2'.format(RA, RM/1000))
            t = simgrp['time']
            for node in simgrp:
                if not node.startswith('v_'):
                    continue
                v = simgrp[node]
                sec = v.attrs['section']
                if sec in syninfo:
                    # print('Synapse', sec)
                    ax.plot(t, v, 'g--')
                elif 'dend_5' in str(sec):
                    ax.plot(t, v, 'y-')
        fig.set_size_inches((8, 6))
        fig.savefig('Vm_multi_syn_20170510_171812_{}.png'.format(ii), transparent=True)
#plt.tight_layout()        
plt.show()        