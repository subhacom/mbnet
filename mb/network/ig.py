# ig.py --- 
# Author: Subhasis Ray
# Created: Mon Jun 25 18:54:16 2018 (-0400)
# Last-Updated: Mon Nov  5 15:31:13 2018 (-0500)
#           By: Subhasis  Ray
# Version: $Id$

# Code:
"""Model of spikes as NHPP and rates taken from an npy file"""
import numpy as np
from scipy import interpolate
from config import Q_, h, logger
import nhpp

class IGSpikes(object):
    def __init__(self, filename, simtime, delta, stabilization_time=0):
        """Create spike times for IG using firing rate over time from
        filename.

        This disregards stimulus onset and duration, and only uses the
        information from the file.

        """
        spikerate = np.load(filename)
        interp = interpolate.interp1d(spikerate['t'], spikerate['Vm'])
        rate_fn = lambda t: interp(t)
        spike_times = nhpp.nhpp_thinning(rate_fn, simtime, delta) * 1e3  # second to millisecond
        spike_times = spike_times[spike_times > stabilization_time].copy()
        self.stimvec = h.Vector(spike_times)
        self.vecstim = h.VecStim()
        self.vecstim.play(self.stimvec)

    def connect(self, sec, pos, synparams):
        self.synlist = []
        self.netcons = []
        count = synparams.get('count', 1)
        for ii in range(count):
            synapse = h.Exp2Syn(sec(pos))
            synapse.e = Q_(synparams.get('e', '-80mV')).to('mV').m
            synapse.tau1 = Q_(synparams.get('tau1', '13.33ms')).to('ms').m
            synapse.tau2 = Q_(synparams.get('tau2', '13.33ms')).to('ms').m
            thresh = Q_(synparams.get('threshold', '-20mV')).to('mV').m
            delay = Q_(synparams.get('delay', '0ms')).to('ms').m
            gmax = Q_(synparams.get('gmax', '1nS')).to('uS').m        
            netcon = h.NetCon(self.vecstim, synapse,  thresh, delay,
                                   gmax)
            self.synlist.append(synapse)
            self.netcons.append(netcon)
            logger.info('Connected IG->{} at {}'.format(sec.name(), pos))
            logger.info('threshold: {}, delay: {}, gmax: {}, e: {}, tau1: {}, tau2: {}'.format(
                netcon.threshold,
                netcon.delay, netcon.weight[0], synapse.e, synapse.tau1,
                synapse.tau2))


class IGIzhi(object):
    def __init__(self, name='IG', inject=70.0, gbarGABA=1.0):
        """IG model based on Izhikevich type neuron."""
        self.sec = h.Section(name)
        self.izhi = h.IzhiGS(self.sec(0.5))
        self.izhi.celltype = 1  # Regular Spiking
        self.izhi.Iin = inject
        self.izhi.gbarGraded = gbarGABA
        self.netcons = []        
        
    def connect_ampa_vecstim(self, src, threshold=-20, delay=5.0, weight=1.0):
        """Connect a VecStim source to AMPA synapse on this neuron"""
        nc = h.NetCon(src, self.izhi, threshold, delay, weight, sec=self.sec)
        self.netcons.append(nc)        

    def connect_ampa_sec(self, sec, pos=0.5, threshold=-20, delay=5.0, weight=1.0):
        """Connect a Section source to AMPA synapse on this neuron"""
        # print('Creating connection from {} to {}'.format(sec.name(), self.sec.name()))
        sec.push()
        nc = h.NetCon(sec(pos)._ref_v, self.izhi, threshold, delay, weight, sec=sec)
        h.pop_section()
        self.netcons.append(nc)

        
        
        
# 
# ig.py ends here
