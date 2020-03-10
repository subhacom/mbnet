# pn_kc_ggn_plot.py --- 
# 
# Filename: pn_kc_ggn_plot.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Fri Feb 16 13:08:41 2018 (-0500)
# Last-Updated: Fri Feb 16 16:11:54 2018 (-0500)
#           By: Subhasis Ray
#     Update #: 263
# 
# Code:

import sys
import os
from timeit import default_timer as timer
import numpy as np
import h5py as h5
import yaml
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

pg.setConfigOptions(antialias=True)
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

default_pen = (0, 0, 0, 100)

def get_event_times(group, nodes=None):
    spike_x = []
    spike_y = []
    if nodes is None:
        nodes = group
    for ii, node in enumerate(nodes):
        st = np.empty(2*len(group[node][:]))
        st[::2] = group[node][:]
        st[1::2] = group[node][:]
        spike_x.append(st)
        sy = np.zeros(st.shape)
        sy[::2] = ii
        sy[1::2] = ii + 0.5
        spike_y.append(sy)
    spike_x = np.concatenate(spike_x)
    spike_y = np.concatenate(spike_y)
    return spike_x, spike_y


def plot_spike_rasters(fname, vm_samples=10):
    """The file `fname` has data from pn_kc_ggn simulation. In the early
    ones I did not record the spike times for KCs.

    """
    start = timer()
    with h5.File(fname, 'r') as fd:
        config = yaml.load(fd['/model/filecontents/mb/network/config.yaml'][0])
        # PN spike raster
        pn_st = fd['/data/event/pn/pn_spiketime']
        gw = pg.GraphicsWindow(title=fd.attrs['description'])
        gw.setWindowTitle(fname)
        pn_plot = gw.addPlot(title='PN spike raster')
        pn_raster = pg.PlotCurveItem()
        spike_x, spike_y = get_event_times(pn_st)
        pn_raster.setData(spike_x, spike_y, connect='pairs', pen=default_pen)
        pn_plot.addItem(pn_raster)
        # LCA KC spike raster, MCA KC spike raster
        lca_kcs = int(config['kc']['number'] * config['kc']['lca_frac'] + 0.5)        
        lca_nodes = [str(ii) for ii in range(lca_kcs)]
        mca_nodes = [str(ii) for ii in range(lca_kcs, config['kc']['number'])]
        try:
            lca_spike_x, lca_spike_y = get_event_times(fd['/data/event/kc/kc_spiketime'],
                                               nodes=lca_nodes)
            mca_spike_x, mca_spike_y = get_event_times(fd['/data/event/kc/kc_spiketime'],
                                               nodes=mca_nodes)
        except KeyError:
            dirname = os.path.dirname(fname)
            fname = 'kc_spikes_' + os.path.basename(fname)
            with h5.File(os.path.join(dirname, fname)) as kc_file:
                lca_spike_x, lca_spike_y = get_event_times(kc_file, nodes=lca_nodes)
                mca_spike_x, mca_spike_y = get_event_times(kc_file, nodes=mca_nodes)
        print('LCA spikes', lca_spike_x.shape[0] // 2)
        gw.nextRow()
        kc_lca_plot = gw.addPlot(title='KC LCA')
        kc_lca_plot.setXLink(pn_plot)
        kc_lca_raster = pg.PlotCurveItem()
        kc_lca_raster.setData(lca_spike_x, lca_spike_y, connect='pairs', pen=default_pen)
        kc_lca_plot.addItem(kc_lca_raster)
        gw.nextRow()
        kc_mca_plot = gw.addPlot(title='KC MCA')
        kc_mca_plot.setXLink(pn_plot)
        kc_mca_raster = pg.PlotCurveItem()
        kc_mca_raster.setData(mca_spike_x, mca_spike_y, connect='pairs', pen=default_pen)
        kc_mca_plot.addItem(kc_mca_raster)
        # LCA KC Vm
        kc_vm_node = fd['/data/uniform/kc/KC_Vm']
        t = np.arange(kc_vm_node.shape[1]) * kc_vm_node.attrs['dt']
        kc_lca_vm = [kc_vm_node[int(ii), :]
                     for ii in np.random.choice(lca_nodes, size=vm_samples,
                                                replace=False)]
        gw.nextRow()
        kc_lca_vm_plot = gw.addPlot(title='KC LCA')
        kc_lca_vm_plot.setXLink(pn_plot)
        for kc_vm in kc_lca_vm:
            kc_vm_curve = kc_lca_vm_plot.plot(t, kc_vm, pen=default_pen)
        # MCA KC Vm
        kc_mca_vm = [kc_vm_node[int(ii), :]
                     for ii in np.random.choice(mca_nodes, size=vm_samples,
                                                replace=False)]
        gw.nextRow()
        kc_mca_vm_plot = gw.addPlot(title='KC MCA')
        kc_mca_vm_plot.setXLink(pn_plot)
        for kc_vm in kc_mca_vm:
            kc_vm_curve = kc_mca_vm_plot.plot(t, kc_vm, pen=default_pen)
        # GGN MCA Vm, GGN LCA Vm
        gw.nextRow()
        ggn_output_vm_plot = gw.addPlot(title='GGN CA')
        ggn_output_vm_plot.setXLink(pn_plot)
        ggn_output_vm_plot.addLegend()
        ggn_output_vm = fd['/data/uniform/ggn_output/GGN_output_Vm']
        pen_lca = (255, 0, 0, 100)
        pen_mca =  (0, 0, 255, 100)
        for ii in np.random.choice(range(ggn_output_vm.shape[0]), size=vm_samples, replace=False):
            sec = ggn_output_vm.dims[0]['source'][ii]            
            sid = sec.rpartition('dend_')[-1].partition('[')[0]
            if sid == '5':
                pen = pen_lca
            else:
                pen = pen_mca
            ggn_output_vm_plot.plot(t, ggn_output_vm[ii, :], pen=pen)
        ggn_output_vm_plot.legend.addItem(pg.PlotDataItem(pen=pen_lca), 'LCA')
        ggn_output_vm_plot.legend.addItem(pg.PlotDataItem(pen=pen_mca), 'MCA')
        # GGN alphaL Vm
        gw.nextRow()
        ggn_alphaL_vm_plot = gw.addPlot(title='GGN alphaL')
        ggn_alphaL_vm_plot.setXLink(pn_plot)
        ggn_alphaL_vm = fd['/data/uniform/ggn_alphaL_input/GGN_alphaL_input_Vm']
        for ii in np.random.choice(range(ggn_alphaL_vm.shape[0]), size=vm_samples, replace=False):
            ggn_alphaL_vm_plot.plot(t, ggn_alphaL_vm[ii,:], pen=default_pen)
        end = timer()
        print('Time for plotting {}s'.format(end - start))
        return gw


if __name__ == '__main__':
    gw = plot_spike_rasters(sys.argv[1])
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
        
        
# 
# pn_kc_ggn_plot.py ends here
