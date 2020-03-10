# single_syn_input_spatial_effect.py --- 
# 
# Filename: single_input_spatial_effect.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Jul 26 16:25:53 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Mon Sep 10 16:22:32 2018 (-0400)
#           By: Subhasis Ray
#     Update #: 609
# URL: 
# Doc URL: 
# Keywords: 
# Compatibility: 
# 
# 

# Commentary: 
# 
# 
# 
# 

# Code:
"""This script computes steady state voltage attenuation through GGN arbor."""
from __future__ import print_function
import os
import sys
import errno
from datetime import datetime
import h5py as h5
import argparse
import numpy as np
from matplotlib import pyplot as plt

#os.environ['NEURONHOME'] = 'c:\\nrn'
#sys.path += ['c:\\nrn\\lib\\python', 'd:\\subhasis_ggn\\model\\common']
from config import Q_, h, logger, timestamp, mypid, myjobid, nrn_version
import nrnutils as nu
import ephys
import neurograph as ng


def make_parser():
    parser = argparse.ArgumentParser(description=""" 
    simulate synaptic input via alpha synapse in `synnode`. The
    parameters `gmax` (uS), `tau` (ms) determine the properties of the
    synapse and `onset` (ms) determines the time for onset of the
    synaptic input. If specified, GGN_B99 template is read from the
    `celltemplate` file.""")
    parser.add_argument('-s', '--synnode', type=str, help='Section no. having voltage clamp',
                        dest='synnode',
                        default='soma')
    parser.add_argument('-f', '--celltemplate', type=str, help='cell template file',
                        dest='celltemplate', default='D:/subhasis_ggn/model/mb/cell_templates/GGN_20170309_sc.hoc')    
    parser.add_argument('-c', '--cell', type=str, help='cell name',
                        dest='cellname', default='GGN_20170309_sc')
    parser.add_argument('-o', '--onset', type=float, help='onset time (ms) for synaptic input',
                        dest='onset',
                        default=50.0)
    parser.add_argument('-n', '--reccount', type=int, help='number of compartments to record from',
                        dest='reccount', default=-1)
    parser.add_argument('-r', '--RM', type=float, help='specific membrane resistance (ohm-cm2)',
                        dest='RM', default=33e3) # 33 Kohm-cm2 according to Laurent, 1990
    parser.add_argument('-a', '--RA', type=float, help='specific axial resistance (ohm-cm)',
                        dest='RA', default=70.0)
    parser.add_argument('-v', '--vclamp', type=float, help='clamping voltage',
                        dest='vclamp', default=70.0)
    parser.add_argument('-e', '--Em', type=float, help='Leak reversal potential (mV)',
                        dest='Em', default=-51.0)
    parser.add_argument('--simtime', type=float,
                        help='simulation time (ms)',
                        dest='simtime', default=500)
    parser.add_argument('--deflection', action='store_true')
    parser.add_argument('--interactive', action='store_true')
    parser.add_argument('-d', '--datadir', type=str, default='.',
                        help='Data directory. The files will saved in this.'
                        )
    parser.add_argument('--outfile', type=str, 
                        help='output file name', default='')
    parser.add_argument('--rec-by-region', action='store_true', help='record `reccount` sections in every region if possible', dest='recbyreg')

    return parser




if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()

    mechs = {'pas': {'g': 1.0/Q_(args.RM, 'ohm*cm**2'), 'e': Q_(args.Em, 'mV')}}
    h.xopen(args.celltemplate)
    ggn = nu.create_cell(args.cellname, mechparams=mechs, Ra=args.RA, filename=args.celltemplate)
    g0 = nu.nrngraph(ggn)
    # for node in g0:
    #     print(node)
    # g, nmap = ng.renumber_nodes(g0, ggn.soma.name())
    reccount = args.reccount
    good_nodes = [n for n in g0.nodes() if g0.node[n]['orig'] is not None]
    ggn.geom_nseg()  # recompute the compartmentalization
    if (reccount <= 0) or (reccount >= len(good_nodes)):
        reccount = len(good_nodes)
        print('Number of sections to record from exceeds total number. Using total:', reccount)
    nodes = [n for n in g0.nodes() if g0.node[n]['orig'] is not None]
    if args.recbyreg:
        rec_regions =  [1, 5, 6, 7, 8]
        rec_nodes_by_sid = ng.select_random_nodes_by_sid(g0.subgraph(nodes), rec_regions, [reccount] * len(rec_regions))
        rec_nodes = list(np.concatenate(rec_nodes_by_sid.values()))
    else:
        rec_nodes = list(np.random.choice(nodes, size=reccount, replace=False))
    recsecs = [g0.node[n]['orig'] for n in rec_nodes]
    synnode = '{}.{}'.format(list(g0.nodes)[0].partition('.')[0],
                             args.synnode)
    if synnode not in g0:
        synnode = np.random.choice(g0.nodes(), size=1)[0]
    # print(synnode)
    synsec = g0.node[synnode]['orig']
    while synsec is None:
        print('No section for', synnode, '.Trying a random one')
        synnode = np.random.choice(g0.nodes(), size=1)
        synsec = g0.node[synnode]['orig']        
    
    if synnode not in rec_nodes:
        recsecs.append(synsec)
        rec_nodes.append(synnode)

    synidx = rec_nodes.index(synnode)
        
    # alpha_syn = insert_alphasynapse(synsec(0.5), onset=args.onset, gmax=args.gsyn, tau=args.tau)
    clamp = ephys.setup_voltage_clamp(synsec,
                                      Q_('{}mV'.format(args.Em)),
                                      Q_('{}ms'.format(args.onset)),
                                      Q_('{}mV'.format(args.Em)),
                                      Q_('0ms'),
                                      Q_('{}mV'.format(args.vclamp)),
                                      Q_('{}ms'.format(args.simtime -
                                                       args.onset)),
                                      pos=0.5,
                                      rs=Q_('0.01ohm'))
    t_vec = h.Vector(np.arange(0, args.simtime, h.dt * 10))  # record data every 10 timesteps
    tabs = [nu.setup_recording(sec, 0.5, 'v', t_vec) for sec in recsecs]
    ts = datetime.now()
    # Initialize and run
    h.tstop = args.simtime
    print('Setting Vm to {}'.format(args.Em))
    h.finitialize(args.Em)
    h.fcurrent()
    while h.t < h.tstop:
        h.fadvance()
    # h.init()
    # h.run()
    
    te = datetime.now()
    delta = te - ts
    print('Time for', h.tstop*1e-3, 's simulation =', delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds)
    if args.interactive:
        plt.plot(t_vec.x, tabs[synidx].x)
        plt.show()
    try:
        os.makedirs(args.datadir)
    except OSError as exc:        
        if exc.errno == errno.EEXIST and os.path.isdir(args.datadir):
            pass
        else:
            raise
    fname = args.outfile
    if fname == '':
        fname = 'Vm_vclamp_{}_JID{}_PID{}.h5'.format(ts.strftime('%Y_%m_%d__%H_%M_%S'),
                                               myjobid,
                                               mypid)
    outfilepath = os.path.join(args.datadir, fname)                               
    # ggn_graph_undirected = g0.to_undirected()
    with h5.File(outfilepath) as fd:
        syninfo = fd.create_group('vclamp')
        syninfo.attrs['section'] = synsec.name()
        syninfo.attrs['clamp_voltage'] = args.vclamp
        syn_region = synsec.name().rpartition('.')[-1].rpartition('_')[-1].partition('[')[0]
        syninfo.attrs['synregion'] = syn_region
        syninfo.attrs['node'] = synnode
        syninfo.attrs['onset'] = args.onset
        fd.attrs['celltemplatefile'] = args.celltemplate
        fd.attrs['RA'] = args.RA
        fd.attrs['g_pas'] = mechs['pas']['g'].to('S/cm**2').m
        fd.attrs['e_pas'] = mechs['pas']['e'].to('mV').m
        fd.attrs['simtime'] = args.simtime
        with open(args.celltemplate, 'r') as cellfile:
            celldef = '\n'.join(cellfile.readlines())
            fd.create_dataset('celltemplate', data=np.string_(celldef))
        if args.deflection:
            secname_array = []
            vclamp_array = []
            for v, recsec, recnode in zip(tabs, recsecs, rec_nodes):
                vclamp_array.append(np.array(v).max())
                secname_array.append(recsec.name())
            ds = fd.create_dataset('v_deflection', data=vclamp_array)
            secs = fd.create_dataset('section', data=secname_array, dtype=h5.special_dtype(vlen=bytes))
            ds.dims[0].label = 'section'
            ds.dims.create_scale(secs, 'section')
            ds.dims[0].attach_scale(secs)
        else:
            tdata = np.array(t_vec)
            time = fd.create_dataset('time', data=tdata)
            for v, recsec, recnode in zip(tabs, recsecs, rec_nodes):
                vdata = np.array(v)

                ds = fd.create_dataset('v_{}'.format(recnode), data=vdata)
                ds.attrs['section'] = recsec.name()
                post_stim = tdata.searchsorted(args.onset)
                rec_region = recsec.name().rpartition('.')[-1].rsplit('_')[-1].split('[')[0]
                ds.attrs['region'] = rec_region
    if args.interactive:
        plt.plot(t_vec.x, tabs[-1].x)
        plt.show()
    print('Saved data in', outfilepath)
    


# 
# single_syn_input_spatial_effect.py ends here
