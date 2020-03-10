# staggered_input.py ---
#
# Filename: sequential_input.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Fri, Nov 18, 2016  5:14:40 PM

# Version:
# Last-Updated: Tue Mar 10 14:51:31 2020 (-0400)
#           By: Subhasis Ray
#     Update #: 530
# Commentary:
#
#
#
#

# Code:
"""Simulate random synaptic input using inhomogeneous Poisson process
at a single area on a GGN model.

Uses double-exponential synapse - by default the two time constants
are same and it has an alpha function form.

The output data is saved as HDF5 file containing
/celltemplate - the NEURON celltemplate definition (dataset containing a string)
/syninfo - dataset containing the section names where synapses were created.
           tau1, tau2, maxrate and gmax are stored as attributes.

if `--deflection` is passed in command line
/v_deflection - two column dataset ('time', 'deflection') containing the
        time of peak voltage and the peak voltage deflection from Em for
        each section.
/section - dataset containing section names corresponding to rows in /v_deflection.
        This is attached as a dimension scale to dimension 0 of /v_deflection

if `--deflection` is NOT specified in command line
/time - time points when voltages were recorded
/v_{section_name} - dataset containing the recorded voltage in section `section_name`.

"""
from __future__ import print_function
import sys
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

plt.rc('figure', figsize=(16, 10))
plt.rc('font', size=32)


def connect_independent_stimuli(spiketrains, synlist, threshold=-20, delay=0, weight=1):
    """Create a list of VecStim object with event times from spiketrains,
    and create NetCon objects to feed this to each synapse.

    spiketrains: list of arrays containing stimulus time for each entry in synlist

    synlist: list of synapses

    Returns list of NetCon objects, VecStim object and the list of Vectors
    containing stimulus times.

    """
    netcons = []
    stimtimes = []
    stimlist = []
    for st, syn in zip(spiketrains, synlist):
        stimvec = h.Vector()
        stimvec.append(*st)
        stimtimes.append(stimvec)
        vecstim = h.VecStim()
        vecstim.play(stimvec)
        stimlist.append(vecstim)
        netcons.append(h.NetCon(vecstim, syn, threshold, delay, weight))
    return netcons, stimlist, stimtimes


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    simulate synaptic input via a synapse in `syncomp`. The parameters
    `gmax` (uS), `tau1` (ms) and `tau2` (ms) determine the properties
    of the synapse and `maxrate` (Hz) determines the global rate of
    the synaptic input. Vm from `reccount` unique compartments from
    each major branch is reccorded. If specified, GGN_B99 template is
    read from the `celltemplate` file.""")
    parser.add_argument('-i', '--input', type=int,
                        help='custom type no. identifying the branch where '
                        'synaptic input should be applied. These are, '
                        '1: soma, 3: basal dendrite, 4: apical dendrite, '
                        '5: lateral calyx, 6: medial calyx, 7: lateral horn, '
                        '8: alpha lobe',
                        dest='inputs',
                        action='append',
                        required=True)
    parser.add_argument('-g', type=float,
                        help='peak conductance of the synapses (in uS)',
                        dest='gmax',
                        default=1e-3)
    parser.add_argument('--Em', type=float,
                        help='Leak reversal potential (mV)',
                        dest='Em',
                        default=-51.0)
    parser.add_argument('--gpas', type=float,
                        help='specific membrane conductance (in S/cm2)',
                        dest='gpas',
                        default=1/33e3)
    parser.add_argument('--ra', type=float,
                        help='specific axial resistance (in Ohm-cm)',
                        dest='RA',
                        default=100.0)
    parser.add_argument('-r', '--rise-t', type=float,
                        help='rise time constants for the synapse (in ms)',
                        dest='rise_t', default=1.0)
    parser.add_argument('-d', '--decay-t', type=float,
                        help='decay time constants for the synapse (in ms)',
                        dest='decay_t', default=1.0)
    parser.add_argument('-s', '--onset', type=float,
                        help='start time of synaptic input (in ms).',
                        dest='onset',
                        default=50)
    parser.add_argument('-n', type=int,
                        help='number of synapses on each input region',
                        dest='syncount', default=10)
    parser.add_argument('-c', type=str,
                        help='cell name',
                        dest='cellname',
                        default='GGN_B99_20161205_d_lambda.hoc')
    parser.add_argument('-f', type=str,
                        help='cell template file',
                        dest='celltemplate',
                        default='GGN_B99_20161205_d_lambda.hoc')
    parser.add_argument('--colormap', type=str,
                        help='Colormap to use for plotting and morphology '
                        'display: {}'.format(ng.colormaps.keys()),
                        dest='colormap',
                        default='5cs3')
    parser.add_argument('--dia-scale', type=float,
                        help='scale diameter by this factor if it is less than'
                        ' `dia-lim`',
                        default=None,
                        dest='diascale')
    parser.add_argument('--deflection', action='store_true', dest='delv',
                        help='Save only peak deflection in Vm')
    parser.add_argument('--dia-lim', type=float,
                        help='minimum diameter over which it should be scaled'
                        ' up (um)',
                        default=None,
                        dest='dialim')
    parser.add_argument('--maxrate', type=float,
                        help='Maximum rate (in Hz) of input on each synapse',
                        default=100, dest='maxrate')
    parser.add_argument('--simtime', type=float,
                        help='simulation time (ms)',
                        dest='simtime', default=500.0)
    parser.add_argument('-w', '--dur', type=float,
                        help='stimulation duration/width (ms)',
                        dest='dur', default=500.0)
    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        help='output filename')
    parser.add_argument('--rec', type=int,
                        help='number of sections to record from in each '
                        'region. <= 0 means record all sections',
                        dest='rec', default=0)

    args = parser.parse_args()
    inputs = args.inputs
    celltemplate = args.celltemplate
    syncount = args.syncount
    onset = args.onset
    gmax = args.gmax
    diascale = args.diascale
    dialim = args.dialim
    maxrate = args.maxrate / 1000.0  # s^-1 to per ms^-1
    mechs = {'pas': {'g': '{} S/cm**2'.format(args.gpas), 'e': '-51.0mV'}}
    h.xopen(celltemplate)
    ggn = nu.create_cell(args.cellname, mechparams=mechs,
                         Ra='{}ohm*cm'.format(args.RA),
                         filename=args.celltemplate)
    if (args.diascale is not None) and (args.dialim is not None):
        count = 0
        ggn.soma.push()
        ordered_tree = h.SectionList()
        ordered_tree.wholetree()
        h.pop_section()
        for sec in ordered_tree:
            sec.push()
            for ii, seg in enumerate(sec):
                if seg.diam < dialim:
                    seg.diam = seg.diam * diascale
                    count += 1
            h.pop_section()
        print('Scaled diameters of', count,
              'sections whose diameters were <', dialim, 'by', diascale)
    g0 = nu.nrngraph(ggn)
    g, nmap = ng.renumber_nodes(g0, ggn.soma.name())
    # This also gets the dummy nodes at the tips
    stype_node_map_all = ng.get_stype_node_map(g)
    stype_node_map = {}
    # Collect only nodes with sections associated.
    for stype, nodes in stype_node_map_all.items():
        good = [node for node in nodes if g.node[node]['orig'] is not None]
        stype_node_map[stype] = good
    synnodes = nu.select_good_nodes_by_sid(g, inputs,
                                           [syncount] * len(inputs),
                                           replace=True)
    synsecs = [g.nodes[n]['orig'] for n in synnodes]
    synsegs = [sec(1.0) for sec in synsecs]
    # Create the exp2syn
    # ereversal = 0 for glutamatergic synapse
    synlist = nu.insert_exp2syn(synsegs, args.rise_t, args.decay_t, 0.0)
    print('Number of synapses', syncount, 'synlist', len(synlist))
    # Create stimulation system
    t_stop = args.simtime
    delta_t = 10.0
    # Build a rate function which linearly goes down from maxrate at
    # t=0 to 0 at t = t_stop.
    rate_fn = lambda t: maxrate * (args.dur - t)/(args.dur)
    spiketimes_list = []
    for syn in synlist:
        spiketimes = nhpp_thinning(rate_fn, args.dur, delta_t)
        spiketimes_list.append(spiketimes + args.onset)
    netcons, stimlist, stimvecs = connect_independent_stimuli(spiketimes_list,
                                                              synlist,
                                                              threshold=-20.0,
                                                              delay=0.0,
                                                              weight=args.gmax)
    if args.rec <= 0:
        rec = [(n, data['orig']) for n, data in g.nodes(data=True)
               if data['orig'] is not None]
    else:
        rec = [(n, g.nodes[n]['orig']) for n in nu.select_good_nodes_by_sid(
            g, stype_node_map.keys(), [args.rec] * len(stype_node_map))]
    syn = zip(synnodes, synsecs)
    rec = list(set(rec + syn))
    recnodes, recsecs = zip(*rec)
    h.tstop = t_stop
    t_vec = h.Vector(np.arange(0, h.tstop, 10 * h.dt))
    # Record data
    tabs = [nu.setup_recording(sec, field='v', t=t_vec) for sec in recsecs]
    ts = datetime.now()
    h.finitialize(args.Em)
    h.fcurrent()
    while h.t < h.tstop:
        h.fadvance()
    # h.run()
    sys.stdout.flush()
    te = datetime.now()
    delta = te - ts
    print('Time for', h.tstop*1e-3, 's simulation =',
          delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds)
    if args.outfile == '':
        outfilename = 'data/Vm_inhpoisson_stim_series_{}.h5'.format(
            ts.strftime('%Y%m%d_%H%M%S'))
    else:
        outfilename = args.outfile
    with h5.File(outfilename) as fd:
        fd.attrs['argv'] = ' '.join(sys.argv)
        sections = [g.node[s]['orig'].name() for s in synnodes]
        syninfo = fd.create_dataset('syninfo', data=sections)
        syninfo.attrs['tau1'] = args.rise_t
        syninfo.attrs['tau2'] = args.decay_t
        syninfo.attrs['maxrate'] = args.maxrate
        syninfo.attrs['gmax'] = gmax
        fd.attrs['celltemplatefile'] = args.celltemplate
        fd.attrs['RA'] = '{}ohm*cm'.format(args.RA)
        fd.attrs['g_pas'] = mechs['pas']['g']
        fd.attrs['e_pas'] = mechs['pas']['e']
        fd.attrs['simtime'] = args.simtime
        fd.attrs['rec_loc'] = 1.0
        tarray = np.asarray(t_vec)
        if not args.delv:
            time = fd.create_dataset('time', data=tarray)
        stiminfo = fd.create_group('stimulus')
        with open(args.celltemplate, 'r') as cellfile:
            celldef = '\n'.join(cellfile.readlines())
            fd.create_dataset('celltemplate', data=np.string_(celldef))
        for ii, (vec, syn) in enumerate(zip(stimvecs, synnodes)):
            stim = stiminfo.create_dataset('stim_{}'.format(ii), data=vec.x)
            stim.attrs['section'] = g.node[syn]['orig'].name()
            stim.attrs['node'] = syn
            # print('Stimulus', ii, 'on', stim.attrs['section'])
        if args.delv:
            t_vpeak_array = []
            secname_array = []
            for v, node, sec in zip(tabs, recnodes, recsecs):
                varray = np.array(v.x)
                maxidx = np.argmax(varray)
                t_vpeak_array.append((tarray[maxidx], varray[maxidx]))
                secname_array.append(sec.name())
            ds = fd.create_dataset('v_deflection', data=np.array(
                t_vpeak_array, dtype=[('time', np.float),
                                      ('vm', np.float)]))
            section = fd.create_dataset('section', data=secname_array,
                                        dtype=h5.special_dtype(vlen=bytes))
            ds.dims[0].label = 'section'
            dscale = ds.dims.create_scale(section, 'section')
            ds.dims[0].attach_scale(section)
        else:
            for v, node, sec in zip(tabs, recnodes, recsecs):
                ds = fd.create_dataset('v_{}'.format(sec.name()), data=v.x)
                ds.attrs['section'] = sec.name()
                ds.attrs['node'] = node
                ds.attrs['sid'] = g.node[node]['s']

#
# staggered_input.py ends here
