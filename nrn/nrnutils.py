# nrnutils.py --- 
# 
# Filename: nrnutils.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Tue Jul 26 16:27:51 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Tue Sep 25 12:37:22 2018 (-0400)
#           By: Subhasis  Ray
#     Update #: 527
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

from __future__ import print_function

import sys
from datetime import datetime
from collections import defaultdict
import numpy as np
import networkx as nx
import h5py as h5
from timeit import default_timer as timer
from config import ur, logger, Q_
import ephys
import neurograph as ng
from neurograph import (tograph, toswc, sorted_edges, branch_points, remove_null_edges, eucd, renumber_nodes)
# from morphoplot import (plot_3d_lines, plot_nodes)


from neuron import h

print('#C')
sys.stdout.flush()

# plt.rcParams['axes.facecolor'] = 'black'
# plt.style.use('dark_background')
h.load_file("stdrun.hoc")

print('#D')
sys.stdout.flush()


def sectype(cid):
    """Get the SWC structure type (sid) from section name `cid`"""
    if '.dend[' in cid:
        stype = 3
    elif '.dend_' in cid:
        stype = int(cid.split('_')[-1].split('[')[0])
    elif 'soma' in cid:
        stype = 1
    elif 'axon' in cid:
        stype = 2
    else:
        stype = 0
    return stype


def nrngraph(cell):
    """Convert a NEURON cell model into networkx graph. The nodes are the
    section names in the cell.

    Each node has the parent section name in attribute `p`. It is -1
    if there is no parent.  x, y, z represent the coordinates of the 0
    end.  orig is SectionRef to the original NEURON section.

    For display purposes, we add the 1 end of the leaf sections as
    dummy nodes. They can be identified by the attribute orig=None.

    """
    g = nx.DiGraph()
    cell.soma.push()   # Assume there is a designated soma
    # This appends all sections accessible from currently
    # accessed section in root to leaves
    ordered_tree = h.SectionList()
    ordered_tree.wholetree()
    h.pop_section()
    for comp in ordered_tree:
        comp.push()
        cid = comp.name()
        stype = sectype(cid)
        g.add_node(cid, x=h.x3d(0),
                   y=h.y3d(0), z=h.z3d(0),
                   r=h.diam3d(0)/2.0, s=stype, orig=comp)
        ref = h.SectionRef(sec=comp)
        if ref.has_parent():
            parent = ref.parent
            g.add_edge(parent.name(), cid, length=parent.L)
            g.node[cid]['p'] = parent.name()
        else:
            g.node[cid]['p'] = -1
        # Insert the 1-end of a section as a dummy node
        if ref.nchild() == 0:
            leaf_node = '{}_1'.format(cid)
            g.add_node(leaf_node, x=h.x3d(1),
                       y=h.y3d(1), z=h.z3d(1),
                       r=h.diam3d(1)/2.0, s=stype, p=comp.name(), orig=None)
            g.add_edge(cid, leaf_node, length=comp.L)
        h.pop_section()
    return g


def get_section_node_map(g):
    ret = {}
    for n, d in g.nodes_iter(data=True):
        if d['orig'] is None:
            continue
        ret[d['orig'].name()] = n
    return ret


def get_alen_pos3d(sec):
    """Get the arclength and 3D position of the poinst in section sec.
    Inspired by
    http://www.neuron.yale.edu/ftp/ted/neuron/extracellular_stim_and_rec.zip:
    interpxyz.hoc

    Returns: length, pos where length is the list of arc3d lengths and
    pos the list of 3D positions (x, y, z), of the 3D points in
    section `sec`.

    """
    npts = int(h.n3d(sec=sec))
    pos = []
    length = []
    for ii in range(npts):
        pos.append((h.x3d(ii, sec=sec),
                    h.y3d(ii, sec=sec),
                    h.z3d(ii, sec=sec)))
        length.append(h.arc3d(ii, sec=sec))
    return length, pos


def get_seg_3dpos(sec, xlist):
    """Obtain the nearest 3D point position on or before the segments
    specified by 1D position in xlist."""
    length, pos = get_alen_pos3d(sec)
    length = np.array(length) / length[-1]  # normalize the length
    ret = []
    for xx in xlist:
        ii = np.searchsorted(length, xx)
        ret.append(pos[ii])
    return ret

                   
def select_good_nodes_by_sid(g, sid_list, counts, replace=False):
    """For each sid in `sid_list` select `count` random nodes with an
    underlying `section` from `g` - a neurongraph.

    Returns a list of selected nodes.

    @seealso: neurograp.select_random_nodes_by_sid

    """
    good_nodes = [n for n, data in g.nodes(data=True) if data['orig'] is not None]
    type_node_map = ng.get_stype_node_map(g.subgraph(good_nodes))
    synnodes = []
    for branch_id, count in zip(sid_list, counts):
        size = len(type_node_map[branch_id])
        if (count > size) and (replace == False):
            print('Changing number of nodes to maximum {} available in branch, since replace={}'.format(size, replace))
            count = size
        synnodes += list(np.random.choice(type_node_map[branch_id],
                                          size=count, replace=replace))
    return synnodes


def select_random_segments_by_sid(g, sid, count, by_length=True, replace=True):
    """Select segments from sections with specified sid.  If by_length is
    True, select with probability proportional to the length of the
    segment.

    """
    good_secs = [data['orig'] for n, data in g.nodes(data=True)
                 if (g.node[n]['orig'] is not None)
                 and (g.node[n]['s'] == sid)]
    seg_list = []
    seg_lengths = []
    for sec in good_secs:
        for seg in sec:
            seg_list.append(seg)
            seg_lengths.append(sec.L/sec.nseg)
    seg_lengths = np.array(seg_lengths)
    probability = None
    if by_length:
        probability = seg_lengths / np.sum(seg_lengths)
    segs = np.random.choice(seg_list, size=count, p=probability, replace=replace)
    # Don't keep the segments around beyond this function, x and sec
    # will suffice to retrieve them
    return [( seg.sec, seg.x) for seg in segs]


def select_random_terminal_segments_by_sid(g, sid, count, by_length=True, replace=True):
    """Select segments from sections with specified sid.  If by_length is
    True, select with probability proportional to the length of the
    segment.

    """
    good_secs = []
    for n in g.nodes():
        sec = g.node[n]['orig']        
        if sec is not None and (g.node[n]['s'] == sid):
            ref = h.SectionRef(sec=sec)
            if len(ref.child) == 0:
                good_secs.append(sec)
    seg_list = []
    seg_lengths = []
    for sec in good_secs:
        for seg in sec:
            seg_list.append(seg)
            seg_lengths.append(sec.L/sec.nseg)
    seg_lengths = np.array(seg_lengths)
    probability = None
    if by_length:
        probability = seg_lengths / np.sum(seg_lengths)
    segs = np.random.choice(seg_list, size=count, p=probability, replace=replace)
    # Don't keep the segments around beyond this function, x and sec
    # will suffice to retrieve them
    return [( seg.sec, seg.x) for seg in segs]  


# h.load_file("nrn_run.hoc")
# h.load_file('cells/ggn_passive.hoc')
# print('Loaded cell')
# dend_8 = h.dend_8 # The sections become objects inside hoc

## Retrieve the A1 and A2 tables for type IIa interneuron
# h.clipboard_retrieve('A1_IIa.csv')
# a1 = h.hoc_obj_[0]
# v1 = h.hoc_obj_[1]
# h.table_A1_caiia(a1._ref_x[0], a1.size(), v1._ref_x[0])
# a2 = h.hoc_obj_[0]
# v2 = h.hoc_obj_[1]
# h.table_A2_caiia(a2._ref_x[0], a2.size(), v2._ref_x[0])

# h.xopen('ggn_b99_corrected_passive.hoc')
# h.xopen('GGN_B99_template_20160719.hoc')
# h.xopen('GGN_B99_template_v0.hoc')


def insert_alphasynapse(segment, onset=200.0, gmax=1e-3, tau=0.1):
    asyn = h.AlphaSynapse(segment)
    asyn.onset = onset
    asyn.gmax = gmax
    asyn.tau = tau
    return asyn


def insert_exp2syn(seglist, tau1, tau2, e):
    """Insert double exponential synapse at the specified segments.

    Returns a list of created synapses. Note that the actual
    conductance is determined from the `weight` of the associated
    NetCon object.

    """
    synlist = []
    for seg in seglist:
        syn = h.Exp2Syn(seg)
        syn.e = e
        syn.tau1 = tau1
        syn.tau2 = tau2
        synlist.append(syn)
    return synlist


def make_exp2syn_one_one(presec, prepos, postsec, postpos, syn_params):
    """Make identical synapses between lists of sections and positions.

    presec: list of presynaptic sections.

    prepos: list of positions in presynaptic sections.

    postsec: list of postsynaptic sections.

    postpos: list of postsynaptic positions.

    syn_params: properties of the synapse. A dict containing tau1,
    tau2, e, threshold, delay, gmax for the synapse.

    Returns

    dict containing `synapses` - a list of synapses,
    `netcons` a list of netcon objects
    `adjacency` a list of tuple (name of presynaptic section, 
                                 presynaptic position, 
                                 name of postsynaptic section, 
                                 postsynaptic position)
    """
    adjacency = []
    synapses  = []
    netcons = []
    tau1 = syn_params.get('tau1')
    tau2 = syn_params.get('tau2')
    e =syn_params.get('e')
    thresh = syn_params.get('threshold')
    delay = syn_params.get('delay')
    gmax = syn_params.get('gmax')
    for presec_, prepos_, postsec_, postpos_ in zip(presec, prepos,
                                                    postsec, postpos):
        syn = h.Exp2Syn(postsec_(postpos_))
        syn.tau1, syn.tau2, syn.e = tau1, tau2, e
        nc = h.NetCon(presec_(prepos_)._ref_v, syn, thresh, delay,
                      gmax, sec=presec_)
        synapses.append(syn)
        netcons.append(nc)
        adjacency.append((presec_.name(), prepos_, postsec_.name(), postpos_))
    return {'synapses': synapses, 'netcons': netcons,
            'adjacency': adjacency}


def make_gradedsyn_one_one(presec, prepos, postsec, postpos, syn_params):
    """Make identical synapses between lists of sections and positions.

    presec: list of presynaptic sections.

    prepos: list of positions in presynaptic sections.

    postsec: list of postsynaptic sections.

    postpos: list of postsynaptic positions.

    syn_params: properties of the synapse. A dict containing the
    attributes for the synapse.

    Returns

    dict containing `synapses` - a list of synapses,
    `adjacency` a list of tuple (name of presynaptic section, 
                                 presynaptic position, 
                                 name of postsynaptic section, 
                                 postsynaptic position)

    """
    adjacency = []
    synapses  = []
    for presec_, prepos_, postsec_, postpos_ in zip(presec, prepos,
                                                    postsec, postpos):
        syn = h.GradedSyn(postsec_(postpos_))
        for attr, val in syn_params.items():
            setattr(syn, attr, val)
        h.setpointer(presec_(prepos_)._ref_v, 'vpre', syn)
        synapses.append(syn)
        adjacency.append((presec_.name(), prepos_, postsec_.name(), postpos_))
    return {'synapses': synapses,
            'adjacency': adjacency}



def connect_stim_syn_one_one(spiketrains, synlist, threshold=-20, delay=0, weight=1):
    """Create a list of VecStim object with event times from spiketrains,
    and create NetCon objects to feed this to each synapse. The
    connections are created one-to-one.

    spiketrains: list of arrays containing stimulus time for each entry in synlist

    synlist: list of synapses

    Returns list of NetCon objects, VecStim object and the list of
    Vectors containing stimulus times. The maximum conductance of the
    synapses are determined by the weight parameter. (See
    NetCon.weight docs)

    """
    netcons = []
    stimtimes = []
    stimlist = []
    if not isinstance(threshold, collections.Iterable):
        threhsold = [threshold] * len(spiketrains)
    if not isinstance(delay, collections.Iterable):
        delay = [delay] * len(spiketrains)
    if not isinstance(weight, collections.Iterable):
        weight = [weight] * len(spiketrains)
    
    for st, syn, th, d, wt in zip(spiketrains, synlist, threhsold,
                                      delay, weight):
        stimvec = h.Vector()
        # stimvec.append(*st)
        stimvec = stimvec.from_python(st)
        stimtimes.append(stimvec)
        vecstim = h.VecStim()
        vecstim.play(stimvec)
        stimlist.append(vecstim)
        netcons.append(h.NetCon(vecstim, syn, th, d, wt))
    return netcons, stimlist, stimtimes


def create_cell(name, mechparams=None, Ra=None, calc_nseg=False, filename=None):
    """Make a copy of the GGN and insert mechanisms based on mechparams.

    mechparams: dict of dicts containing parameters for
    mechanisms. Each entry should be
    
    name: { param: value, param: value, ...}
    
    Other than name, the param keys depend on the mechanism
    definition. For example, the predefined passive mechanism has name
    `pas` and it defines g and e as parameters. The corresponding dict entry
    will be `'pas': {'g': 2.5e-4, 'e'=-65.0}` for inserting
    passive conductance with density 2.5e-4 S/cm2 on each compartment
    and the reversal potential -65 mV.

    """
    if not hasattr(h, name):
        if filename is None:
            raise Exception('Cell not preloaded. Specify template filename')
        h.xopen(filename)
    cell = eval('h.{}()'.format(name))
    if mechparams is None:
        mechparams = {}
    print(mechparams)
    mech_dict = ephys.create_mechs(mechparams)
    insert_mechanisms(cell, mech_dict.values())
    if (Ra is not None) or calc_nseg:
        cell.soma.push()
        ordered_tree = h.SectionList()
        ordered_tree.wholetree()
        h.pop_section()
        if not isinstance(Ra, float ):  # Assume pint.UnitRegistry().Quantity
            Ra = ephys.Q_(Ra).to('ohm*cm').m
        for sec in ordered_tree:
            sec.push()
            if Ra is not None:
                sec.Ra = Ra
            # lambda_f = np.sqrt(sec.diam/(4 * sec.Ra * sec.g_pas))
            # nseg = int(10 * sec.L / lambda_dc + 0.5)
            nseg = int((sec.L/(0.1 * h.lambda_f(100)) + 0.999)/2) * 2 + 1
            sec.nseg = nseg
            h.pop_section()
    return cell


def setup_dynamic_clamp(sec, pos, vm, t, dt=h.dt):
    """Create a dynamic clamp on section `sec` at position `pos` and play
    the voltage `vm` over time `t`. The sample points are selected at
    `dt` interval.

    """
    clamp = h.SEClamp(sec(pos))
    clamp.rs = 1e-3
    clamp.dur1 = 1e9
    vsample = np.interp(np.arange(t[0], t[-1], dt), t, vm)
    vplay = h.Vector(vsample)
    vplay.play(clamp._ref_amp1, dt)
    return clamp, vplay


def record_vm(sec, pos=0.5, t=None):
    vec = h.Vector()
    if t is None:
        vec.record(sec(pos)._ref_v)
    else:
        vec.record(sec(pos)._ref_v, t)        
    return vec


def record_spiketimes(sec, pos=0.5, threshold=-20.0):
    if not hasattr(h, 'nil'):
        h('objref nil')
    sec.push()
    nc = h.NetCon(sec(0.5)._ref_v, h.nil, threshold, 0.0, 1.0)
    h.pop_section()
    vec = h.Vector()
    nc.record(vec)
    return nc, vec                  


def setup_recording(sec, pos=0.5, field='v', t=None):
    """Convenience function for recording an arbitrary field (`field`)
    from section `sec` at position `pos`.
    
    If `t` is None, recording is done every integration timestep
    (h.dt), if a float, this is the time interval between recordings,
    if another vector, recording is done at the time points in `t`
    (see neuron Vector.record documentation).

    """
    vec = h.Vector()
    if t is None:
        vec.record(getattr(sec(pos), '_ref_{}'.format(field)))
    else:
        vec.record(getattr(sec(pos), '_ref_{}'.format(field)), t)
    return vec


def setup_current_clamp(cable, delay=100.0, duration=10.0, amplitude=100e-3,
                        pos=0.0):
    """Insert a current clamp electrode to deliver `amplitude` nA current
    at `pos` fraction of length for `duration` ms starting at `delay`
    ms"""
    print(cable, type(cable))
    trode = h.IClamp(pos, sec=cable)
    setattr(trode, 'del', delay)
    trode.dur = duration
    trode.amp = amplitude
    print('Inserted current clamp electrode at', pos,
          'fraction of cable to deliver', amplitude,
          'nA current starting at', delay,
          'ms for', duration, 'ms')
    return trode


def insert_mechanisms(cell, mechanisms, ek=-80.0 * ur.mV,
                      ena=55.0 * ur.mV, eca=160.0 * ur.mV):
    """Insert a mechanism in all sections for a given conductance density"""
    for sec in cell.allsec():
        for mech in mechanisms:
            mech.insert_into(sec)

            
def block_run(Dt=100.0, logger=None):
    """Run for tstop in blocks of Dt time"""
    if h.tstop < Dt:
        Dt = h.tstop
    if logger is not None:
        logger.info('Starting simulation for {}'.format(h.tstop))
    print('Starting simulation for {}'.format(h.tstop))
    start = timer()
    while h.t < (h.tstop - Dt):
        h.continuerun(h.t+Dt)
        end = timer()
        if logger is not None:
            logger.info('Finished till {} ms of simulation in {} s.'.format(
                h.t, end - start))
    h.continuerun(h.tstop)
    end = timer()
    if logger is not None:
        logger.info('Finished {} ms of simulation in {} s.'.format(
            h.tstop, end - start))

            
# 
# nrnutils.py ends here
