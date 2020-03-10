"""Check the effect of synchronized synaptic input at random locations
of specific branches of the GGN.

"""
from __future__ import print_function
import sys
import argparse
from datetime import datetime
import h5py as h5
import nrnutils as nu
import neurograph as ng

def plot_Vm(fname):
    with h5.File(fname, 'r') as fd:
        syninfo = [s for s in fd['syninfo']]
        for nodename in fd:
            if nodename.startswith('v_'):
                segment = fd[nodename].attrs['section']
                if segment in syninfo:
                    color = 'gray'
                    alpha = 0.3
                    lw = 1.0
                    ls = '--'
                else:
                    color=None
                    alpha = 0.7
                    lw = 2.0
                    ls = '-'
                plt.plot(fd['time'], fd[nodename], color=color, alpha=alpha, ls=ls, lw=lw, label=segment)
    plt.legend()
    plt.show()        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    simulate synaptic input via alpha synapse in `syncomp`. The
    parameters `gmax` (uS), `tau` (ms) determine the properties of the
    synapse and `onset` (ms) determines the time for onset of the
    synaptic input. Vm from `reccount unique compartments from each
    major branch is reccorded. If specified, GGN_B99 template is read
    from the `celltemplate` file.""")
    parser.add_argument('-i', '--input', type=int,
                        help='custom type no. identifying the branch where synaptic input should be applied. These are, 1: soma, 3: basal dendrite, 4: apical dendrite, 5: lateral calyx, 6: medial calyx, 7: lateral horn, 8: alpha lobe',
                        dest='inputs',
                        action='append',
                        required=True)
    parser.add_argument('-g', type=float,
                        help='peak conductance of the synapses (in uS)',
                        dest='gmax',
                        default=1e-3)
    parser.add_argument('-t', type=float,
                        help='time constant of the alpha synapse (in ms)',
                        dest='tau', default=1.0)
    parser.add_argument('-s', type=float,
                        help='onset time of synaptic input (in ms)',
                        dest='onset', default=50.0)
    parser.add_argument('-n', type=int,
                        help='number of synapses',
                        dest='syncount', default=10)
    parser.add_argument('--synfile', type=str,
                        help='File containing list of synapse sections',
                        dest='synfile', default=None)
    parser.add_argument('--reccount', type=int,
                        help='number of sections to record from in each region',
                        dest='reccount',
                        default=10)
    parser.add_argument('--recfile', type=str,
                        help='File containing list of recording sections',
                        dest='recfile', default=None)
    parser.add_argument('-f', type=str,
                        help='cell template file',
                        dest='celltemplate',
                        default='GGN_B99_20160725.hoc')
    parser.add_argument('-c', type=str,
                        help='cell name',
                        dest='cellname',
                        default='GGN_B99')
    parser.add_argument('--dia-scale', type=float,
                        help='scale diameter by this factor if it is less than `dia-lim`',
                        default=None,
                        dest='diascale')
    parser.add_argument('--dia-lim', type=float,
                        help='minimum diameter over which it should be scaled up (um)',
                        default=None,
                        dest='dialim')
    
    args = parser.parse_args()
    inputs = args.inputs
    celltemplate = args.celltemplate
    syncount = args.syncount
    onset = args.onset
    tau = args.tau
    gmax = args.gmax
    diascale = args.diascale
    dialim = args.dialim
    mechs = {'pas': {'g': 3e-5, 'e': -51.0}}
             # {'name': 'caiiag', 'gbar': 2.5e-5},
             # {'name': 'ka', 'gbar': 0.03}]
             
    h.xopen(celltemplate)
    ggn = nu.create_cell(args.cellname, filename=args.celltemplate, mechparams=mechs)
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
                    # print('Changed diameter of segment', ii, seg.x, 'of section', sec.name(), 'to', seg.diam)
                    count += 1
            h.pop_section()
        print('Scaled diameters of', count, 'sections whose diameters were <', dialim, 'by', diascale)
    g0 = nu.nrngraph(ggn)
    g, nmap = ng.renumber_nodes(g0, ggn.soma.name())
    type_node_map = defaultdict(list)
    # Reverse lookup table for node by stype
    for n, d in g.nodes_iter(data=True):
        stype = d['s']
        if g.node[n]['orig'] is not None:   # Skip the dummy nodes at terminal
            type_node_map[stype].append(n)
    synnodes = []
    if args.synfile is None:
        for branch_id in inputs:
            nodes = type_node_map[branch_id]
            size = syncount
            if len(nodes) < syncount:
                size = len(nodes)
            for n in np.random.choice(nodes, size=size, replace=False):
                 # print(n, g.node[n])
                 synnodes.append(n)
    else:
        sec_node_map = get_section_node_map(g)
        with open(args.synfile) as fd:
            for line in fd:
                try:
                    secname = line.strip()
                    synnodes.append(sec_node_map[secname])
                except KeyError:
                    print('Could not fidn section "{}"'.format(secname))
                    raise
        
    # Select one random node from each stype
    recnodes = []
    if args.recfile is not None:
        sec_node_map = get_section_node_map(g)
        with open(args['recfile']) as fd:
            for line in fd:
                try:
                    secname = line.strip()
                    recnodes.append(sec_node_map[secname])
                except KeyError:
                    print('Could not fidn section "{}"'.format(secname))
                    raise        
    else:
        for stype, nodes in type_node_map.items():
            size = args.reccount
            if args.reccount > len(nodes):
                size = len(nodes)
            recnodes += list(np.random.choice(nodes, size=size, replace=False))
            # print(recnodes)
    recnodes = list(set(recnodes + synnodes))
    recsecs = [g.node[n]['orig'] for n in recnodes]
    synapses = []
    for snode in synnodes:
        synsec = g.node[snode]['orig']
        alpha_syn = insert_alphasynapse(synsec(1.0), onset=onset, gmax=gmax, tau=tau)
        synapses.append(alpha_syn)
    t_vec = h.Vector()
    t_vec.record(h._ref_t)    
    tabs = setup_recording(recsecs)
    h.tstop = onset + 100 * tau
    ts = datetime.now()
    h.finitialize(mechs['pas']['e'])
    h.fcurrent()
    while h.t < h.tstop:
        h.fadvance()
    
    # h.run()
    te = datetime.now()
    delta = te - ts
    print('Time for', h.tstop*1e-3, 's simulation =', delta.days * 86400 + delta.seconds + 1e-6 * delta.microseconds)
    outfilename = 'data/A_Vm_multi_syn_{}.h5'.format(ts.strftime('%Y%m%d_%H%M%S'))
    with h5.File(outfilename) as fd:
        sections = [g.node[s]['orig'].name() for s in synnodes]
        syninfo = fd.create_dataset('syninfo', data=sections)
        syninfo.attrs['tau'] = tau
        syninfo.attrs['onset'] = onset
        syninfo.attrs['gmax'] = gmax
        
        fd.attrs['g_pas'] = mechs['pas']['g']
        fd.attrs['e_pas'] = mechs['pas']['e']
        fd.attrs['RA'] = ggn.soma.Ra
        time = fd.create_dataset('time', data=np.asarray(t_vec))
        for v, sec, node in zip(tabs, recsecs, recnodes):
            ds = fd.create_dataset('v_{}'.format(node), data=np.array(v))
            ds.attrs['section'] = sec.name()
        model = fd.create_group('model')
        sec_names  = []
        sec_RA = []
        sec_g = []
        sec_len = []
        sec_dia = []
        sec_parent = []
        for sec in ggn.allsec():
            sec_names.append(sec.name())
            sec_RA.append(sec.Ra)
            sec_g.append(sec(0.5).pas.g)
            sec_len.append(sec.L)
            sec_dia.append(sec.diam)
            ref = h.SectionRef(sec=sec)
            if ref.has_parent():
                sec_parent.append( ref.parent.name())
            else:
                sec_parent.append('')

        name_ds = model.create_dataset('section', data=sec_names)
        ra_ds = model.create_dataset('RA', data=sec_RA)
        gpas_ds = model.create_dataset('gpas', data=sec_g)
        len_ds = model.create_dataset('length', data=sec_len)
        dia_ds = model.create_dataset('diameter', data=sec_dia)
        parent_ds =  model.create_dataset('parent', data=sec_parent)
    print('Data saved in', outfilename)
    plot_Vm(outfilename)
        
        

    


# 
# localized_input_output.py ends here
