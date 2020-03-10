"""Check the effect of synchronized synaptic input at random locations
of specific branches of the GGN.

"""
from __future__ import print_function
import sys
import argparse
from datetime import datetime
import h5py as h5
import nrnutils as nu
from nrnutils import *
import ephys
from localized_input_output import plot_Vm

def update_section(sec, args):
    sec.Ra = args.get('RA', sec.Ra)
    if 'gpas' in args:
        ephys.set_mech_param(sec, 'pas', 'g', args['gpas'])
        ephys.set_mech_param(sec, 'pas', 'e', args['Em'])
    if ('gk' in args) and (args['gk'] > 0):
        ephys.set_mech_param(sec, 'ka', 'gbar', args['gk'])
    if ('gca' in args) and (args['gca'] > 0):
        ephys.set_mech_param(sec, 'cam', 'gbar', args['gca'])


def setup_model(args):
    """Setup the model with synapses and recording tables."""
    cell_template = args['celltemplate']
    input_regions = args['inputs']
    syn_counts = args['syncounts']
    onset = args['onset']
    tau = args['tau']
    gsyn = args['gsyn']
    mechs = {'pas': {'g': Q_(args['gpas'], 'S/cm**2'), 'e': Q_(args['Em'], 'mV')}}
    if args['gk'] > 0:
        mechs['ka'] = {'gbar': Q_(args['gk'], 'S/cm**2')}
    if args['gca'] > 0:
        mechs['cam'] = {'gbar': Q_(args['gca'], 'S/cm**2')}             
    h.xopen(cell_template)
    ggn = nu.create_cell(args['cellname'], filename=args['celltemplate'], mechparams=mechs)
    for sec in ggn.allsec():
        sec.Ra = args['RA']
    g0 = nu.nrngraph(ggn)        
    g, nmap = ng.renumber_nodes(g0, ggn.soma.name())
    # Pick out nodes with sections associated with them (not terminal dummy nodes)
    sec_nodes = [n for n in g.nodes() if g.node[n]['orig'] is not None]
    # Select the nodes for synapse insertion
    if args['synfile'] is not None:  # This option is to allow replication
        syn_nodes = []
        sec_node_map = get_section_node_map(g)
        with open(args['synfile']) as fd:
            for line in fd:
                try:
                    secname = line.strip()
                    syn_nodes.append(sec_node_map[secname])
                except KeyError:
                    print('Could not fidn section "{}"'.format(secname))
                    raise
    else:
        syn_nodes_by_sid = ng.select_random_nodes_by_sid(g.subgraph(sec_nodes), args['inputs'], args['syncounts'])
        syn_nodes = list(np.concatenate([v for v in syn_nodes_by_sid.values() if len(v) > 0]))
        # If np.concatenate arg contains empty list, the result gets converted to float
    syn_secs = []
    synapses = []
    for syn_node in syn_nodes:
        sec = g.node[syn_node]['orig']
        syn_secs.append(sec)
        syn = insert_alphasynapse(sec(1.0), onset=args['onset'], gmax=args['gsyn'], tau=args['tau'])
        synapses.append(syn)
    if args['recfile'] is not None:
        rec_nodes = []
        sec_node_map = get_section_node_map(g)
        with open(args['recfile']) as fd:
            for line in fd:
                try:
                    secname = line.strip()
                    rec_nodes.append(sec_node_map[secname])
                except KeyError:
                    print('Could not fidn section "{}"'.format(secname))
                    raise        
    else:
        if args['recbyreg']:
            rec_regions =  [1, 3, 4, 5, 6, 7, 8]
            rec_nodes_by_sid = ng.select_random_nodes_by_sid(g.subgraph(sec_nodes), rec_regions, [args['reccount']] * len(rec_regions))
            rec_nodes = np.concatenate([v for v in rec_nodes_by_sid.values() if len(v) > 0])
            rec_nodes = list(rec_nodes)
        else:
            rec_nodes = list(np.random.choice(sec_nodes, size=args['reccount'], replace=False))
    rec_nodes = list(set(rec_nodes + syn_nodes))
    rec_secs = [g.node[n]['orig'] for n in rec_nodes]
    t_vec = h.Vector()
    t_vec.record(h._ref_t)    
    v_vecs = setup_recording(rec_secs)
    syn_vecs = []
    for syn in synapses:
        vec = h.Vector()
        vec.record(syn._ref_i)
        syn_vecs.append(vec)
    return {'cell': ggn,
            'cellgraph': g,
            'synapses': synapses,
            'synsecs': syn_secs,
            'synnodes': syn_nodes,
            'recsecs': rec_secs,
            'recnodes': rec_nodes,
            'tvec': t_vec,
            'vvecs': v_vecs,
            'synvecs': syn_vecs}


def init_and_run(args, model_dict=None):
    """Initialize the model and run simulation using parameters from
    args."""
    if model_dict is None:
        model_dict = setup_model(args)
    else:
        for sec in model_dict['cell'].allsec():
            update_section(sec, args)
        model_dict['cell'].geom_nseg() 
    simtime = args['simtime']
    if simtime < args['onset']:
        simtime = args['onset'] + 100.0  # Arbitrarily make it at least up to 100 ms after input
        args['simtime'] = simtime
    h.tstop = simtime
    print('Starting simulation with argumets:\n')
    print(args)    
    start = datetime.now()
    if ('Erest' in args) and (args['Erest'] is not None):
        h.finitialize(args['Erest'])
    else:
        h.finitialize(args['Em'])
    h.fcurrent()
    while h.t < h.tstop:
        h.fadvance()
    end = datetime.now()
    dt = end - start
    print('Finished simulation # {} of {} ms in {} s'.format(args['pid'], h.tstop,
                                                             dt.days * 86400 + dt.seconds + dt.microseconds * 1e-6))
    return model_dict
        

def make_parser():
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
    parser.add_argument('-g', '--gsyn', type=float,
                        help='peak conductance of the synapses (in uS)',
                        dest='gsyn',
                        default=1e-3)
    parser.add_argument('-t', '--tau', type=float,
                        help='time constant of the alpha synapse (in ms)',
                        dest='tau', default=1.0)
    parser.add_argument('-s', '--onset', type=float,
                        help='onset time of synaptic input (in ms)',
                        dest='onset', default=50.0)
    parser.add_argument('-n', '--syncount', type=int,
                        help='number of synapses',
                        dest='syncounts', action='append', required=True)
    parser.add_argument('--synfile', type=str,
                        help='File containing list of synapse sections',
                        dest='synfile', default=None)
    parser.add_argument('--reccount', type=int,
                        help='number of compartments to record from',
                        dest='reccount', default=10)
    parser.add_argument('--recfile', type=str,
                        help='File containing list of recording sections',
                        dest='recfile', default=None)
    parser.add_argument('-c', '--cell', type=str,
                        help='cell name',
                        dest='cellname',
                        default='GGN_B99')
    parser.add_argument('-f', '--cellfile', type=str,
                        help='cell template file',
                        dest='celltemplate',
                        default='GGN_B99_20161205_d_lambda.hoc')
    parser.add_argument('--RMmin', type=float,
                        help='Specific membrane resistance (ohm-cm2) minimum',
                        dest='RMmin', default=33e3)  # 33 Kohm-cm2 according to Laurent, 1990
    parser.add_argument('--RMmax', type=float,
                        help='Specific membrane resistance (ohm-cm2) maximum',
                        dest='RMmax', default=33e3)
    parser.add_argument('--RMdivs', type=int,
                        help='Number of divisions in RM parameter search',
                        dest='RMdivs', default=1)
    parser.add_argument('--RAmin', type=float,
                        help='Specific axial resistance (ohm-cm) minimum',
                        dest='RAmin', default=100.0)
    parser.add_argument('--RAmax', type=float,
                        help='Specific axial resistance (ohm-cm) maximum',
                        dest='RAmax', default=100.0)
    parser.add_argument('--RAdivs', type=int,
                        help='Number of divisions in RA parameter search',
                        dest='RAdivs', default=1)
    parser.add_argument('-e', '--Em', type=float,
                        help='Leak reversal potential (mV)',
                        dest='Em', default=-51.0)
    parser.add_argument('--gca', type=float, help='conductance density of Ca channel in S/cm2',
                        dest='gca', default=6.1e-3)
    parser.add_argument('--gk', type=float, help='conductance density of K channel in S/cm2',
                        dest='gk', default=2.1e-3)
    parser.add_argument('-p', '--pid', type=int,
                        help='Process ID used for keeping file name unique',
                        dest='pid', default=np.random.randint(1000))
    parser.add_argument('--simtime', type=float,
                        help='simulation time (ms)',
                        dest='simtime', default=0)
    parser.add_argument('-d', '--datadir', type=str, 
                        help='Data directory. The files will saved in its subdirectories by date of simulation.',
                        dest='datadir', default='/data/rays3/ggn/single_syn_input/')
    parser.add_argument('--rec-by-region', action='store_true',
                        help='record `reccount` sections in every region if possible',
                        dest='recbyreg')
    
    # parser.add_argument('--dia-scale', type=float,
    #                     help='scale diameter by this factor if it is less than `dia-lim`',
    #                     default=None,
    #                     dest='diascale')
    # parser.add_argument('--dia-lim', type=float,
    #                     help='minimum diameter over which it should be scaled up (um)',
    #                     default=None,
    #                     dest='dialim')
    return parser


def save_data(model_dict, arg_dict, outfilename):
    with h5.File(outfilename) as fd:
        fd.attrs['dt'] = h.dt
        fd.attrs['tunit'] = 'ms'
        simgroup = fd.create_group('simulation_{}'.format(arg_dict['pid']))
        simgroup.attrs['timestamp'] = arg_dict['ts'].strftime('%Y%m%d_%H%M%S')
        sections = [sec.name() for sec in model_dict['synsecs']]
        syninfo = simgroup.create_dataset('syninfo', data=sections)
        syninfo.attrs['tau'] = arg_dict['tau']
        syninfo.attrs['onset'] = arg_dict['onset']
        syninfo.attrs['gmax'] = arg_dict['gsyn']
        simgroup.attrs['g_pas'] = arg_dict['gpas']
        simgroup.attrs['RM'] = arg_dict['RM']
        simgroup.attrs['RA'] = arg_dict['RA']
        simgroup.attrs['e_pas'] = arg_dict['Em']
        time = simgroup.create_dataset('time', data=np.asarray(model_dict['tvec']))
        for v, sec, node in zip(model_dict['vvecs'], model_dict['recsecs'], model_dict['recnodes']):
            ds = simgroup.create_dataset('v_{}'.format(node), data=np.array(v))
            ds.attrs['section'] = sec.name()
            ds.attrs['dt'] = h.dt
            ds.attrs['tunit'] = 'ms'
        for g, sec, node in zip(model_dict['synvecs'], model_dict['synsecs'], model_dict['synnodes']):
            ds = simgroup.create_dataset('isyn_{}'.format(node), data=np.array(g))
            ds.attrs['section'] = sec.name()
            ds.attrs['dt'] = h.dt
            ds.attrs['tunit'] = 'ms'
        model = simgroup.create_group('model')
        sec_names  = []
        sec_RA = []
        sec_g = []
        sec_len = []
        sec_dia = []
        sec_parent = []
        for sec in model_dict['cell'].allsec():
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
    

if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()
    if len(args.syncounts) != len(args.inputs):
        print('A synapse count must be specified for each input region')
        sys.exit(0)
    arg_dict = vars(args)
    RA = np.linspace(args.RAmin, args.RAmax, args.RAdivs)
    RM = np.linspace(args.RMmin, args.RMmax, args.RMdivs)
    model_dict = None
    simid = 0
    ts = datetime.now()
    outfilename = os.path.join(args.datadir, 'B_Vm_multi_syn_{}.h5'.format(ts.strftime('%Y%m%d_%H%M%S')))
    for ra in RA:        
        ts = datetime.now()
        arg_dict['ts'] = ts
        arg_dict['RA'] = ra        
        for rm in RM:
            arg_dict['pid'] = simid
            arg_dict['gpas'] = 1.0 / rm
            arg_dict['RM'] = rm
            model_dict = init_and_run(arg_dict, model_dict)
            save_data(model_dict, arg_dict, outfilename)
            simid += 1
    # plot_Vm(outfilename)
        
        

    


# 
# localized_input_output.py ends here
