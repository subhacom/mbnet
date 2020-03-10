# -*- coding: utf-8 -*-
"""
Display a NEURON celltemplate in 3D
"""


import sys
import os
import argparse
from neuron import h
import nrnutils as nu
import neurograph as ng
import networkx as nx

sys.path.append('D:/subhasis_ggn/model/morphutils')
sys.path.append('c:/nrn/lib/python')
os.environ['NEURONHOME'] = 'c:\nrn'

import morph3d


def cell_to_morph3d(cell):
    """Convert a neuron cell to a morphological tree.
    
    Here we create nodes by recovering the 3d point coordinates in the
    section.

    Returns: g, root where g is a digraph rooted at soma, root is the
    node label for the root node ({soma.name()}_0).

    This will need a cleanup for nrn7.5.

    """
    g = nx.DiGraph()
    stack = [cell.soma]
    while len(stack) > 0:
        sec = stack.pop()
        # This is roundabout way is required because nrn7.4 does not
        # provide explicit equivalent of `access {section}`. In nrn7.5
        # the 3d functions are available as Section methods.
        h('access {}'.format(sec.name()))
        stype = nu.sectype(sec.name())
        pt3d = int(h.n3d())
        # pt3d = int(sec.n3d()):  # only nrn >= 7.5
        for ii in range(pt3d):        
            name = '{}_{}'.format(sec.name(), ii)
            x = h.x3d(ii)
            y = h.y3d(ii)
            z = h.z3d(ii)
            d = h.diam3d(ii)
            g.add_node(name, x=x, y=y, z=z, r=d/2.0, s=stype, orig=sec)
        for ii in range(1, pt3d):
            n1 = '{}_{}'.format(sec.name(), ii-1)
            n2 = '{}_{}'.format(sec.name(), ii)
            length = ng.eucd(g, n1, n2)
            g.add_edge(n1, n2, length=length)
        current = h.SectionRef(sec=sec)
        if current.has_parent():
            h('access {}'.format(current.parent.name()))            
            n1 = '{}_{}'.format(current.parent.name(),
                                int(h.n3d()-1))
            g.add_edge(n1, '{}_0'.format(sec.name()),
                       length=0)
        for child in current.child:
            # print('Adding', child.name())
            stack.append(child)
    return g, '{}_0'.format(cell.soma.name())


if __name__ == '__main__':
    parser = argparse.ArgumentParser('Display a neuron cell template in 3D')
    parser.add_argument('-f', '--filename', type=str, required=True,
                        dest='filename',
                        help='file containing cell template definition')
    parser.add_argument('-c', '--cellname', type=str, required=True,
                        dest='cellname',
                        help='name of the cell template in the template file')
    parser.add_argument('-s', '--stick', action='store_true', dest='stick',
                        help='Ball and stick model. If false, the 3D locations using pt3d are used')
    parser.add_argument('-o', '--origlabel', action='store_true', dest='origlabel',
                        help='If the original neuron section names should be used.')
    args = parser.parse_args()
    h.xopen(args.filename)
    cellclass = getattr(h, args.cellname)
    cell = cellclass()
    if args.stick:
        graph = nu.nrngraph(cell)
        for n in graph.nodes:
            print(n)
        root = '{}.soma'.format(list(graph.nodes)[0].partition('.')[0])
    else:
        graph, root = cell_to_morph3d(cell)
    numeric, nmap = ng.renumber_nodes(graph, start=root)
    print(len(graph.nodes))
    print(max(numeric.nodes))
    print(len(numeric.nodes))
    if args.origlabel:
        reverse = {node: name.split('.')[-1] for name, node in nmap.items() }
        label_nodes = reverse.keys()
        labels = reverse.values()
        
    else:
        label_nodes = list(numeric.nodes)
        labels = label_nodes
        
    morph3d.show_neurons({args.cellname: numeric}, label_nodes_dict={args.cellname: label_nodes},
                         labels_dict={args.cellname: labels}, lines=args.stick)
