# displaycell.py ---
#
# Filename: displaycell.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Wed Jul 13 17:04:41 2016 (-0400)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Mar 10 15:14:09 2020 (-0400)
#           By: Subhasis Ray
#     Update #: 562
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

# Commentary:
#
#
#
#

# Code:
from __future__ import print_function
import sys
import numpy as np
import argparse
import ast
import networkx as nx
import neurograph as ng
import morphoplot as mp
from matplotlib import pyplot as plt
# plt.style.use('ggplot')


def parse_num_list(arg):
    ret = ast.literal_eval(arg)
    return ret


def make_tmatrix(xshift, yshift, zshift,
                 thetax, thetay, thetaz):
    """Create a rigid transformation matrix given the parameters.

    theta are angles in degrees.

    The transformations are aplied in the sequence: xrotation,
    yrotation, z rotation, translation.

    """
    tmatrix = np.identity(4)
    tmatrix[:-1, -1] = xshift, yshift, zshift
    degtorad = np.pi / 180.0
    thetax *= degtorad
    thetay *= degtorad
    thetaz *= degtorad
    xcos = np.cos(thetax)
    ycos = np.cos(thetay)
    zcos = np.cos(thetaz)
    xsin = np.sin(thetax)
    ysin = np.sin(thetay)
    zsin = np.sin(thetaz)
    xrot = np.identity(4)
    yrot = np.identity(4)
    zrot = np.identity(4)
    xrot[1, 1] = xcos
    xrot[2, 2] = xcos
    xrot[1, 2] = - xsin
    xrot[2, 1] = xsin
    yrot[0, 0] = ycos
    yrot[0, 2] = ysin
    yrot[2, 0] = - ysin
    yrot[2, 2] = ycos
    zrot[0, 0] = zcos
    zrot[0, 1] = - zsin
    zrot[1, 0] = zsin
    zrot[1, 1] = zcos
    print('Xrot')
    print(xrot)
    print('Yrot')
    print(yrot)
    print('Zrot')
    print(zrot)
    transform = tmatrix.dot(zrot.T).dot(yrot.T).dot(xrot.T)
    return transform


def rigid_transform_graph(cellgraph, transform):
    """Update graph by applying rigid transform on every node of
    cellgraph. The original cellgraph is modified."""
    print('Rigid transformation matrix')
    print(transform)
    for node, attr in cellgraph.nodes(data=True):
        pos = np.array([attr['x'], attr['y'], attr['z'], 1.0])
        newpos = transform.dot(pos.T)
        attr['x'] = newpos[0]
        attr['y'] = newpos[1]
        attr['z'] = newpos[2]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Display morphology from' +
                                     ' an swc file')
    parser.add_argument('-f', '--filename', type=str, dest='filenames',
                        action='append',
                        help='SWC morphology file to display.' +
                        ' It can be specified multiple times.')
    parser.add_argument('-t', '--translate', type=parse_num_list,
                        dest='translate', default=[], action='append',
                        help='shift in X, Y and Z position of morphology.' +
                        ' Should be a list inside quotes,' +
                        ' e.g., "(100, 50, 10)". It can be specified' +
                        ' multiple times, each applied to ' +
                        'corresponding file.')
    parser.add_argument('-r', '--rotate', type=parse_num_list, dest='rotate',
                        default=[], action='append',
                        help='Rotations in degrees around x, y and z axes.' +
                        ' Should be a list inside quotes,' +
                        'e.g., "(100, 50, 10)". Can be specified multiple' +
                        ' times, each applied to the corresponding input' +
                        ' morphology file.')
    parser.add_argument('-m', '--module', type=str, dest='module',
                        default='vtk',
                        help='module for 3d visualization: ' +
                        'mplot3d, vpython, vispy or vtk, default is vtk.')
    parser.add_argument('-l', '--lines', default=0.0, dest='lines',
                        help='whether to use wire diagram' +
                        ' or tubes for morphology display.' +
                        ' If > 0, that width lines are used.')
    parser.add_argument('-c', '--colormap', type=str,
                        default='3cd2', dest='colormap',
                        help='colormap for displaying the different ' +
                        'structures, 3cd2, 3cp, 3cs, 4cp, 5cd2, 5cs3,' +
                        '10cp, 10cs3 are the available values. ' +
                        'They have 3 to 10 different colors for ' +
                        'cells with as many different SWC structure ids')
    parser.add_argument('-b', '--black-background', action='store_true',
                        help='Set background color to black', dest='bg')
    parser.add_argument('--branches', action='store_true',
                        help='Label branching nodes',
                        default=False, dest='branches')
    parser.add_argument('--leaves', action='store_true',
                        help='Label leaf nodes',
                        default=False, dest='leaves')
    parser.add_argument('-s', '--struct-id', type=int, nargs='+',
                        dest='sid', help='label nodes with specified' +
                        ' structure ids', default=[])
    parser.add_argument('-a', '--scalebar', action='store_true',
                        help='Show scale bar', dest='scalebar')
    parser.add_argument('-F', '--fullscreen', action='store_true',
                        help='Display fullscreen (only for VTK)',
                        dest='fullscreen')
    args = parser.parse_args()
    module = args.module
    lines = float(args.lines)
    print(args.rotate)
    try:
        colormap = ng.colormaps[args.colormap]
    except KeyError:
        print('Colormap must be one of', ng.colormaps.keys())
        sys.exit(0)
    if args.bg:
        background = (0, 0, 0)
    else:
        background = (1, 1, 1)
    if module == 'vispy':
        from morph3d_vispy import neuron3d
    elif module == 'vpython':
        from morph3d_vpython import neuron3d
    else:
        from morph3d_vtk import neuron3d
    # Build the graph
    combined_cellgraph = nx.DiGraph()
    for ii, fname in enumerate(args.filenames):
        print('Opening', fname)
        cellgraph = ng.tograph(fname)
        xshift, yshift, zshift = (0, 0, 0)
        print('##', args.translate)
        if len(args.translate) > ii:
            xshift, yshift, zshift = args.translate[ii]
        xrot, yrot, zrot = (0, 0, 0)
        if len(args.rotate) > ii:
            xrot, yrot, zrot = args.rotate[ii]
        tmatrix = make_tmatrix(xshift, yshift, zshift, xrot, yrot, zrot)
        print('Transformation matrix for', fname)
        print(tmatrix)
        rigid_transform_graph(cellgraph, tmatrix)
        if ii > 0:
            maxnode = max(combined_cellgraph.nodes())
            for node in cellgraph.nodes():
                cellgraph.node[node]['p'] += maxnode
            mapping = {n: n + maxnode for n in cellgraph.nodes()}
            cellgraph = nx.relabel_nodes(cellgraph, mapping)
        combined_cellgraph = nx.union(combined_cellgraph, cellgraph)
    label_nodes = []
    if args.branches:
        label_nodes, degrees = zip(*ng.branch_points(combined_cellgraph))
    label_nodes = list(label_nodes)
    if args.leaves:
        label_nodes += [node for node in combined_cellgraph.nodes()
                        if combined_cellgraph.degree(node) == 1]
    if len(args.sid) > 0:
        sid_node_map = ng.get_stype_node_map(combined_cellgraph)
        for sid in args.sid:
            label_nodes += sid_node_map[sid]
    if module == 'mplot3d':
        ax = mp.plot_3d_lines(combined_cellgraph)
        plt.show()
    elif module == 'vtk':
        neuron3d(combined_cellgraph, lines=lines, label_nodes=label_nodes,
                 labels=label_nodes, nodecolor=colormap,
                 background=background, axes=args.scalebar,
                 fullscreen=args.fullscreen)
    elif module == 'vispy':
        neuron3d(combined_cellgraph)
    elif module == 'vpython':
        neuron3d(combined_cellgraph)
    else:
        print('Unknown module', module)

#
# displaycell.py ends here
