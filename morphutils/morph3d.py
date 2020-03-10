# morph3d_vtk.py ---
#
# Filename: morph3d_vtk.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Thu Jul 14 18:34:15 2016 (-0400)
# Version:
# Package-Requires: ()
# Last-Updated: Thu Sep  6 16:58:04 2018 (-0400)
#           By: Subhasis Ray
#     Update #: 1389
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

"""Reference example from:
http://www.vtk.org/Wiki/VTK/Examples/Python/GeometricObjects/Display/ColoredLines

Modified with ref to:
http://www.vtk.org/Wiki/VTK/Examples/Cxx/VisualizationAlgorithms/TubesWithVaryingRadiusAndColors

"""

from __future__ import print_function
from collections import defaultdict
import sys
import numpy as np
import argparse
import ast
import networkx as nx
import neurograph as ng
import vtk
from vtk.util import numpy_support as vtknp


fndict = {'x': vtk.vtkTransform.RotateX,
          'y': vtk.vtkTransform.RotateY,
          'z': vtk.vtkTransform.RotateZ}


def create_labels(graph, label_nodes, labels=[], priorities=[],
                  transform=None):
    label_points = vtk.vtkPoints()
    label_str = vtk.vtkStringArray()
    label_str.SetName('Labels')
    for ii, n in enumerate(label_nodes):
        label_points.InsertNextPoint((graph.node[n]['x'],
                                      graph.node[n]['y'],
                                      graph.node[n]['z']))
        if ii < len(labels):
            label = str(labels[ii])
        else:
            label = str(n)
        label_str.InsertNextValue(label)
    label_polydata = vtk.vtkPolyData()
    label_polydata.SetPoints(label_points)
    label_polydata.GetPointData().AddArray(label_str)
    hierarchy = vtk.vtkPointSetToLabelHierarchy()
    hierarchy.SetLabelArrayName(label_str.GetName())
    if (len(priorities) > 0) and (len(priorities) == len(label_nodes)):
        parray = vtk.vtkIntArray()
        parray.SetName('Priorities')
        for priority in priorities:
            parray.InsertNextValue(priority)
        hierarchy.SetPriorityArrayName(parray.GetName())
    label_mapper = vtk.vtkLabelPlacementMapper()
    label_mapper.SetInputConnection(hierarchy.GetOutputPort())
    if transform is not None:
        transfilter = vtk.vtkTransformFilter()
        transfilter.SetInputData(label_polydata)
        transfilter.SetTransform(transform)
        transfilter.Update()
        hierarchy.SetInputConnection(transfilter.GetOutputPort())
    else:
        hierarchy.SetInputData(label_polydata)
    label_actor = vtk.vtkActor2D()
    label_actor.SetMapper(label_mapper)
    return label_actor


def nrngraph2polydata(ngraph, nodecolor=ng.nodecolor_4cp, lines=False):
    """Convert a graph representation of a neuron into a VTK PolyData
    object."""
    node_map = {n: ii for ii, n in enumerate(ngraph.nodes())}
    nodes = vtk.vtkPoints()
    nodes.SetNumberOfPoints(len(ngraph))
    for ii, n in enumerate(ngraph.nodes()):
        nodes.SetPoint(ii, (ngraph.node[n]['x'],
                            ngraph.node[n]['y'],
                            ngraph.node[n]['z']))
    # Create the representation of the edges
    edges = vtk.vtkCellArray()
    colors = vtk.vtkUnsignedCharArray()
    colors.SetName('Colors')
    colors.SetNumberOfComponents(3)
    colors.SetNumberOfTuples(len(ngraph.edges()))
    for ii, (n0, n1) in enumerate(ngraph.edges()):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, node_map[n0])
        line.GetPointIds().SetId(1, node_map[n1])
        edges.InsertNextCell(line)
        colors.SetTuple3(ii, *nodecolor[ngraph.node[n0]['s'] % len(nodecolor)])
    polydat = vtk.vtkPolyData()
    polydat.SetPoints(nodes)
    polydat.SetLines(edges)
    polydat.GetCellData().SetScalars(colors)
    if not lines:  # Render tree with cylinders
        radius = np.array([ngraph.node[n]['r'] for n in ngraph])
        radius = vtknp.numpy_to_vtk(radius, deep=True,
                                    array_type=vtk.VTK_FLOAT)
        radius.SetName('Radius')
        polydat.GetPointData().AddArray(radius)
        polydat.GetPointData().SetActiveScalars(radius.GetName())
        polydat.GetCellData().AddArray(colors)
    return polydat


def get_cylinder_view(polydat):
    tube_filt = vtk.vtkTubeFilter()
    tube_filt.SetNumberOfSides(10)
    tube_filt.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
    tube_filt.SetInputData(polydat)
    return tube_filt


def get_transformed(inputobj, soma, pos, rot, relative):
    """pos: position of the object (if relative, that of soma).

    rot: list of tuples [('x', xdegrees), ('y', ydegrees), ('z', zdegrees)]"""
    print('ROTATION', rot)
    transform = vtk.vtkTransform()
    transform.Identity()
    transform.PostMultiply()   # this sets order of transformations
    transform.Translate(-soma['x'], -soma['y'], -soma['z'])
    if rot is not None:
        for axis, angle in rot:
            fndict[axis](transform, angle)
    if not relative:
        transform.Translate(soma['x'], soma['y'], soma['z'])
    if pos is not None:
        transform.Translate(*pos)
    transfilter = vtk.vtkTransformFilter()
    transfilter.SetTransform(transform)
    if isinstance(inputobj, vtk.vtkPolyData):
        transfilter.SetInputData(inputobj)
    else:
        transfilter.SetInputConnection(inputobj.GetOutputPort())
    return transfilter, transform


def neuron3d(neurotree, label_nodes=[], labels=[], priorities=[],
             pos=[0, 0, 0], rot=[], nodecolor=ng.nodecolor_4cp,
             lines=False, relative=True):
    """Create 3D display for neurons up to the level of vtk actors for
    neurons in neurotree, with also actors for labels if specified.

    neurotree: graph specifying the morphology.

    label_nodes: list of nodes in neurotree that should be labeled

    labels: list of strings specifying actual label text

    priorities: display priorities of the labels

    pos: position of the object (if relative, that of soma).

    rot: list of tuples [('x', xdegrees), ('y', ydegrees), ('z', zdegrees)]

    lines: if False, show neuron as tubes (with diameters specified in
    morphology), if True, show it as lines of width 1, if a positive
    number, lines of that width.

    relative: if the positions of the neuron nodes are to be moved
    such that soma is origin (as opposed to absolute).

    """
    mapper = vtk.vtkPolyDataMapper()
    soma = neurotree.node[1]
    polydata = nrngraph2polydata(neurotree, nodecolor=nodecolor,
                                 lines=lines)
    if lines <= 0.0:
        ggn = get_cylinder_view(polydata)
    else:
        ggn = polydata
    # Rotation is always with soma as origin
    transfilter, transform = get_transformed(ggn, soma, pos, rot, relative)
    mapper.SetInputConnection(transfilter.GetOutputPort())
    mapper.ScalarVisibilityOn()
    mapper.SetScalarModeToUseCellFieldData()
    mapper.SelectColorArray('Colors')
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    if lines > 0:
        actor.GetProperty().SetLineWidth(lines)
        print('Set linewidth to', lines)
    if len(label_nodes) > 0:
        label_actor = create_labels(neurotree, label_nodes,
                                    labels=labels,
                                    priorities=priorities,
                                    transform=transform)
    else:
        label_actor = None
    return actor, label_actor


def show_neurons(tree_dict, label_nodes_dict={}, labels_dict={},
                 priorities_dict={}, pos_dict={}, rot_dict={},
                 nodecolor=ng.nodecolor_4cp, background=(0, 0, 0),
                 lines=False, relative=True, stereo=False, axes=True,
                 fullscreen=True):
    """
    tree_dict: dict containing {cell_name: morphology tree}
    label_nodes_dict: {cell_name: list of label nodes}
    labels_dict: {cell_name: list of labels corresponding to label nodes}
    priorities_dict: {cell_name: list of priorities of the labels}
    axes: If true, show axis with scale bar info
    """
    label_actors = {}
    neuron_actors = {}
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(*background)

    for cell, neurotree in tree_dict.items():
        label_nodes = label_nodes_dict.get(cell, [])
        labels = labels_dict.get(cell, [])
        priorities = priorities_dict.get(cell, [])
        pos = pos_dict.get(cell, [0, 0, 0])
        rot = rot_dict.get(cell, [])
        neuron_actor, label_actor = neuron3d(neurotree,
                                             label_nodes=label_nodes,
                                             labels=labels,
                                             priorities=priorities,
                                             pos=pos, rot=rot,
                                             nodecolor=nodecolor,
                                             lines=lines,
                                             relative=relative)
        neuron_actors[cell] = neuron_actor
        label_actors[cell] = label_actor
        renderer.AddActor(neuron_actor)
        if label_actor is not None:
            renderer.AddActor2D(label_actor)
    if axes:  # show axis with scale bar info
        _ax = vtk.vtkCubeAxesActor2D()
        _ax.SetLabelFormat('%3.0f')
        _ax.SetNumberOfLabels(0)
        # _ax.SetYAxisVisibility(False)
        # _ax.SetZAxisVisibility(False)
        _ax.SetBounds(0, 200, 0, 200, 0, 200)
        _ax.SetXLabel('X')
        _ax.SetYLabel('Y')
        _ax.SetZLabel('Z')
        _ax.SetXOrigin(0)
        _ax.SetYOrigin(0)
        _ax.SetZOrigin(0)
        color = (1.0-background[0], 1.0-background[1], 1.0-background[2])
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(color)
        _ax.SetAxisLabelTextProperty(tprop)
        _ax.GetProperty().SetColor(*color)
        # _ax.SetFlyModeToClosestTriad()
        _ax.SetCamera(renderer.GetActiveCamera())
        renderer.AddActor(_ax)
    renderer.ResetCamera()
    win = vtk.vtkRenderWindow()
    win.AddRenderer(renderer)
    if fullscreen:
        win.FullScreenOn()
    if stereo:
        win.GetStereoCapableWindow()
        win.StereoCapableWindowOn()
        win.SetStereoRender(1)
        win.SetStereoTypeToCrystalEyes()
        win.SetFullScreen(1)
        #win.SetStereoTypeToRedBlue()

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(win)

    win.Render()
    interactor.Initialize()
    interactor.Start()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Display morphology from' +
                                     ' an swc file')
    parser.add_argument('-f', '--filename', type=str, dest='filenames',
                        action='append',
                        help='SWC morphology file to display.' +
                        ' It can be specified multiple times.')
    parser.add_argument('-t', '--translate', type=str,
                        dest='translate', default=[], action='append',
                        help='shift in X, Y and Z position of soma.' +
                        ' Should be a list inside quotes,' +
                        ' e.g., "(100, 50, 10)". It can be specified' +
                        ' multiple times, each applied to ' +
                        'corresponding file.')
    parser.add_argument('-r', '--rotate', type=str, dest='rotate',
                        default=[], action='append',
                        help='Rotations in degrees around x, y and z axes.' +
                        ' Should be a axis and degrees inside quotes in the ' +
                        'order to be applied,' +
                        'e.g., "x 100 y 50 z 10". Can be specified multiple' +
                        ' times, each applied to the corresponding input' +
                        ' morphology file with soma as the pivot point.')
    parser.add_argument('-m', '--mirror', type=str, dest='mirror',
                        default=[], action='append',
                        help='Mirror along the axes (on the plane normal to' +
                        ' the axis through origin, i.e. x for reflecting on ' +
                        'YZ plane. Applied before rotation or translation.' +
                        ' Should be one or more axis names  in the order to' +
                        ' be applied')
    parser.add_argument('-l', '--lines', type=float, dest='lines',
                        default=0.0, help='If > 0, use wire diagram of ' +
                        'specified linewidth else tubes for morphology ' +
                        'display')
    parser.add_argument('-c', '--colormap', type=str,
                        default='3cd2', dest='colormap',
                        help='colormap for displaying the different ' +
                        'structures, 3cd2, 3cp, 3cs, 4cp, 5cd2, 5cs3,' +
                        '10cp, 10cs3 are the available values. ' +
                        'They have 3 to 10 different colors for ' +
                        'cells with as many different SWC structure ids')
    parser.add_argument('-b', '--black-background', type=str,
                        default='0 0 0',  dest='bg',
                        help='Set background color. Should be a string' +
                        ' specifying "r g b" where each entry is in [0, 1].')
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
    parser.add_argument('--relative', action='store_true',
                        help='Make soma origin for each cell',
                        dest='relative')
    args = parser.parse_args()
    try:
        colormap = ng.colormaps[args.colormap]
    except KeyError:
        print('Colormap must be one of', ng.colormaps.keys())
        sys.exit(0)
    background = [0, 0, 0]
    for ii, fstr in enumerate(args.bg.split()):
        background[ii] = float(fstr)
        if ii >= 3:
            break
    cell_graphs = {}
    cell_pos = {}
    cell_rot = {}
    label_nodes = defaultdict(list)
    for ii, fname in enumerate(args.filenames):
        cg = ng.tograph(fname)
        cell_graphs[fname] = cg
        if ii < len(args.mirror):
            maxes = args.mirror[ii].split()
            for ax in maxes:
                for n in cg.nodes():
                    cg.node[n][ax] = -cg.node[n][ax]
        if ii < len(args.translate):
            cell_pos[fname] = [int(x) for x in args.translate[ii].split()]
        else:
            cell_pos[fname] = [0, 0, 0]
        if ii < len(args.rotate):
            rot = args.rotate[ii].split()
            for jj in range(0, len(rot), 2):
                rot[jj+1] = float(rot[jj+1])
            cell_rot[fname] = zip(rot[0::2], rot[1::2])
        else:
            cell_rot[fname] = []
        if args.branches:
            lnodes, edges = zip(*ng.branch_points(cg))
            label_nodes[fname] = list(lnodes)
        if args.leaves:
            label_nodes[fname] += [node for node in cg.nodes()
                                   if cg.degree(node) == 1]
        if len(args.sid) > 0:
            sid_node_map = ng.get_stype_node_map(cg)
            for sid in args.sid:
                label_nodes[fname] += sid_node_map[sid]
    show_neurons(cell_graphs, lines=args.lines, label_nodes_dict=label_nodes,
                 labels_dict={}, pos_dict=cell_pos, rot_dict=cell_rot,
                 nodecolor=colormap,
                 background=background, axes=args.scalebar,
                 fullscreen=args.fullscreen,
                 relative=args.relative)




#
# morph3d_vtk.py ends here
