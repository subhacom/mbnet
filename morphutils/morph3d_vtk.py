# morph3d_vtk.py --- 
# 
# Filename: morph3d_vtk.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Jul 14 18:34:15 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Tue Aug 15 18:16:02 2017 (-0400)
#           By: Subhasis Ray
#     Update #: 695
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

from collections import defaultdict
import vtk
import numpy as np
import networkx as nx
from neurograph import nodecolor_10cp, nodecolor_4cp

def create_labels(graph, label_nodes, labels=[], priorities=[]):
    label_points = vtk.vtkPoints()
    label_str = vtk.vtkStringArray()
    label_str.SetName('Labels')
    for ii, n in enumerate(label_nodes):
        label_points.InsertNextPoint((graph.node[n]['x'],
                                      graph.node[n]['y'],
                                      graph.node[n]['z']))
        if len(labels) == len(label_nodes):
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
    if vtk.VTK_MAJOR_VERSION <= 5:
        hierarchy.SetInput(label_polydata)
    else:
        hierarchy.SetInputData(label_polydata)
    label_actor = vtk.vtkActor2D()
    label_actor.SetMapper(label_mapper)
    return label_actor



def nrngraph2vtk(neuron_graph, label_nodes=[], labels=[], priorities=[], nodecolor=nodecolor_4cp, background=(0,0,0), lines=0, stereo=False):
    """Builds a VTK pipeline to render the neuron_graph in 3D.
    label_nodes: list of nodes to be labeled

    labels: list of labels for the nodes (must match label_nodes)

    priorities: display priorities of the label nodes. This determines
    which labels will be displayed when not enough space is available
    to show them all.

    nodecolor: dict specifying color for different node structure ids
    {sid: (r, g, b)}

    background: background colour specified as (r, g, b) 

    lines: whether to display the neuron with fixed width lines or as
    tapered cylinders.

    stereo: whether to use stereoscopic display mode.

    Returns (renderer, actor)

    """
    ncyl = len(neuron_graph.edges())
    nodes = vtk.vtkPoints()
    nodes.SetNumberOfPoints(len(neuron_graph))
    node_map = {}
    radius = vtk.vtkFloatArray()
    radius.SetName('Radius')
    # radius.SetNumberOfValues(len(neuron_graph))
    radius.SetNumberOfValues(len(neuron_graph))
    for ii, n in enumerate(neuron_graph.nodes()):
        nodes.SetPoint(ii,(neuron_graph.node[n]['x'],
                                 neuron_graph.node[n]['y'],
                                 neuron_graph.node[n]['z']))
        node_map[n] = ii
        radius.SetValue(ii, neuron_graph.node[n]['r'])
    
    cell_arr = vtk.vtkCellArray()
    colors = vtk.vtkUnsignedCharArray()
    colors.SetName('Colors')
    colors.SetNumberOfComponents(3)
    colors.SetNumberOfTuples(ncyl)    
    for ii, (n0, n1) in enumerate(neuron_graph.edges()):
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, node_map[n0])
        line.GetPointIds().SetId(1, node_map[n1])
        cell_arr.InsertNextCell(line)
        colors.SetTuple3(ii, *nodecolor[neuron_graph.node[n0]['s'] % len(nodecolor)])
        # radius.SetValue(ii, 0.5*(neuron_graph.node[n0]['r']+neuron_graph.node[n1]['r']))
    polydat = vtk.vtkPolyData()
    polydat.SetPoints(nodes)
    polydat.SetLines(cell_arr)
    mapper = vtk.vtkPolyDataMapper()
    if lines > 0:
        polydat.GetCellData().SetScalars(colors)
        if vtk.VTK_MAJOR_VERSION <= 5:
            mapper.SetInput(polydat)
        else:
            mapper.SetInputData(polydat)    
    else:
        polydat.GetPointData().AddArray(radius)
        polydat.GetPointData().SetActiveScalars(radius.GetName())
        polydat.GetPointData().AddArray(colors)
        tube_filt = vtk.vtkTubeFilter()
        tube_filt.SetNumberOfSides(10)
        tube_filt.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        if vtk.VTK_MAJOR_VERSION <= 5:        
            tube_filt.SetInput(polydat)
        else:
            tube_filt.SetInputData(polydat)
        mapper.SetInputConnection(tube_filt.GetOutputPort())    
        mapper.ScalarVisibilityOn()
        mapper.SetScalarModeToUsePointFieldData()
        mapper.SelectColorArray(colors.GetName())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    if lines > 0:
        actor.GetProperty().SetLineWidth(lines)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(*background)
    renderer.AddActor(actor)
    if len(label_nodes) > 0:
        label_actor = create_labels(neuron_graph, label_nodes, labels=labels, priorities=priorities)
        renderer.AddActor2D(label_actor)

    return renderer, actor
    

def neuron3d(neuron_graph, label_nodes=[], labels=[], priorities=[], nodecolor=nodecolor_4cp, background=(0,0,0), lines=False, stereo=False, axes=True, fullscreen=True):
    """axes: If true, show axis with scale bar info
    """
    renderer, actor = nrngraph2vtk(neuron_graph,
                                   label_nodes=label_nodes, labels=labels, priorities=priorities,
                                   nodecolor=nodecolor, background=background, lines=lines)
    if axes:  # show axis with scale bar info
        _ax = vtk.vtkCubeAxesActor2D()
        _ax.SetLabelFormat('%3.0f')
        _ax.SetNumberOfLabels(0)
        _ax.SetYAxisVisibility(False)
        _ax.SetZAxisVisibility(False)
        _ax.SetBounds(0, 200, -600, -400, 0, 200)
        _ax.SetXLabel('')
        _ax.SetXOrigin(0)
        _ax.SetYOrigin(-600)
        _ax.SetZOrigin(0)
        color = (1.0-background[0], 1.0-background[1], 1.0-background[2])
        tprop = vtk.vtkTextProperty()
        tprop.SetColor(color)
        _ax.SetAxisLabelTextProperty(tprop)
        _ax.GetProperty().SetColor(*color)
        _ax.SetFlyModeToClosestTriad()        
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

    # exporter = vtk.vtkX3DExporter()
    # exporter.SetFileName('ggn.x3d')
    # exporter.SetRenderWindow(win)
    # exporter.Write()
    
    interactor.Initialize()
    interactor.Start()


    
    
        
                        

# 
# morph3d_vtk.py ends here
