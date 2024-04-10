# vtk_dump_movie.py ---
#
# Filename: vtk_dump_movie.py
# Description:
# Author: Subhasis Ray
# Maintainer:
# Created: Mon Nov  7 16:56:31 2016 (-0500)
# Version:
# Package-Requires: ()
# Last-Updated: Tue Apr  9 20:44:05 2024 (+0530)
#           By: Subhasis Ray
#     Update #: 198

# Code:

from __future__ import print_function
import sys
import numpy as np
import neurograph as ng
import vtk
from morph3d_vtk import nrngraph2vtk
from displaycell import make_tmatrix, rigid_transform_graph

import argparse

def dump_movie(filename, neuron_graph, label_nodes=[], labels=[],
               priorities=[], nodecolor=ng.nodecolor_4cp,
               background=(0,0,0), lines=False,
               xrot=0.0, yrot=0.0, zrot=0.0,
               xangle=0, yangle=0, zangle=0,
               frames_per_degree=10, framerate=25, size=(800, 600)):
    """Save the 3D display as a movie"""
    x = [neuron_graph.nodes[n]['x'] for n in neuron_graph]
    y = [neuron_graph.nodes[n]['y'] for n in neuron_graph]
    z = [neuron_graph.nodes[n]['z'] for n in neuron_graph]
    center_x = (min(x) + max(x)) * 0.5
    center_y = (min(y) + max(y)) * 0.5
    center_z = (min(z) + max(z)) * 0.5
    tmatrix = make_tmatrix(-center_x, -center_y, -center_z, xrot, yrot, zrot)
    rigid_transform_graph(neuron_graph, tmatrix)
    # Translate origin to center
    renderer, actor = nrngraph2vtk(neuron_graph,
                                   label_nodes=label_nodes, labels=labels, priorities=priorities,
                                   nodecolor=nodecolor, background=background, lines=lines)
    x = [neuron_graph.nodes[n]['x'] for n in neuron_graph]
    y = [neuron_graph.nodes[n]['y'] for n in neuron_graph]
    z = [neuron_graph.nodes[n]['z'] for n in neuron_graph]
    center_x = (min(x) + max(x)) * 0.5
    center_y = (min(y) + max(y)) * 0.5
    center_z = (min(z) + max(z)) * 0.5
    print('Center:', center_x, center_y, center_z)
    actor.SetOrigin(center_x, center_y, center_z)
    win = vtk.vtkRenderWindow()
    win.SetSize(*size)
    win.AddRenderer(renderer)

    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(win)

    win.Render()
    interactor.Initialize()
    xframes = xangle * frames_per_degree
    dx = 0 if xframes == 0 else 1.0 / frames_per_degree
    yframes = yangle * frames_per_degree
    dy = 0 if yframes == 0 else 1.0 / frames_per_degree
    zframes = zangle * frames_per_degree
    dz = 0 if zframes == 0 else 1.0 / frames_per_degree
    frames = int(max((xframes, yframes, zframes)))
    print('X-rotation', dx, 'Y-rotation', dy, 'Z-rotation', dz, 'Frames', frames)

    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(win)
    # windowToImageFilter.SetInputBufferTypeToRGBA()
    windowToImageFilter.ReadFrontBufferOff()
    windowToImageFilter.Update()
    writer = vtk.vtkAVIWriter()
    writer.SetRate(framerate)
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.SetFileName(filename)
    writer.Start()
    for frame in range(frames):
        actor.RotateX(dx)
        actor.RotateY(dy)
        actor.RotateZ(dz)
        interactor.GetRenderWindow().Render()
        windowToImageFilter.Modified() # This is crucial
        writer.Write()
    writer.End()
    interactor.Start()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', type=str,
                        help='swc file to be rendered',
                        dest='infile', required=True)
    parser.add_argument('--xrot', type=float,
                        default=0.0,
                        help='rotation around x axis (before movie start)')
    parser.add_argument('--yrot', type=float,
                        default=0.0,
                        help='rotation around y axis (before movie start)')
    parser.add_argument('--zrot', type=float,
                        default=0.0,
                        help='rotation around z axis (before movie start)')
    parser.add_argument('-x', type=float,
                        help='angle (in degrees) for x rotation',
                        dest='xangle', default=0.0)
    parser.add_argument('-y', type=float,
                        help='angle (in degrees) for y rotation',
                        dest='yangle', default=0.0)
    parser.add_argument('-z', type=float,
                        help='angle (in degrees) for z rotation',
                        dest='zangle', default=0.0)
    parser.add_argument('-f', '--frames', type=int,
                        help='number of frames per degree',
                        dest='fpd', default=100)
    parser.add_argument('-r', '--rate', type=int,
                        help='frame rate',
                        dest='framerate', default=25)
    parser.add_argument('-o', '--output', type=str,
                        help='output filename', dest='outfile',
                        required=True)
    args = parser.parse_args()
    infile = args.infile
    outfile = args.outfile
    fpd = args.fpd
    framerate = args.framerate
    print(args.xangle, args.yangle, args.zangle, fpd, framerate)
    graph = ng.swc2graph(infile)

    dump_movie(outfile, graph, xrot=args.xrot, yrot=args.yrot, zrot=args.zrot,
               xangle=args.xangle, yangle=args.yangle, zangle=args.zangle,
               frames_per_degree=fpd, framerate=framerate)
    print('Saved movie in', outfile)


#
# vtk_dump_movie.py ends here
