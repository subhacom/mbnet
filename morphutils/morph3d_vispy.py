# morph3d.py --- 
# 
# Filename: morph3d.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Jul 14 11:12:07 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Tue Apr  9 20:45:31 2024 (+0530)
#           By: Subhasis Ray
#     Update #: 92
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
"""Attempt at displaying 3D morphology using VPython"""

import numpy as np
from vispy import scene as scene
from vispy.scene import visuals as vis
import networkx as nx

import neurograph as ng

def neuron3d(neuron_graph, color='blue', name='Neuron'):
    canvas = scene.SceneCanvas(title=name)
    view = canvas.central_widget.add_view()
    camera = scene.cameras.TurntableCamera(fov=45, azimuth=-45, parent=view.scene)
    view.camera = camera
    comps = []
    for n0, n1 in neuron_graph.edges():
        pos0 = np.array((neuron_graph.nodes[n0]['x'],
                neuron_graph.nodes[n0]['y'],
                neuron_graph.nodes[n0]['z']))
        pos1 = np.array((neuron_graph.nodes[n1]['x'],
                neuron_graph.nodes[n1]['y'],
                neuron_graph.nodes[n1]['z']))
        mid = (pos0 + pos1)/2.0
        comp = vis.Tube(name='{}-{}'.format(n0, n1), points=np.vstack((pos0, mid, pos1)), radius=neuron_graph.nodes[n0]['r'], color=color)
        view.add(comp)
        comps.append(comp)
    canvas.show()
    canvas.app.run()



# 
# morph3d.py ends here
