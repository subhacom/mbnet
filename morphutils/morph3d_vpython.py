# morph3d.py --- 
# 
# Filename: morph3d.py
# Description: 
# Author: Subhasis Ray
# Maintainer: 
# Created: Thu Jul 14 11:12:07 2016 (-0400)
# Version: 
# Package-Requires: ()
# Last-Updated: Tue Apr  9 20:45:48 2024 (+0530)
#           By: Subhasis Ray
#     Update #: 44
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
import vis
import networkx as nx

import bodies
import neurograph as ng

def neuron3d(neuron_graph, color=vis.color.white):    
    # vis.scene.stereo = 'redcyan'
    vis.scene.center = vis.vector((neuron_graph.nodes[1]['x'], neuron_graph.nodes[1]['y'], neuron_graph.nodes[1]['z']))
    vis.scene.width = 1024
    vis.scene.height = 768
    for n0, n1 in neuron_graph.edges():
        pos0 = np.array((neuron_graph.nodes[n0]['x'],
                neuron_graph.nodes[n0]['y'],
                neuron_graph.nodes[n0]['z']))
        pos1 = np.array((neuron_graph.nodes[n1]['x'],
                neuron_graph.nodes[n1]['y'],
                neuron_graph.nodes[n1]['z']))
        comp = vis.cylinder(pos=pos0, axis=pos1-pos0, radius=neuron_graph.nodes[n0]['r'], color=color, display=vis.scene)        

    
        

# 
# morph3d.py ends here
