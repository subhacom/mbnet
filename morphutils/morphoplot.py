# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 15:53:10 2016

@author: Subhasis
"""


from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d


def plot_3d_points(graph):
    """Plot the nodes of the tree as points in 3D"""
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for node, attr in graph.nodes(data=True):
        ax.plot([attr['x']], [attr['y']], [attr['z']], 'o') # This works - will need to plot cylinders for edges.
    return fig

    
def closest(G, x, y, z):
    """Return the node closest to the position x, y, z."""
    closest = -1
    dmin = 1e9
    for n, attr in G.nodes(data=True):
        d = np.sqrt((attr['x'] - x)**2 + (attr['y'] - y)**2 + \
            (attr['z'] - z)**2)
        if d < dmin:
            dmin = d
            closest = n
    return closest, dmin
    
    
def plot_3d_lines(graph, ax=None, color='k', alpha=0.7):
    """Plot the neuronal morphology tree in 3D
    if arrow is True use arrows from parent to child node.
    
    If color is a dict then each segment is color is looked up by the structure type attribute of its child node.
    """
    def show_node(event):
        ln = event.artist
        xdata = ln.get_xdata()
        ydata = ln.get_ydata()
        zdata = ln.get_zdata()
        ind = event.ind
        node, d = closest(graph, xdata[ind], ydata[ind], zdata[ind])
        print('onpick points:', node, d)
    
    if ax is None:
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.figure.canvas.mpl_connect('pick_event', show_node)
        
    for node, attr in graph.nodes(data=True):
        if isinstance(color, dict):
            c=color[attr['s']]
        else:
            c = color
        try:
            px = graph.nodes[attr['p']]['x']
            py = graph.nodes[attr['p']]['y']
            pz = graph.nodes[attr['p']]['z']
            ax.plot([attr['x'], px], [attr['y'], py], [attr['z'], pz], 
                    color=c, ls='-', alpha=alpha)            
        except KeyError:
            ax.plot([attr['x']], [attr['y']], [attr['z']], 'o', color=c,
                    alpha=alpha)            
    return ax
    
def mark_leaf_nodes(graph, ax, color='r', marker='o'):
    for node in graph.nodes():
        if graph.degree(node) == 1:
            ax.plot([graph.nodes[node]['x']], [graph.nodes[node]['y']], [graph.nodes[node]['z']], marker=marker, color=color)
            ax.text(graph.nodes[node]['x'], graph.nodes[node]['y'], graph.nodes[node]['z'], str(node))
    return ax


def plot_nodes(g, nodes, ax, color='r', marker='^'):
    """Mark the nodes of G plotted on axes ax.
    
    Useful for labeling specific nodes in two stages:

    ax = plot_3d_lines(g)
    ax = plot_nodes(g, nodes, ax)
    """
    for n in nodes:
        x = g.nodes[n]['x']
        y = g.nodes[n]['y']
        z = g.nodes[n]['z']
        ax.plot([x], [y], [z], marker=marker, color=color)
        ax.text(x, y, z, str(n))

        
def plot_nodes_2d(g, nodes, ax, label=False, proj='xy', color='r', marker='^',  alpha=1.0):
    for n in nodes:
        x = g.nodes[n][proj[0]]
        y = g.nodes[n][proj[1]]
        ax.plot([x], [y], marker=marker, color=color, alpha=alpha)
        if label:
            ax.text(x, y, str(n))

        
def plot_edges_2d(g, ax, label=False, proj='xy', color='r', marker='^', ls='-', alpha=0.8):
    for n1, n2 in g.edges():
        x1 = g.nodes[n1][proj[0]]
        y1 = g.nodes[n1][proj[1]]
        x2 = g.nodes[n2][proj[0]]
        y2 = g.nodes[n2][proj[1]]
        ax.plot([x1, x2], [y1, y2], marker=marker, ls=ls, color=color, alpha=alpha)
        if label:
            ax.text(x1, y1, str(n1))
        
        
        
