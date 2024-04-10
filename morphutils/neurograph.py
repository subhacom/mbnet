# -*- coding: utf-8 -*-
"""Data type definitions, maps and functions for morphology
processing.

Created on Thu Mar 17 15:45:19 2016

@author: Subhasis

"""

from __future__ import print_function

from collections import defaultdict
import itertools

import numpy as np
import networkx as nx


swc_dtype = np.dtype([('n', int), ('s', int),
                      ('x', float), ('y', float), ('z', float),
                      ('r', float), ('p', int)])


# The structure id is defined for 1-4, and anything above is custom
# Use neurite as the struct name for the custom ones
swc_sid = defaultdict(lambda: itertools.repeat('neurite').next)
swc_sid.update({
    0: 'unknown',
    1: 'soma',
    2: 'axon',
    3: 'dend_b',
    4: 'dend_a',
})

# These are custom segment types for GGN. It does not follow swc
# standard
custom_types = {
    'LCA': 5,
    'MCA': 6,
    'LH': 7,
    'alphaL': 8
    }

all_sid = {}   # Mapping from struct id to human readable name
all_sid.update(swc_sid)
for k, v in custom_types.items():
    all_sid[v] = k

# Mapping string names to numeric SID
name_sid = {v: k for k, v in all_sid.items()}


# This is colorblind safe and printer friendly dark color scheme named
# 3-class Dark2 from colorbrewer2.org
nodecolor_3cd2 = {
    0: (27, 158, 119),
    1: (217, 95, 2),
    2: (117, 112, 179),
}

nodecolor_3cs2 = {
    0: (102, 194, 165),
    1: (252, 141, 98),
    2: (141, 160, 203),
}

nodecolor_3cp = {
    0: (166, 206, 227),
    1: (31, 120, 180),
    2: (178, 223, 138)
}

# 4-class paired is the only color-blind safe and printer friendly
# qualitative colour scheme in colorbrewer.org


nodecolor_4cp = {
    0: (166, 206, 227),
    1: (31, 120, 180),
    2: (178, 223, 138),
    3: (51, 160, 44),
}

# 5-class dark2 is print friendly but not colorblind safe scheme from
# colorbrewer.org
nodecolor_5cd2 = {
    0: (27, 158, 119),
    1: (217, 95, 2),
    2: (117, 112, 179),
    3: (231, 41, 138),
    4: (102, 166, 30)
}

nodecolor_5cs3 = {
    0: (141, 211, 199),
    1: (255, 255, 179),
    2: (190, 186, 218),
    3: (251, 128, 114),
    4: (128, 177, 211)
}

# 10-class paired is dark scheme from colorbrewer.org
nodecolor_10cp = {
    0: (166, 206, 227),
    1: (31, 120, 180),
    2: (178, 223, 138),
    3: (51, 160, 44),
    4: (251, 154, 153),
    5: (227, 26, 28),
    6: (253, 191, 111),
    7: (255, 127, 0),
    8: (202, 178, 214),
    9: (106, 61, 154),
}

nodecolor_10cs3 = {
    0: (141, 211, 199),
    1: (255, 255, 179),
    2: (190, 186, 218),
    3: (251, 128, 114),
    4: (128, 177, 211),
    5: (253, 180, 98),
    6: (179, 222, 105),
    7: (252, 205, 229),
    8: (217, 217, 217),
    9: (188, 128, 189)
}

#: Alternative Scheme for Qualitative Data from
#: https://personal.sron.nl/~pault/
nodecolor_7q = { 
    0: (187, 187, 187),
    1: (102, 204, 238),
    2: (68, 119, 170),
    3: (170, 51, 119),
    4: (238, 102, 119),
    5: (204, 187, 68),
    6: (34, 136, 51)
}

#: 7-class paired from colorbrewer2.org
nodecolor_7cp = {  
    0: (166, 206, 227),
    1: (31, 120, 180),
    2: (178, 223, 138),
    3: (51, 160, 44),
    4: (251, 154, 153),
    5: (227, 26, 28),
    6: (253, 191, 111),
}

#: 7 class dark 2
nodecolor_7cd2 = { 
    0: (127, 158, 119),
    1: (217, 95, 2),
    2: (117, 112, 179),
    3: (231, 41, 138),
    4: (102, 166, 30),
    5: (230, 171, 2),
    6: (166, 118, 29)
}

nodecolor_7ca = {
    0: (127, 201, 127),
    1: (190, 174, 212),
    2: (253, 192, 134),
    3: (255, 255, 153),
    4: (56, 108, 176),
    5: (240, 2, 127),
    6: (191, 91, 23),
}

#: https://personal.sron.nl/~pault/colourschemes.pdf, SRON EPS technical note
nodecolor_9q = { 
    0 : (136, 204, 238),
    1 : (221, 204, 119),
    2 : (68, 170, 153),
    3 : (153, 153, 51),
    4 : (204, 102, 119),
    5 : (170, 68, 153),
    6 : (51, 34, 136),
    7 : (17, 119, 51),
    8 : (136, 34, 85),
}

#: colorblind friendly
nodecolor_15cb = { 
    0: (0, 0, 0),     #000000
    1: (0, 73, 73),     #004949   # soma
    2: (0, 146, 146),     #009292
    3: (255, 109, 182),     #ff6db6   # basal dendrite
    4: (255, 182, 219),     #ffb6db
    5: (219, 109, 0),     #db6d00 <- 12      # LCA
    6: (0, 109, 219),     #006ddb   # MCA
    7: (182, 109, 255),     #b66dff    # LH
    8: (109, 182, 255),     #6db6ff  # alphaL
    9: (182, 219, 255),     #b6dbff
    10: (146, 0, 0),     #920000
    11: (146, 73, 0),     #924900
    12: (73, 0, 146),     #490092
    13: (36, 255, 36),     #24ff24
    14: (255, 255, 109),     #ffff6d
}


colormaps = {'3cd2': nodecolor_3cd2,
             '3cs2': nodecolor_3cs2,
             '3cp': nodecolor_3cp,
             '4cp': nodecolor_4cp,
             '5cd2': nodecolor_5cd2,
             '5cs3': nodecolor_5cs3,
             '7q': nodecolor_7q,
             '7cp': nodecolor_7cp,
             '7ca': nodecolor_7ca,
             '7cd2': nodecolor_7cd2,
             '9q': nodecolor_9q,
             '15cb': nodecolor_15cb,
             '10cp': nodecolor_10cp,
             '10cs3': nodecolor_10cs3}


def swc2numpy(fname):
    """Read an swc file as a numpy record array.

    The record array has the following fields corresponding to swc:

    n: int
        numeric id of the point

    s: int
        structure id of the point, 0 - unknown, 1 - soma, 2 - axon,
        3 - basal dendrite, 4 - apical dendrite, 5 - custom

    x, y, z, r: float
        x, y and z position and radius

    p: int
        numeric id of the parent node

    """
    data = np.loadtxt(fname, dtype=swc_dtype)
    return data


def eucd(G, n1, n2):
    """Eucledian distance between n1 and n2.

    Assumes there exist attributes x, y and z for each node.
    """
    return np.sqrt((G.nodes[n1]['x'] - G.nodes[n2]['x'])**2
                   + (G.nodes[n1]['y'] - G.nodes[n2]['y'])**2 +
                   (G.nodes[n1]['z'] - G.nodes[n2]['z'])**2)


def swc2graph(source):
    """Convert an SWC file into a networkx DiGraph.

    Parameters

    source: str or numpy.recarray.
            If str, it is the path of the swc file. If recarray, it should
            have swc_dtype as the dtype.


    Returns

    networkx.DiGraph: Each node corresponds to a point in the trace, with same
                      attributes as swc_dtype (other than `n`). Each edge is
                      assigned the attribute `length` - the euclidean distance
                      bewteen the nodes it connects.

    Each compartment has an edge to each of its children.

    """
    if isinstance(source, str):
        data = swc2numpy(source)
        print('Loaded ', source)
    else:
        data = source
    g = nx.DiGraph()
    for row in data:
        if row['p'] >= 0:
            g.add_edge(row['p'], row['n'])
        else:
            g.add_node(row['n'])
        g.nodes[row['n']]['x'] = row['x']
        g.nodes[row['n']]['y'] = row['y']
        g.nodes[row['n']]['z'] = row['z']
        g.nodes[row['n']]['r'] = row['r']
        g.nodes[row['n']]['p'] = row['p']
        g.nodes[row['n']]['s'] = row['s']
    # Store physical length of each edge
    # UPDATE: networkx2.0 removed Graph.edges_iter
    if nx.__version__.startswith('1.') or nx.__version__.startswith('0.'):
        for n0, n1 in g.edges_iter():
            g.edge[n0][n1]['length'] = eucd(g, n0, n1)
    else:
        for n0, n1 in g.edges:
            g.edges[n0, n1]['length'] = eucd(g, n0, n1)
    return g


def graph2swc(G, filename):
    """Save the morphology in graph G as an swc file"""
    with open(filename, 'w') as fd:
        for n in sorted(G.nodes()):
            try:
                fd.write('{} {s} {x:.6f} {y:.6f} {z:.6f} {r:.6f} {p}\n'.format(
                    n, **G.nodes[n]))
            except KeyError as e:
                print('Error with node', n, ':', e)
                print('Node attributes:', G.nodes[n])

                
def numpy2swc(swc, filename):
    """Convert numpy array into SWC.

    Parameters
    ----------
    swc: np.ndarray or np.recarray
        The columns in swc array (recarray) must be:
        id, type, x, y, z, r, parentid
    filename: path of the SWC file to save data in.
    """
    np.savetxt(filename, swc, '%d %d %.6f %.6f %.6f %.6f %d')
    

def sorted_edges(G, attr='length', reverse=False):
    """Sort the edges by attribute `attr`, defaults to `length`"""
    if nx.__version__.startswith('1.') or nx.__version__.startswith('0.'):
        return sorted(G.edges_iter(), key=lambda x: G[x[0]][x[1]][attr],
                      reverse=reverse)
    return sorted(list(G.edges), key=lambda x: G.edges[x[0], x[1]][attr],
                  reverse=reverse)


def branch_points(G):
    """Retrieve a list of branch points and number of branches at those
    points"""
    if nx.__version__.startswith('1.') or nx.__version__.startswith('0.'):
        return [(node, degree) for node, degree in G.out_degree().items()
                if degree > 1]
    else:
        return [(node, degree) for node, degree in G.out_degree()
                if degree > 1]


def n_branch(G):
    """Number of branches in cell graph G"""
    # +1 For the tree trunk from soma
    return sum([d for n, d in branch_points(G)]) + 1 


def get_stype_node_map(g):
    """Return a reverse lookup table `{sid: list of nodes}` of nodes by
    stype for neuron-graph `g`
    """
    ret = defaultdict(list)
    for n, d in g.nodes(data=True):
        ret[d['s']].append(n)
    return ret


def select_random_nodes_by_sid(g, sid_list, counts, replace=False):
    """For each sid in `sid_list` select `count` random nodes from g - a
    neurongraph where `count` is the corresponding entry in `counts`
    list.

    Returns a dict of selected nodes mapped to the sids.

    @seealso: nrnutils.select_random_nodes_by_sid - to be deprecated

    """
    type_node_map = get_stype_node_map(g)
    synnodes = select_random_items_by_key(type_node_map, sid_list,
                                          counts, replace)
    return synnodes


def select_random_items_by_key(listdict, keys, counts, replace=False):
    """Make a dict of random selections from the dict items.

    listdict: {k: L} a dict mapping keys to lists

    keys: list of keys for which we select entries from L.

    counts: list of numbers specifying how many elements to pick from
    L for the corresponding key in `keys`. If the entry is bigger than
    the list length, then the entire list is taken instead of throwing
    error (something that numpy.random.choice does).

    replace: select with replacement (default False)

    returns: dict {k: L2} where L2 is a list containing a random
    selection of items from {k: L} with size specified by
    corresponding entry in counts.

    """
    ret = {}
    for key, count in zip(keys, counts):
        size = len(listdict[key])
        if count < size:
            size = count
        ret[key] = list(np.random.choice(listdict[key],
                                         size=size, replace=False))
    return ret


def remove_null_edges(G, n0, n1=None, attr='length', lim=0.1,
                      rtol=1e-5, atol=1e-8):
    """Recursively remove edges in G starting with [n0->n1] if the edge has
    attribute attr=0.

    WARNING: it modifies G itself and the nodes should renumbered.
    """
    for n in list(nx.all_neighbors(G, n0)):
        if n != G.nodes[n0]['p']:
            remove_null_edges(G, n, n0, rtol=rtol, atol=atol)
    if n1 is None:
        return
    if np.isclose(G[n0][n1][attr], 0.0, rtol=rtol, atol=atol):
        for n in list(nx.all_neighbors(G, n0)):
            if n not in G.neighbors(n0):
                G.nodes[n]['p'] = n1
                G.add_edge(n, n1, G[n][n0])
        G.remove_node(n0)


def remove_shorter_edges(G, n0=1, lim=0.1, verbose=False):
    """Remove edges in G starting with n0 if the edge is shorter than lim.

    WARNING: it modifies G itself and the nodes should renumbered.

    """
    todo = {n0: None}
    while todo:
        keys = list(todo.keys())
        for n0 in keys:
            todo.pop(n0)
            for n in list(nx.neighbors(G, n0)):
                if G[n0][n]['length'] < lim:
                    if verbose:
                        print('removing', n0, '->', n, 'of length',
                              G[n0][n]['length'])
                    for n2 in list(G.neighbors(n)):
                        G.nodes[n2]['p'] = n0
                        G.add_edge(n0, n2, length=G[n][n2]['length'])
                    G.remove_node(n)
                    todo[n2] = None
                else:
                    todo[n] = None


def to_undirected_nrn(G):
    """Workaround to convert graph with nodes whose attribute 'orig'
    is neuron section"""
    G_ = G.__class__()
    G_.add_nodes_from(G)
    for edge in G.edges():
        print('#', edge)
        G_.add_edge(*edge)
    for node in G.nodes:
        G_.nodes[node].update(G.nodes[node])
        if G.nodes[node]['orig'] is not None:
            G_.nodes[node]['orig'] = G.nodes[node]['orig'].name()
        else:
            G_.nodes[node]['orig'] = None
    G2 = G_.to_undirected()
    for node in G.nodes:
        G2.nodes[node]['orig'] = G.nodes[node]['orig']
    return G2
    
                    
def remove_longer_edges(G, lim=100.0, verbose=False):
    """Delete the edges which are longer than lim. Returns

    c, e

    where c is a list of connected components. If no edges were
    removed, a list containingc contains the original graph G.

    e is a list of node pairs defining the removed edges.

    """
    edges = sorted_edges(G, reverse=True)
    long_edges = []
    for n0, n1 in edges:
        if G[n0][n1]['length'] > lim:
            if verbose:
                print('long edge: {} -- {}: {} um'.format(
                    n0, n1, G[n0][n1]['length']))
            long_edges.append((n0, n1))
        else:
            # edges are already sorted in descending order, skip shorter edges
            break
    components = []
    # NEURON7.7/NetworkX 3.x sim incompatible - sections cannot be deepcopied as they cannot be pickled
    G_ = G.__class__()
    G_.add_nodes_from(G)
    G_.add_edges_from(G)
    for node in G.nodes:
        G_.nodes[node]['orig'] = G.nodes[node]['orig'].name()
    if len(long_edges) > 0:
        G2 = G_.to_undirected()
        for (n0, n1) in long_edges:
            G2.remove_edge(n0, n1)
        for node in G2:
            G2.nodes[node]['orig'] = G.nodes[node]['orig']
        for sub in nx.connected_component_subgraphs(G2):
            components.append(nx.DiGraph(G.subgraph(sub.nodes())))
    else:
        components = [G]
    if verbose:
        print('Removed long edges:', len(long_edges))
    return components, long_edges


def cleanup_morphology(G, start=1, lmin=0.1, rmin=0.1, rdefault=0.5,
                       verbose=False):
    """Remove edges that are shorter than lmin and with radius less than
    rmin.

    """
    Gtmp = G.copy()
    remove_shorter_edges(Gtmp, start, lim=lmin, verbose=verbose)
    update_thin_segs(Gtmp, start, lim=rmin, default=rdefault, verbose=verbose)
    Gnew, nmap = renumber_nodes(Gtmp, start=start)
    return Gnew, nmap


def renumber_nodes(G, start=1, undirected=True):
    """Renumber nodes so that NeuronLand converter does not mess up.

    undirected: bool, if True convert G to undirected graph
    first. This helps avoid problems when the order of the
    parent-child nodes is incorrect but connections are fine.
        

    Returns: (new_graph, node_map)

    node_graph: a new DiGraph with the nodes numbered by bfs sequence.

    node_map: dict mapping node number in G to that in node_graph

    """
    ret = nx.DiGraph()
    node_map = {}
    if undirected:
        G = to_undirected_nrn(G)
    # The nodes are connected as child->parent, we are starting from root,
    # hence reverse    
    for ii, (n1, n2) in enumerate(nx.dfs_edges(G, start)):
        m1 = ii+1
        m2 = ii+2
        if n1 in node_map:
            m1 = node_map[n1]
        else:
            node_map[n1] = m1
            ret.add_node(m1)
            ret.nodes[m1].update(G.nodes[n1])
        node_map[n2] = m2
        ret.add_edge(m1, m2, length=G[n1][n2])
        ret.nodes[m2].update(G.nodes[n2])
        ret.nodes[m2]['p'] = m1
    ret.nodes[1]['p'] = -1
    return ret, node_map


def update_thin_segs(G, n=1, lim=0.1, default=None, verbose=False):
    """Change segments with r < lim to have r = (max(r_children) +
    r_parent)/2. If there is no parent node, then max(r_children) is
    the radius. This is somewhat arbitrary.

    """
    if default is None:
        default = lim
    stack = [n]
    while len(stack) > 0:
        n = stack[-1]
        stack.pop()
        attr = G.nodes[n]
        if attr['r'] < lim:
            rcmax = 0.0
            for c in list(G.neighbors(n)):
                rc = G.nodes[c]['r']
                if rcmax < rc:
                    rcmax = rc
            parent = attr['p']
            if parent > 0:
                rp = G.nodes[parent]['r']
            else:
                rp = rcmax
            r = (rp + rcmax) * 0.5
            if r > lim:
                G.nodes[n]['r'] = r
            else:
                G.nodes[n]['r'] = default
                if verbose:
                    print('Node {} radius set to default {}.'.format(
                        n, default))
            stack += list(G.neighbors(n))


def soma_distance(g):
    """Physical distance along the path from soma to each node"""
    soma = 1
    distdict = {soma: 0.0}
    for n0, n1 in nx.dfs_edges(g, soma):
        if n1 not in distdict:
            distdict[n1] = distdict[n0] + g.edge[n0][n1]['length']

    return distdict


def soma_pathlen(g):
    """The number of edges between soma and each node"""
    soma = 1
    distdict = {soma: 0}
    for n0, n1 in nx.dfs_edges(g, soma):
        if n1 not in distdict:
            distdict[n1] = distdict[n0] + 1
    return distdict


def soma_branch_len(g):
    """Number of branchpoints between each node and soma."""
    soma = 1
    distdict = {soma: 0}
    for n0, n1 in nx.dfs_edges(g, soma):
        if n1 not in distdict:
            if nx.degree(g, n0) > 1:
                dbranch = 1
            else:
                dbranch = 0
            distdict[n1] = distdict[n0] + dbranch

    return distdict


def eleclen(g, rm_sp=1000.0, ra_sp=100.0, avg_dia=False, inplace=False):
    """Compute electrotonic length of each edge in `g` assuming specific
    membrane resistance rm_sp (default 1000 Ohm-cm2) and specific
    axial resistance ra_sp (default 100 Ohm-cm).

    if avg_dia is True, each edge is considered a cylinder with
    diameter that is average of the two nodes. Otherwise the total
    membrane resistance Rm and total axial resistance Ra are computed
    for the frustum (cut cone) and then the electrotonic length
    constant of the compartment is computed as:

    lambda = l_physical * sqrt(Rm / Ra).

    The electrotonic length is: l_physical / lambda = sqrt(Ra / Rm)

    Returns a copy of g with added edge attribute `L` containing
    this value.

    You can use networkx.shortest_path_length() to obtain the
    electrotonic distance between node pairs in the same format, i.e.,

    if source and target are specified, returns the electrotonic
    distance L between source and target.

    if only the source is specified, returns a dict {target: L} for
    all possible target nodes.

    if only the target is specified, returns a dict {source: L} for
    all possible source nodes.

    if neither the source nor the target is specified, returns a dict
    {source: {target: L}}

    """
    g2 = g if inplace else g.copy()
    rratio = 1e4 * rm_sp / ra_sp    # 1e4 to convert cm to um
    for n0, n1 in nx.dfs_edges(g2):
        rad = ((0.5 * (g2.nodes[n0]['r'] + g2.nodes[n1]['r']))
               if avg_dia else g2.nodes[n0]['r'])
        g2[n0][n1]['L'] = g2[n0][n1]['length'] / np.sqrt(rad * rratio)
    return g2


def join_neurites(left, leftnode, right, rightnode, leftroot=None,
                  translate=False):
    """Join `leftnode` of cellgraph `left` with `rightnode` of cellgraph
    `right`.

    translate: if True, `right` will be translated to make `leftnode` and `rightnode`
    overlap.
    All the nodes will be renumbered.
    """
    dn = max(left.nodes()) + 1
    label_map = {node: node + dn for node in right.nodes()}
    relabeled = nx.relabel_nodes(right, label_map)
    if translate:
        dx, dy, dz = (left.nodes[leftnode]['x'] - right.nodes[rightnode]['x'],
                      left.nodes[leftnode]['y'] - right.nodes[rightnode]['y'],
                      left.nodes[leftnode]['z'] - right.nodes[rightnode]['z'])
    else:
        dx, dy, dz = 0, 0, 0
    for node in relabeled:
        relabeled.nodes[node]['p'] += dn
        relabeled.nodes[node]['x'] += dx
        relabeled.nodes[node]['y'] += dy
        relabeled.nodes[node]['z'] += dz
    relabeled.nodes[label_map[rightnode]]['p'] = leftnode
    combined = nx.union(left, relabeled)
    dl = 0 if translate else np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
    combined.add_edge(leftnode, label_map[rightnode], length=dl)
    if leftroot is None:
        leftroot = min(left.nodes())
    return renumber_nodes(combined, leftroot)[0]


def length(cellgraph):
    """Total length of the neurites"""
    return sum([cellgraph.edges[n0, n1]['length']
                for n0, n1 in cellgraph.edges])


def n_leaves(cellgraph):
    """Number of leaf nodes"""
    return len([n for n in cellgraph.nodes()
                if (cellgraph.out_degree(n) == 0) and
                (cellgraph.in_degree(n) == 1)])


# def union_cellgraphs(*args):
#     start = 1
#     combined_cellgraph = nx.DiGraph()
#     for ii, cellgraph in enumerate(args):
#         if ii > 0:
#             maxnode = max(combined_cellgraph.nodes())
#             for n, attr in cellgraph.nodes(data=True):
#                 attr['p'] += maxnode
#             for jj, (n1, n2) in enumerate(cellgraph.edges
#             cellgraph = nx.relabel_nodes(cellgraph, mapping)
#         combined_cellgraph = nx.union(combined_cellgraph, cellgraph)
