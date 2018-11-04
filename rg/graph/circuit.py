#temp: later pick your objects - just is just to prototype
from sympy import *

import operator
import functools
import math
import pandas as pd
import numpy as np

#work in progress - how do we safely set the flows the way we want on all diagrams
def apply_flow(inc_mat, func, optional_tree_to_do=None):
    """
    apply all flows in the diagram given the supplied func from this module 
    e.g. 
    apply_flow(sunset,set_circuit_flow) # will set the circut flow
    apply_flow(sunset,set_source_sink_flow) # will set the source to sink flow
    
    and these *should* be done at every level but I am still working on the theory for this
    
    an optional callback could do something like determine the omega signs but need to be careful
    
    """
    #requires testing - this can work by change but generally I think the rule is to always start the cicuits by something already in the tree
    #that way any permutations for edge ordering will not break the path and will continue whatever flow exists
    #the source and sink will surely be respected by something in that edge ordering 
    
    def _permute_circuit_start_with_existing(c,existing=[]):
        if len (existing) ==0:return c
        while c[0] not in existing:
            h = c.pop(0)
            c = c + [h]
        return c

    all_circuits = circuits(inc_mat,inc_mat.spanning_trees.iloc[0].values)
    existing_edges = []
    for c in all_circuits:
        c = _permute_circuit_start_with_existing(c,existing_edges)
        existing_edges = list(set(c + existing_edges))
        res = func(inc_mat,c)
        #print("source,sink",res)

def _ensure_edge_path_(edges, edge_ids, h=[], edge_permutation=[]):
    """there should be a path from vertex to vertex in the correct order - we need to chain the vertices along the edges
        eg = [[0, 1], [2, 1], [4, 3], [0, 3], [4, 2]]
       _ensure_edge_path_(eg)  =>  [[0,1],[1,2],[2,4],[4,3],[3,0]]
    """
    #copy the list
    edge_ids = list(edge_ids)
    def _match_(h,e):
        if e[0] == h[-1]:return e
        if e[1] == h[-1]: return list(reversed(e))
        return None
    
    if len(edges)>0:
        if len(h) == 0:
            h = [edges.pop(0)]
            edge_permutation = [edge_ids.pop(0)]
            return _ensure_edge_path_(edges,edge_ids,h,edge_permutation)
        for i,e in enumerate(edges):
            _e = _match_(h[-1],e)
            if _e is not None:
                edges.pop(i)
                h.append( _e )
                edge_permutation.append(edge_ids.pop(i))
                return _ensure_edge_path_(edges,edge_ids,h,edge_permutation)
    return h,edge_permutation




def set_circuit_flow(inc_mat, circ, entrainment={} ):
    circ_evs = [inc_mat.directed_edge_as_vertices (e) for e in circ]
    circ_evs,edge_permutation = _ensure_edge_path_(circ_evs,edge_ids=circ)
    inc_mat.ensure_edge_directions(dict(zip(edge_permutation, circ_evs)))
    
def set_source_sink_flow(inc_mat, circ, known_source_sinks = [], maxout=100):
    adj = _adj_struct_from_edges_(inc_mat,circ)
    source,sink = tuple(adj[adj["e"].isin(circ)][["v","d"]].drop_duplicates().sort_values("d")[:2]["v"].values)
    #check if we already have a known source and sink and treat it as a negative score - if we choose our own, add to the list
    
    #partition the circuit between source and sink
    circ_evs = [inc_mat.directed_edge_as_vertices (e) for e in circ]
    #print("before",circ_evs)
    circ_evs,circ = _ensure_edge_path_(circ_evs,circ)   
    #print("after", circ_evs)
    S,s,cursor,p1 = -1,-1,0,[]
    #print("source, sink:",source,sink)
    while s == -1:
        index = cursor%len(circ_evs)
        v = circ_evs[index]
        #print(index, v)
        if v[0] == source:S = index
        if v[0] == sink and S != -1:s = index  
        if S != -1 and s == -1:p1.append(index)
        cursor+= 1
        if cursor == maxout:return
        
    p2 = list(set(range(len(circ))) - set(p1))
    #return [circ_evs[e] for e in p1],[list(reversed(circ_evs[e])) for e in p2]
    #keep one side and reverse the other
    new_edges = [circ_evs[e] for e in p1] + [list(reversed(circ_evs[e])) for e in p2]
    edge_ids = [circ[i] for i in p1+p2]
    d = dict(zip(edge_ids, new_edges))
    #print(d)
    inc_mat.ensure_edge_directions(d)
    #print(circ)
    return source,sink
    
    
# def _ensure_edge_path2_(edges, h=[]):
#     """there should be a path from vertex to vertex in the correct order - we need to chain the vertices along the edges
#         eg = [[0, 1], [2, 1], [4, 3], [0, 3], [4, 2]]
#        _ensure_edge_path_(eg)  =>  [[0,1],[1,2],[2,4],[4,3],[3,0]]
#     """
#     def _match_(h,e):
#         if e[0] == h[-1]:return e
#         if e[1] == h[-1]: return list(reversed(e))
#         return None
    
#     if len(edges)>0:
#         if len(h) == 0:
#             h = [edges.pop(0)]
#             return _ensure_edge_path_(edges,h)
#         #pop the next guy
#         e = edges.pop(0)
#         #make sure she is the right way around
#         print(h,e)
#         h.append( _match_(h[-1],e) )
#         return _ensure_edge_path_(edges,h)
    
#     return h

def _adj_struct_from_edges_(inc_mat, eset=None):
    tuples = []
    if eset == None: eset = list(range(inc_mat._num_internal_edges))
    for e in eset:
        t = list(np.where(inc_mat[:,e]!=0)[0])
        tuples.append(t+[e])
        tuples.append(list(reversed(t))+[e])
    df = pd.DataFrame(tuples,columns=["a","b", "e"]).sort_values("a")
    tds = pd.DataFrame([[k,v] for k,v in inc_mat.tree_distance_map.items()],columns=["v", "d"])
    return pd.merge(df,tds, left_on="a", right_on="v")

def random_walk(adj, nedges,i=0,j=0, max_length=100, min_path_length=1, set_edge_sense=None):
    """Random walk on edge boolean space until we construct a path of length greater than 1 using a clever backtrack trick"""
    
    l = 0 #safety
    epath = np.zeros(nedges, bool)
    _i = i   
    def _random_choice_(adj,i):
        mat = adj[adj["a"]==i][["e", "b"]].as_matrix()
        return tuple(mat[np.random.randint(len(mat))])
    ###################################################
    loop_formed = lambda : epath.sum() > min_path_length and _i == i 
    while l < max_length and not loop_formed():
        e, _i = _random_choice_(adj, _i )
        #print(e,_i)
        epath[e] = np.logical_not(epath[e])
        l+=1   
    if not loop_formed():epath=False 
    ####################################################
    #return epath
    return list(np.where(epath)[0])
    
def circuits(inc_matrix, spanning_tree=None,sort=True):  
    if spanning_tree is None:spanning_tree = inc_matrix.spanning_trees.iloc[0].values
    st_comp = inc_matrix.edge_complement(spanning_tree)
    #this spanning tree must have L missing edges
    all_circuits = []
    for e in st_comp:#add one edge and determine it's incident vertex to start a circuit trace
        ev = inc_matrix.directed_edge_as_vertices(e)
        #print("edge",e, "is", ev)
        sub_tree = [e] + list(spanning_tree)
        #print(sub_tree)
        adj = _adj_struct_from_edges_(inc_matrix,sub_tree)
        circuit = list(random_walk(adj,len(inc_matrix.edges), ev[0],ev[0]))
        all_circuits.append(circuit)   
        
    def circuit_sorter(inc_mat):
        graph = _adj_struct_from_edges_(inc_mat)
        def _sorter_(a):
            distinct_vertices = np.unique(np.array([inc_mat.directed_edge_as_vertices(e) for e in a ]))
            return graph[graph["a"].isin(distinct_vertices)][["a","d"]].drop_duplicates()["d"].mean()
        return _sorter_
    
    if sort:
        s = circuit_sorter(inc_matrix)
        all_circuits.sort(key=lambda x : s(x))
        
    return all_circuits


# class cycle_finder(object):
#     """
#     Adapted from
#     https://stackoverflow.com/questions/12367801/finding-all-cycles-in-undirected-graphs
#     See also
#     http://dspace.mit.edu/bitstream/handle/1721.1/68106/FTL_R_1982_07.pdf
#     #graph = [[1, 2], [1, 3], [1, 4], [2, 3], [3, 4], [2, 6], [4, 6], [8, 7], [8, 9], [9, 7]]
#     graph = [[0,1], [0,2], [1,2], [1,3], [2,3]]

#     for p in cycle_finder(graph): print(p) #returns the edges in the circuit indexed from the graph
#     for p in cycle_finder(graph, "vertices"): # returns the vertices in the cicuit

#     """
#     def __init__(self, graph_edges, c_type="edges"):
#         #i need them sorted here - this is an undirected graph too
#         self.graph = graph_edges#[sorted(l) for l in list(graph_edges)]
#         self.cycles = []  
#         self.c_type = c_type
#         for edge in self.graph:
#             for node in edge:  
#                 self.findNewCycles([node])
    
#     def __expand_edges__(self, vlist):
#         for i,v in enumerate(vlist): 
#             yield sorted([vlist[i-1],v])
            
#     @property
#     def cycle_basis(self):
#         l = list(self)
#         b = np.zeros((len(l),len(self.graph)),np.int)
#         for i,r in enumerate(b): r[l[i]] = 1
#         return b
    
#     def __iter__(self):
#         for cy in self.cycles:
#             res = [node for node in cy]
#             if self.c_type == "vertices": yield res
#             else: yield [self.graph.index(e) for e in self.__expand_edges__(res)]
               
#     def findNewCycles(self,path):
#         def isNew(path): return not path in self.cycles
#         def visited(node, path):  return node in path
#         def invert(path):return rotate_to_smallest(path[::-1])
#         def rotate_to_smallest(path):
#             n = path.index(min(path))
#             return path[n:]+path[:n]   
        
#         start_node,next_node,sub = path[0],None,[]
        
#         for edge in self.graph:
#             node1, node2 = edge
#             if start_node in edge:
#                     next_node = node2 if node1 == start_node else node1
#             if not visited(next_node, path):
#                     sub = [next_node]
#                     sub.extend(path)
#                     self.findNewCycles(sub);
#             elif len(path) > 2  and next_node == path[-1]:
#                     p = rotate_to_smallest(path);
#                     inv = invert(p)
#                     if isNew(p) and isNew(inv):  
#                         self.cycles.append(p)
                        