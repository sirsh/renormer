colours = ["black", "orange", "grey" , "blue", "green"]
species_labels = ["a", "b", "c", "d"]

from .incidence_matrix import incidence_matrix, _delta_
from .circuit import circuits,set_circuit_flow,set_source_sink_flow,apply_flow
from .viz import show_spanning_trees, ring_diagram, star,arc, debug_graph, tabulate_graphs,simple_ring_graph
from .propagators import scalar_propagator,pole_descriptor
from . import circuit
from . import perturbation 


import os
import numpy as np
def _edge_data_(self):
   
    vinf = self.shape[0] - 1
    def map_inf(e):
        return [
            e[0] if e[0] != vinf else -1,
            e[1] if e[1] != vinf else -1,
        ]
    v = np.array([map_inf(self.directed_edge_as_vertices(e)) for e in range(self.shape[-1])])
    w =np.array(self.species_vector).reshape(len(self.species_vector),1)
    return np.hstack([v,w])  

def save_graph_collection(col, file):
    with open(file, "w") as f:
        for s in col: f.write(str(_edge_data_(s).tolist())+ os.linesep)

def _graph_from_saved_format_(l):
    l = np.array(l)
    edges, species = l[:,:2].tolist(), l[:,-1].tolist()
    #return edges,species
    return incidence_matrix(edges=edges, species_vector=species)

def load_graph_collection(file, as_inc_matrix=True):
    import ast
    lst = []
    with open(file) as f:
        for l in f:
            l = l.rstrip('\n')
            if l != "":
                l = ast.literal_eval(l)
                if as_inc_matrix: l = _graph_from_saved_format_(l)
                lst.append( l )
    return lst

def convergent_diagrams_and_loops(col):
    """
    plot these with
    
    tabulate_graphs([loops[k]["graph"] for k in keys], diagram_class=simple_ring_graph)
    
    """
    convergent_diags = [c for c in col if c.num_ex_edges < 5]
    loops = {}
    for c in convergent_diags:   
        lh = c.loop_hash()
        if lh not in loops: loops[lh] = { "graph" : c.loop_part, "counter" : 0 }
        loops[lh]["counter"] += 1
    
    keys = list(loops.keys())   
    
    return convergent_diags, loops, keys