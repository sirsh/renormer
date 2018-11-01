#this should be top level one so noone should import it
#the main innovation i can make is use the global representation in the graph theoretic sense for the graphs
#and then have a vetex collection which shows the tensor decompositions in terms of residual and internal half edges
#i can just keep 2 of three of these set types e.g. residual and internal (which can be summed)
#graph.vertex_residuals, graph.vertex_flags

#todo - helpers for incident matrix such as enumerate things we care about e.g. external edges, external vertices etc.
#if graphs have extra edge meta data for species type, field type, they can be slwoly replace interaction

from sympy import latex,init_printing,Matrix,collect,Symbol,expand,solve
import numpy as np
from functools import reduce
from operator import mul,add

import itertools

from . propagator import *
#watch circular
from . interaction import interaction as J, composite_interaction as G
from . diagrams import diagram_set

class composite_interaction_graph(object):
    def __init__(self, comp_interaction):
        self.ci = comp_interaction
        temp=self.ci.graph_dataframes()[1]
        self._edges = [d["edge"] for d in self.ci.edges]
        #todo - get the external vertices which will be connected to the infinity node
        self._externals = list(self.__gen_ex_vert__())
        
        #vids are the non 0 source or tarket in the external edges dataframe
        l = list(zip( temp[temp.internal==False].source, temp[temp.internal==False].target))
        l = [ l_[0] if l_[0] >= 0 else l_[1] for l_ in l ]  
        self._inc = self.__edges_to_incidence__(self._edges, l)
        self._internal_edges =  np.where(self._inc[-1]==0)[0]
        self._external_edges =  np.where(self._inc[-1]!=0)[0]
        
    def __gen_ex_vert__(self):
        """Im not sure what the definition of external vertices in our diagrams is yet - it may be either source and sink or those in the residual"""
        for i, v in enumerate(self.ci.tensor):
            if v.sum() > 0:  yield i
            
    def __repr__(self):  return self._repr_latex_()
    
    def _repr_latex_(self):
        init_printing(use_latex='mathjax')    
        return latex(Matrix(self._inc),  mode='inline')
      
    def __edges_to_incidence__(self, edges, non_distinct_external_vertices):
        distinct_vs = np.unique(reduce(add, [list(i) for i in edges],[]))
        n_nodes=len(distinct_vs)
        #this is a standard graph theory thing to create a planar extension connecting to infinity node
        #here i got to me edge dataframe and take everythig in the residual
        #print(edges, non_distinct_external_vertices)
        for ex in non_distinct_external_vertices: 
            infinity_node = distinct_vs.max() +1
            edges.append((ex,infinity_node))
        if len(non_distinct_external_vertices) > 0: n_nodes+=1
        inc = np.zeros((n_nodes,len(edges)),dtype=np.int)
        for i, e in enumerate(edges):
            inc[e[1],i] = 1
            inc[e[0],i] = -1    

        return inc
    
    @property
    def gamma_integral(self):
        propagators = cpropagator.from_edges(self.ci.edges,k_node_index=1)
        K = cpropagator.integrate(propagators)
        PI = kpropagator(K)
        return PI.gamma_integral()
            
    @property
    def incident_matrix(self):return self._inc
    
    @property
    def laplacian(self):return np.dot(self._inc,self._inc.T)
    
    def laplacian_from_adj(self, adj): return adj.sum(axis=1) * np.eye(adj.shape[0],dtype=np.int) + (adj * -1)
    
    @property
    def is_1PI(self):  
        inc = self._inc
        ex_cols = [ c for c in range(inc.shape[1]) if inc[:,c][-1]== 1]
        external_vertices = list(set([ list(np.where(inc[:,c]==-1)[0])[0] for c in ex_cols]))
        internal_valance = [  len(np.where(inc[v][:-len(ex_cols)] != 0)[0]) for v in external_vertices]
        return not np.any(np.array(internal_valance)<2)

    
    def betti_number(self,n=1):
        """Assuming the diagram has one connected component E - V"""
        #for now retur the one we calculated in the construction
        return self.ci.loops

    
    def sym_i_cofactor(self, mat,r=-1,c=-1):
        """remove ith row and ith column"""
        temp = Matrix(mat)
        if r==-1: r = temp.shape[0]-1
        if c==-1: c = temp.shape[1]-1
        #print("removing rows and cols", r,c)
        temp.col_del(c)
        temp.row_del(r)
        return temp


    def symbol_laplacian(self):
        """
        Given the incidence matrix, we have edges (columns) joining vertices (rows)
        Here the edge index determines the symbol while the laplacian is a matrix of size |V|x|V|
        We go through the inc matrix for each vertex, and we get the edge id and the location of the other vertex
        """
        inc = self._inc
        external_edges =  np.where(inc[-1]!=0)[0]
        """This computes the input of Kirchoff polymial (reduced laplacian) in symbol form"""
        def _int_to_dummy(i):  
            if i in external_edges: return Symbol("z")
            return Symbol("x"+str(i+1))
        
        def corresponding_vertices(inc, v):
            """returns a tuple (edge_id, other_vertex_id)"""
            for e in np.nonzero(inc[v])[0]:  yield (e, list(set(np.nonzero(inc[:,e])[0]) - set([v]))[0])

        def _symbol_row_(vec,v):
            _ret = [0 for i in range(inc.shape[0])]
            others = list(corresponding_vertices(inc,v))
            ##all the edges are added to the diagnal where they point to another vertex
            _ret[v] = reduce(add,[_int_to_dummy(t[0]) for t in others],0)
            #and the other vertex location is marked with edge symbol too
            for e,v in others: _ret[v] = -1*_int_to_dummy(1*(e))         
            return _ret

        return Matrix([_symbol_row_(row,v) for v, row in enumerate(inc)])

    def kirchhoff_poly(self):
        lap = self.symbol_laplacian()
        reduced = self.sym_i_cofactor(lap)
        P= reduced.det(method='berkowitz')
        return collect(expand(P),Symbol("z"))

    def _polys_as_edge_index_tuples_(self, Ps):
        def reduction(P):  return [tuple([ int(a.name[1:]) for a in term.atoms() ])  for term in P.args ]
        return (reduction(Ps[0]),reduction(Ps[1]) )

    def graph_polynomials(self, as_tuple_index=False, cut_momenta=None,masses=None):
        S = self.kirchhoff_poly()
        U = self.__W__(S, order=1)
        F0 = self.__W__(S, order=2)

        #this is just for convenience of mapping to other things e.g. plotting forests
        if as_tuple_index: return self._polys_as_edge_index_tuples_( (U,F0) )
        if masses == None:masses = [1 for i in range(len(self._internal_edges))]
        try:
            if cut_momenta == None:cut_momenta = [1 for i in range(len(F0.as_ordered_terms()))]
            F0 = reduce(add,[F0.as_ordered_terms()[i]*m for i, m in enumerate(cut_momenta) ],0)
        except:
            #tis is temporary because i do not know how to handle yet
            pass
        
        mass_sum = reduce(add, [Symbol("x"+str(i+1)) * masses[i]  for i in self._internal_edges], 0)
        
        #multiple F0 terms by the corresponding cut momentum
        
        F = F0 + U*mass_sum
        return U, F0, F

    def __W__(self, P, order = 1):
        """This is a helper function that computes the inverted coefficients for the kirchoff polymial
        """
        ##for the polynomials this is important - the divisor is the number of internal edges we know about
        max_index = len(np.where(self._inc[-1]==0)[0])
        #this is something i need to think about - there is a coefficient of the polynomial actually which i know how to deal with
        def _remove_constants_(s):
            first = s.args[0]
 
            return s if not first.is_constant() else s/first
        
        def _inversion_(S):
            """This is actually a set complement"""
            pre_term = reduce(mul, [Symbol("x"+str(i)) for i in range(1,max_index+1)],1)
         
            return reduce(add , [_remove_constants_(pre_term/t) for t in S.as_ordered_terms() ], 0)

        return  _inversion_(P.coeff(Symbol("z"),order))

    def compute_neg_cut_momenta(self):
        return 1
    
    def should_toggle_edge_columns(arr, axis=0, prefactor=-1):
        mask = arr!=0
        nonzero_index = np.where(mask.any(axis=axis), mask.argmax(axis=axis), 0) #has column is the invalid value
        #prefactor means that if we want the start to be -1, we multiply by -1 to swithc the truth table
        #that way, we can act on the matrix by 1 or -1 depending on if we want to toggle the edges
        return prefactor* np.array([ arr[cs,index] for index, cs in enumerate(nonzero_index)])


    @property
    def edge_system(self):return self.__edge_system__()
        
    def __edge_system__(self):
        inc = self._inc
        #the first off the toggle result - assumed exactly one here. think about this
        k_edge = list(np.where(composite_interaction_graph.should_toggle_edge_columns(inc)==-1)[0])
        edges = composite_interaction_graph.edge_labels(inc, k_edge)
        eq = []
        external_vertex_ids = np.where(inc[-1]==1)[0]  
        #return external_vertex_ids, k_edge

        for index, r in enumerate(inc[:-1]):
            vertex_incs = np.where(r!=0)[0] #get the row entries where there is something going on
            #create a list of terms where we know the symbol label and the sign
            s = [r[i]*edges[i] for i in vertex_incs ]
            res = reduce(add, s, 0)
            eq.append(res)

        eq.append(reduce(add, [ edges[i] for i in external_vertex_ids ], 0))
        return eq, solve(eq), edges


    def internal_edges_index(inc): return [  c for c in range(inc.shape[1]) if inc[-1,c] != 1 ] 

    def shortest_path(graph, start, end, path=[]):
        path = path + [start]
        if start == end:  return path
        if start not in graph:  return None
        shortest = None
        for node in graph[start]:
            if node not in path:
                newpath = composite_interaction_graph.shortest_path(graph, node, end, path)
                if newpath:
                    if not shortest or len(newpath) < len(shortest):
                        shortest = newpath
        return shortest
    
    def expand_path(path):
        for enum, i in enumerate(path[1:]):
            edge = (path[enum], i)
            #find this edge in the incident matrix and invert it (assuming started from canon)
            yield list(edge)

    def edges_to_adjaceny_dict(edges):
        edges = [e["edge"] for e in edges]
        d= {}
        for e in edges:
            s = e[0]
            if s not in d:d[s] = []
            d[s].append(e[1])
        return d

    def __vertex_vertex_edges__(inc):
        for c in range(inc.shape[1]):
            col = inc[:,c]
            col = list(np.nonzero(col)[0])
            yield col

    def invert_path(inc, path):
        vve = composite_interaction_graph.__vertex_vertex_edges__(inc)
        inversion = []
        for index, v in enumerate(vve):
            for p in path:
                multiplier = -1 if np.array_equal(p,v) else 1
                if multiplier == -1:continue
            inversion.append(multiplier)
        return  inversion * inc

    def edge_labels(inc,short_edge=[1]):
        #how to designate k

        externals = np.where(inc[-1]!=0)[0]
        internals = list(set(range(inc.shape[1])) - set(externals) - set(short_edge))
        edges = {}
        for index, i in enumerate(internals): edges[i] = Symbol("q"+str(index+1))
        for index, i in enumerate(externals): edges[i] = Symbol("p"+str(index+1))

        if len(short_edge)> 0 : edges[short_edge[0]] = Symbol("k")
        return edges

    def draw_decomp(self):
        i,r = self.tensors_from_incidence()
        vs = i + r
        return diagram_set([G(J(v)) for v in vs])
        

    def tensors_from_incidence(self, num_species=2):
        IN_FIELD = 0
        OUT_FIELD = 1
        G = self
        temp = G.ci.graph_dataframes()[1]
        num_vertices = G._inc.shape[0] - 1
        species = temp.species.values   
        #a convention is needed- this work for externals, but for internals we need convention - its in field if the vid for the source is less than for the target
        #then later, one has to interpret; check your vids, and check polarity
        field_types = (temp.source<temp.target).astype(int).values 

        #for each edge get the index, species, IO_field type and vertex id
        external_edges = [(c,species[c], field_types[c], list(np.where(G._inc[:,c]==-1)[0])[0]) for c in range(G._inc.shape[1]) if G._inc[:,c][-1]==1]
        residual = np.zeros((num_vertices, num_species,2 ), dtype=np.int)
        for e in external_edges:  residual[e[3]][e[1]][e[2]] +=1

        map_list = lambda x: [IN_FIELD if c == 1 else OUT_FIELD for c in x]
        #for each edge get the index, species, IO_field type and vertex id
        internal_edges = [(c,species[c], map_list(list(G._inc[:,c][G._inc[:,c]!=0])), list(np.where(G._inc[:,c]!=0)[0])) for c in range(G._inc.shape[1]) if G._inc[:,c][-1]==0]

        internals = np.zeros((num_vertices, num_species,2 ), dtype=np.int)
        for e in internal_edges: 
            #print(e)
            #this is a bit weird but I have split the edge in two leaving stubs (flags) on the respective vertices
            #im making a big assumption here and that is that we always go out to the larger vertex so the last index is simple 0,1
            internals[e[3][0]][e[1]][0] +=1
            internals[e[3][1]][e[1]][1] +=1

        #return internal_edges
        return internals, residual



##############################################
def _combine(set_of_interactions):
        def _merge(a,b): 
            if a is None:return b
            if b is None:return a
            return a * b
        return reduce(_merge, set_of_interactions, None)
    
def graph_permutations(primitives, loop_orders = [0,1], max_k=3):
    if not isinstance(loop_orders, list):loop_orders = [loop_orders]
    l = []   
    for i in range(2,max_k+1):
        for tup in list(itertools.permutations(primitives,i)):
            res = _combine(list(tup))
            if res.loops in loop_orders:
                l.append(res)
    #todo : define uniqueness and validity
    return list(set(l))

def distinct_loops(fgset,insist_on_same_residual=False):
    #I think only the same residual makes sense because otherwise the external momenta are different? but then its the same graph so maybe not
    d = {}
    for f in fgset:
        if f.loops < 1 or not composite_interaction_graph(f).is_1PI:continue
        inv_res = f.residual_complement
        inv_res.flags.writeable = False
        h = hash(str(inv_res))
        if h not in d:d[h] = []
        d[h].append(f)
    return d

######################################################
# some interesting things for later
#https://www.math.ucdavis.edu/~amoore/code/kauffmanstates/kauffmanstates.py