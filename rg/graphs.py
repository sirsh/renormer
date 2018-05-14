#the main innovation i can make is use the global representation in the graph theoretic sense for the graphs
#and then have a vetex collection which shows the tensor decompositions in terms of residual and internal half edges
#i can just keep 2 of three of these set types e.g. residual and internal (which can be summed)
#graph.vertex_residuals, graph.vertex_flags

from sympy import latex,init_printing,Matrix,collect,Symbol,expand,solve
import numpy as np
from functools import reduce
from operator import mul,add

import itertools

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
        if f.loops < 1:continue
        inv_res = f.residual_complement
        inv_res.flags.writeable = False
        h = hash(str(inv_res))
        if h not in d:d[h] = []
        d[h].append(f)
    return d

######################################################

class composite_interaction_graph(object):
    def __init__(self, comp_interaction,connect_inf=False):
        self.ci = comp_interaction
        temp=self.ci.graph_dataframes()[1]
        self._edges = [d["edge"] for d in self.ci.edges]
        #todo - get the external vertices which will be connected to the infinity node
        self._externals = [] if connect_inf == False else list(self.__gen_ex_vert__())
        #vids are the non 0 source or tarket in the external edges dataframe
        l = list(zip( temp[temp.internal==False].source, temp[temp.internal==False].target))
        l = [ l_[0] if l_[0] >= 0 else l_[1] for l_ in l ]  
        self._inc = self.__edges_to_incidence__(self._edges, l)
     
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
    def incident_matrix(self):return self._inc
    
    @property
    def laplacian(self):return np.dot(self._inc,self._inc.T)
    
    def laplacian_from_adj(self, adj): return adj.sum(axis=1) * np.eye(adj.shape[0],dtype=np.int) + (adj * -1)
    
    @property
    def is_1PI(self):  return len(np.where(np.diagonal(self.laplacian)<2)[0]) == 0
    
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
        external_start_index=-1
        edges = self._edges
        inc = self._inc
        external_edges = self._externals
        if len(external_edges): external_start_index = len(edges) +1 - len(external_edges)
        #print("external edges start at", external_start_index)
        """This computes the input of Kirchoff polymial (reduced laplacian) in symbol form"""
        def _int_to_dummy(i):  
            if i >= external_start_index and external_start_index > 0: return Symbol("z")# +str(i)for testing we might care about which z
            return Symbol("x"+str(i))
        def _symbol_row_(vec,v):
            _ret = [0 for i in range(inc.shape[0])]
            entries = list(np.where(vec!=0)[0])
            """take the incident matrix and use the zero indexed enries to create dummy variables"""
            _ret[v] = reduce(add,[_int_to_dummy(i+1) for i in entries],0)
            #foreach entries, look at the inc column and find my adjoint
            for e in entries:
                adj = list(set(np.nonzero(inc[:,e])[0]) - set([v]))[0]
                _ret[adj] = -1*_int_to_dummy(1*(e+1))         
            return _ret

        return Matrix([_symbol_row_(row,v) for v, row in enumerate(inc)])

    def kirchhoff_poly(self):
        inc = self._inc
        edges = self._edges, 
        infinity_connect_list=self._externals
        external_edges = edges[-len(infinity_connect_list):] if len(infinity_connect_list) > 0 else []
        lap = self.symbol_laplacian()
        reduced = self.sym_i_cofactor(lap)
        P= reduced.det(method='berkowitz')
        return collect(expand(P),Symbol("z"))

    def graph_polynomials(self, include_masses=False):
        max_ind = len(self._edges)
        S = self.kirchhoff_poly()
        U = self.__W__(S, order=1, max_index=max_ind)
        F0 = self.__W__(S, order=2, max_index=max_ind)
        NCM = self.compute_neg_cut_momenta()
        SUM_MASS = 0 if include_masses == False else reduce(add, [Symbol("x"+str(i)) *Symbol("m"+str(i))**2 for i in range(1,max_index+1)],1)
        F = F0
        return U, NCM * F + U * SUM_MASS

    def __W__(self, P, max_index, order = 1):
        """This is a helper function that computes the inverted coefficients for the kirchoff polymial
        """
        def _inversion_(S):
            """This is actually a set complement"""
            pre_term = reduce(mul, [Symbol("x"+str(i)) for i in range(1,max_index+1)],1)
            return reduce(add , [((1/t)*pre_term).as_content_primitive()[-1] for t in S.args ], 0)

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



