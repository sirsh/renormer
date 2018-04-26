from sympy import latex,init_printing,Matrix,collect,Symbol,expand
import numpy as np
from functools import reduce
from operator import mul,add
class composite_interaction_graph(object):
    def __init__(self, comp_interaction,connect_inf=False):
        self.ci = comp_interaction
        self._edges = [d["edge"] for d in self.ci.edges]
        #todo - get the external vertices which will be connected to the infinity node
        self._externals = [] if connect_inf == False else list(self.__gen_ex_vert__())
        self._inc = self.__edges_to_incidence__(self._edges,self._externals)
     
    def __gen_ex_vert__(self):
        """Im not sure what the definition of external vertices in our diagrams is yet - it may be either source and sink or those in the residual"""
        for i, v in enumerate(self.ci.tensor):
            if v.sum() > 0:  yield i
            
    def __repr__(self):  return self._repr_latex_()
    
    def _repr_latex_(self):
        init_printing(use_latex='mathjax')    
        return latex(Matrix(self._inc),  mode='inline')
      
    def __edges_to_incidence__(self, edges,connect_to_inf=[]):
        distinct_vs = np.unique(reduce(add, [list(i) for i in edges],[]))
        n_nodes=len(distinct_vs)
        for ex in connect_to_inf: #this is a standard graph theory thing to create a planar extension connecting to infinity node
            infinity_node = distinct_vs.max() +1
            edges.append((ex,infinity_node))
        if len(connect_to_inf) > 0: n_nodes+=1
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

        
        