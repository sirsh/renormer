#TODO
#1. add colour (integrate with drawing), residual operators, loop hash,
# RG can have global species colours, species labels, etc.
#having dones this it is possible to create a propagator matrix format to construct propagators the vector way

#implement edge sets by wrapping set and overloading /set notation

#temp: later pick your objects - just is just to prototype
from sympy import *

import operator
import functools
import math
import pandas as pd
import numpy as np

from .. import _sum_, _product_


EXITING = -1
ENTERING = 1
#todo: many objects that could be preloaded once that arre not such as the dataframe adjacency with metadata (filterable by trees and cicuit edge ids)

class _delta_(object):
    def __init__(self, name, ids, signs):
        self._name = name
        self._d = Symbol("delta_"+str(ids[0]))
        self._vec = symbols(" ".join([name+"_"+str(abs(v)) for v in ids]))
        signed = lambda v : [int(signs[i]) * Symbol(name+"_"+str(v)) for i,v in enumerate(ids)]
        self._signed = signed(self._vec)
        self._exp = functools.reduce(operator.add, self._signed, 0)
        self._rep = Function("delta_"+str(ids[0]))(self._exp)
    def __repr__(self):  return self._repr_latex_()
    def _repr_latex_(self): return latex(self._rep,  mode='inline')
    def __call__(self, i):
        return solve(self._exp, Symbol(self._name+"_"+(str(i))))[0]
    
    @property
    def value(self):return self._rep
        
    def apply(self,exp):
        term = self._signed[0]
        s = solve(self._exp, term)[0]
        return exp.subs({term : s})
    
class internal_edge_set(object):
    def __init__(self, inc):
        self.inc = inc
        self.T = list(inc.spanning_trees.iloc[0].values)
        self.Tc = inc.edge_complement(self.T)
        self.G = list(range(len(inc.edges)))   
        
    @property 
    def I_S(self):
        mat = self.inc.copy()
        mat[:,self.Tc] = 0
        return mat[:-1,self.G]
    
    @property 
    def IC_S(self):
        mat = self.inc.copy()
        mat[:,self.T] = 0
        return mat[:-1,self.G]
    
    @property 
    def II_S(self):
        mat = self.inc.copy()
        return mat[:-1,self.T]
    
    @property 
    def IIC_S(self):
        mat = self.inc.copy()
        return mat[:-1,self.Tc]
    
    @property
    def X(self):  return np.round(-1*np.matmul( np.linalg.pinv(self.I_S),(self.IC_S) )).astype(int)
    
    @property
    def X_delta(self): 
        for i, col in enumerate(self.X.T):
            m = list(np.nonzero(col)[0])
            if len(m) > 0: yield _delta_("alpha", [i]+m, [1]+list(col[m]))
    
class incidence_matrix(np.ndarray):  
    #things to watch out for - handling multiple edges
    #todo allow a dictionary of letters maybe for vertices to make it easier to init and maybe consider not using sizes
    #alternative loading model where we do not compute the k_matrix?
    #alternative where we do not add the edge count and vertex count - seems a bit inconvenient
    def __new__(self, edge_count=None, vertex_count=None, edges=[],external_vertices=[], species_vector=None):
        #parse
        internal_edges = [e for e in edges if  -1 not in e]
        external_edges =  [e for e in edges if -1 in e]
        if len(external_edges) ==0 and len (external_vertices) > 0: 
            external_edges += [[-1,v] for v in external_vertices]
            edges = internal_edges + external_edges
        vertices = [v for v in np.unique(np.array(edges).flatten()) if v != -1]
        external_vertices = [v for v in np.unique(np.array(external_edges).flatten()) if v != -1]
        edge_count = len(edges)
        vertex_count = len(vertices)    
        size = (len(vertices)+1,len(edges))
        
        #build
        obj = super(incidence_matrix, self).__new__(self, size, np.int)  
        obj[:] = 0 
        obj._num_internal_edges = len(internal_edges)
        
        obj.add_edges(edges)
        
        obj._external_vertices = external_vertices
        obj._internal_vertices = list(set(list(range(vertex_count))) - set(external_vertices))
        
        obj._kmatrix = obj.kirchhoff
        obj._betti_1 = obj._kmatrix.degree()#~ whats the relatioship?

        obj.species_vector = species_vector
        
        return obj
    
    #add external edge and internal edges sets for convenience
    
    @property #the last row determines if the edge is connected to vinty i.e. external
    def num_internal_edges(self):len(np.where(self[-1] == 0)[0])
        
    #for testing it would be good to return this in incidence matrix format 
    def get_star(self,vid):
        star_links = {}
        for i,e in enumerate(self.T):
            if e[-1] == 0: star_links[i] = list(e)
        star_links
        m = np.zeros(self.shape,np.int)
        for k in star_links: m[:,k] = star_links[k]

        my_edges = self.get_incident_edges(vid)
        targets = []
        vs = self.shape[0] -1
        for v in self[:,my_edges].T:
            targets.append([ c if c < vs else -1 for c in np.nonzero(v)[0] if c != vid][0])

        return np.array(list(zip(my_edges,list(self[:,my_edges][vid]),targets))) 

    
    def residual_matrix(self):
        redges = []
        species = []
        for i,e in enumerate(self.T):
            if e[-1] == 1: redges.append([0,-1])
            if e[-1] == -1: redges.append([-1,0])
            if self.species_vector != None:
                species.append(self.species_vector[i])
        return incidence_matrix(edges=redges, species_vector=species if len(species)>0 else None)
        
        
    def residual(self, exiting=True, entering=True, sort=False, as_tensor=False):
        
        def _residual_tensor_(r,n=2):
            """
              Create a tensor for this direction convention <----
                              [ [ species0_out, species0_in],
                                [ species1_out, species1_in],
                                ...
                              ]
             #for now i harcode the maximum species which may not be known to graph but for my scope 2 is enough for now
            """
            species = np.unique(r[:,-1])
            n = np.maximum(n, len(species))
            t = np.zeros((n,2),np.int)
            for row in r:
                index = 1 if row[1] == 1 else 0
                t[row[-1]][index] += 1
            return t

        external_edges = self[-1][self.num_internal_edges:]
        residual_info = np.array(list( zip(  range(len(external_edges)), external_edges*-1, self.species_vector if self.species_vector is not None else [0 for i in range(len(external_edges))])))
        if sort: residual_info=residual_info[residual_info[:,-1].argsort()]
        if exiting == False:residual_info=residual_info[residual_info[:,1] != -1]
        if entering == False:residual_info=residual_info[residual_info[:,1] !=  1]
        return residual_info if not as_tensor else _residual_tensor_(residual_info)


    def vedges(self, voffset=0):
        """
        This coverts all edges into [v,w] and it is smart enough to treat v_infty
        Also, if we want to choose an ID offset for vertices when naming them, it is convenient to do so here
        """
        vinfty = self.shape[0] -1

        for e in self.T:
            f,t = np.where(e==EXITING)[0][0],np.where(e==ENTERING)[0][0]
            if f == vinfty:f=-1
            else: f = f+voffset
            if t == vinfty:t=-1   
            else: t = t+ voffset
            yield [f,t]

    @property
    def _internal_edge_set_(self):return internal_edge_set(self)
        
    def ensure_edge_directions(self,d):
        """Supply a dictoinary with edge id and vertex e.g. { 0: [1,2]} means edge 0 takes vertex 1 to vertex 2"""
        for k,v in d.items():self[:,k][v[0]],self[:,k][v[1]] = EXITING,ENTERING
        
    def flip_edges(self,edge_indecies):
        if not isinstance(edge_indecies,list): edge_indecies =[edge_indecies]
        for e in edge_indecies:self[:,e]*= -1
        return self
        
    def reset_edge_orders(self, edges):pass
    def reset_edge_orders_from_hierarchy(self):pass
    
    def get_all_connections(self,a,b):
        for i,r in enumerate(self.T):
            if(r[a]==ENTERING) and (r[b]==EXITING):yield i, ENTERING
            if(r[a]==EXITING) and (r[b]==ENTERING):yield i, EXITING
            
    def try_directed_edge_index(self,a,b):
        """
        I think I should reverse this - postive end should be arriving
        checks all the columns (edges) to see if there is an edge that starts at a(-1) and ends at b(1)
        """
        for i,r in enumerate(self.T):
            if(r[a]==EXITING) and (r[b]==ENTERING):
                return i
        return None
    
    def exists_directed_edge(self,a,b):
        return self.try_directed_edge_index != None

    #dont override these on a class that overrides ndarray - crazy!
    #def __mul__(self, B):
    #    return self.dprod(B)
    #def __rmul__(self, B):
    #    return self.dprod(B)
    
    def dprod(self, B, species_chord_set):
        #get my ins
        #todo figure out a maximal chord set if nothing passed in
        inputs = self.residual(exiting=False)
        outputs = B.residual(entering=False)
        ind,outd = {},{}
        for i in inputs:
            if i[-1] not in ind: ind[i[-1]] = [] #check species type
            ind[i[-1]].append(i[0]) #add eid to correct species

        for i in outputs:
            if i[-1] not in outd: outd[i[-1]] = [] #check species type
            outd[i[-1]].append(i[0]) #add eid to correct species

        #return ind,outd
        bs = []
        for i,bonds in enumerate(species_chord_set): #this is a tuple e.g. (2,0,1) for species A,B,C in the theory
            for b in range(bonds):
                bs.append([ind[i].pop(), outd[i].pop()])

        all_my_edges = list(self.vedges())
        num_my_vertices =  self.shape[0]-1
        all_their_edges= list(B.vedges(voffset=num_my_vertices))

        #then we must determine what target is on our edge and what target is on their edge and make new

        new_edges = []
        for b in bs:
            re = all_my_edges.pop(b[0])
            le = all_their_edges.pop(b[1])
            new_edges.append([le[0], re[-1]])

        new_edges = all_my_edges + all_their_edges + new_edges
        return incidence_matrix(edges=new_edges)

    @property
    def cycle_basis(self): return cycle_finder(self.edges).cycle_basis
    
    #this is a bit problematic, I have not decided on the best conv. unclude or exclude external by default?

    @property
    def edges(self): return [self.directed_edge_as_vertices(e) for e in range(self._num_internal_edges)]
        
    @property
    def internal_vertices (self):return self._internal_vertices
    
    @property
    def external_vertices (self):return self._external_vertices
    
    def get_coincident_vertices_any(self, vset,exclusions = []):
        items = []
        for v in vset:
            for _v in self.get_coincident_vertices(v,exclusions):
                if _v not in items: items.append(_v)
        return items
        
    def get_incident_edges(self, v,internal=False):
        all_edges = np.nonzero(self[v])[0]
        if internal:        
            return [e for e in all_edges if e < self._num_internal_edges]
        return all_edges
    
    def get_coincident_vertices(self, v,exclusions = []):
        #these are the data columns for which v has edge data
        v_inc_frame = self[:,np.where(self[v]!=0)[0]]
        exclusions = exclusions + [v]
        #now get the row index for things that are co-incident using the trick but just checking non zero rows
        #this is a bit odd and probably better way but count non zero columns uses the count and nonzero combined 
        #and we ensure that v is not in the list which is requred when it has only one edge for example
        #I check if it is incoming or outgoing - we could change it to in and out by changing != 0 to +1
        #the 1 is important for directionality bias e.g. when drawing we only want to draw the directed edge or handle unidirectionality
        incident_vertices = [l for l in list(np.nonzero((np.count_nonzero(v_inc_frame,axis=1)>0).astype(int))[0]) if l not in exclusions]
        return incident_vertices
        
    def add_edges(self, eset):
        for i, v in enumerate(eset): 
            self[v[0],i],self[v[1],i] = EXITING,ENTERING
            
            
    @property
    def tree_distance_map(self):
        dmap = {}
        for h in self.vertex_hierarchy():
            for l in h[1]:dmap[l]= h[0]
        return dmap

    def vertex_hierarchy(self, vset=[-1], d=0, h={} ):
        """
        recursion from external nodes or other seed set
        """
        if d ==0: h= {}
        #need to re-think the indexing a little
        #this is so that I can store the last vertex (inf) as either -1 or the length of the thing
        if -1 in vset: h[len(self)-1] = 2
        for _v in vset: h[_v] = d
        if len(vset) > 0:yield d, vset
        #look for any coincident vertices that are not on my level or above
        next_layer = self.get_coincident_vertices_any(vset ,vset+ list(h.keys()))
        if(len(next_layer)>0):
            for t in self.vertex_hierarchy(next_layer,d+1,h):
                yield t
    
   
    def __coincident_vertex__(self, v, e):
        col = self[:,e]
        [e for e in np.where(self[:,1]!=0)[0] if e != 2][0]
        return [i for i in np.where(col != 0)[0] if i != v][0]
   
    def get_alt_adjacency(self):
        d = {}
        for i,r in enumerate(self):
            eset = self.get_incident_edges(i)
            l = []
            for e in eset: l.append(( self.__coincident_vertex__(i,e),e, "I" if e < len(self.edges) else "E")) 
            d[i] = l
        return d
    
    @property
    def flat_adjacency(self): return [ list(np.nonzero(v)[0]) for v in  self.adjacency_with_inf]
      
    @property
    def adjacency(self):
        #[ np.nonzero(v)[0]for v in  INC.adjacency_with_inf]
        return self.adjacency_with_inf()[:-1,:-1]
  
    def adjacency_with_inf(self,restricted_edges=None,use_edge_labels_1indexed=False):
        ad = np.zeros((len(self),len(self)),np.int)
        for i, v in enumerate(self):  
            if restricted_edges is None or i in restricted_edges:
                #does not work for multiple edges - what do we want to do about this?
                #todo - count how many times it is connected and then create a lookup of edges
                ad[i][self.get_coincident_vertices(i)]  =1 if not use_edge_labels_1indexed else self.get_incident_edges(i)+ 1
        return ad#[:-1,:-1] if len(self._external_vertices) > 0 else ad
    
    @property
    def edge_labels(self):
        #the last row is the external vertex so check which edges are connected
        ex_edges = np.nonzero(self[-1])[0]
        return [ "p"+str(i) for i in range(self.shape[-1]) if i not in ex_edges]   +   [ "q"+str(i) for i in range(len(ex_edges))]
       
    @property
    def vertex_constraints(self,display=False):
        edge_labels = self.edge_labels
        for v in self:
            yield delta_constraint([edge_labels[i] for i in np.nonzero(v)[0]])
            
    @property
    def edge_label_matrix(self): return Matrix(self) * diag( *self.edge_labels)
    
    def constraint_system(self,solve_it=False): 
        equations = []
        incmat = self.edge_label_matrix
        for i in range(incmat.shape[0]): 
            equations.append(functools.reduce(operator.add, incmat.row(i), 0))
        mat= Matrix(equations)
        return mat if not solve_it else solve(mat)
        
    @property
    def symbol_adjacency(self):pass
    
    def __sym_i_cofactor__(mat,r=-1,c=-1):
        """remove ith row and ith column"""
        temp = Matrix(mat)
        if r==-1: r = temp.shape[0]-1
        if c==-1: c = temp.shape[1]-1
        temp.col_del(c)
        temp.row_del(r)
        return temp
    
    def _internal_symbols_(self): 
        ex_edges = np.nonzero(self[-1])[0]
        return [Symbol("p"+str(i)) for i in range(self.shape[-1]) if i not in ex_edges]
    
    def _external_symbols_(self): 
        ex_edges = np.nonzero(self[-1])[0]
        return [Symbol("q"+str(i)) for i in range(len(ex_edges))]
    
    @property
    def kirchhoff(self):
        #get the laplacian
        L = self.edge_label_matrix * self.edge_label_matrix.T
        #scale the terms because of the way we did it
        for s in self._internal_symbols_(): L=L.subs({s**2:s})
        #treat all external edges as z so we can get a polynomial in z for forest order
        for s in self._external_symbols_(): L=L.subs({s**2:Symbol("z")}) 
        #take the reduced laplacian
        reduced =  incidence_matrix.__sym_i_cofactor__(L)  
        #solve determinant
        P= reduced.det(method='berkowitz')
        #make it pretty
        res = collect(expand(P),Symbol("z"))
        res = Poly(res, Symbol("z"))
        return res
    
    @property
    def first_betti_number(self): return self._betti_1
    
    @property
    def spanning_trees(self): return pd.DataFrame(self.k_forests(-1))
    
    def k_forests(self, order): 
        p = self._kmatrix.coeffs()[order]
        #this must be considered temporary - one thing might be to use symbols that evaluate to useful things instead of bare symbols
        for monomials in p.as_ordered_terms(): 
            edges = [int(str(e).replace("p","")) for e in monomials.free_symbols]
            yield edges
            
    def edge_complement(self,edges):
        return list(set(range(len(self.edges))) - set(edges))
    
    def directed_edge_as_vertices(self,e):
        """
        the edge goes from cell marked 1 to cell marked -1
        """
        e = self[:,e]
        return [np.where(e==EXITING)[0][0],np.where(e==ENTERING)[0][0]]
    
        
    def __sample__():
        return incidence_matrix(9, 6,
                       [ [0,1], [1,2], [2,0], [3,4],[4,3],[3,0],[4,2], [4,5], [5,3] ],
                       external_vertices=[0,1,2])
        
    
    #this is a work in progress
    def edges_to_loop_basis(self):
        from . import circuits,set_source_sink_flow
        """
        each edge is described as a linear combination of loop and external momenta 
        but exactly how to do this properly is not clear to me - it basically comes down to proper signs but the signs are important
        presume that if we take any circuit we can find something that splits omega signs on the loop or the integral would vanish
        (nice to have) would be to ne consistent with how we choose sink source vertices and which edge we choose as the loop carrier
        """
        inc_mat = self
        def _split_(p): return str(p)[:1], int(str(p)[1:])  

        #this may be important and also it may be important to control which edge is the loop edge but not sure yet
        for c in circuits(inc_mat): set_source_sink_flow(inc_mat,c)    
        sol = inc_mat.constraint_system(True)

        missing= set(inc_mat._internal_symbols_()) - set(sol.keys())
        loop_edges = [_split_(p)[-1] for p in missing]
        L = len(loop_edges)
        N = L + len(inc_mat._external_symbols_())
        edge_basis_map = {}

        def _to_edge_basis_(s,loop_edges):
            def _sign_of_(p):return s.as_coefficients_dict()[p]      
            b = np.zeros(N,np.int)
            for l in loop_edges:b[loop_edges.index(l)] = _sign_of_(Symbol("p"+str(l)))
            for i in range(N-len(loop_edges)): b[i+len(loop_edges)] = _sign_of_(Symbol("q"+str(i)))
            return b

        for p in sol:
            t,e = _split_(p)
            if t == "p":edge_basis_map[e] = _to_edge_basis_(sol[p], loop_edges)      

        for e in loop_edges:
            b = np.zeros(N,np.int)
            b[loop_edges.index(e)] = 1
            edge_basis_map[e] = b

        return edge_basis_map

    
    #TODO check for vanishing momentum on 2-forests
    def _symanzik_pols_(self, symbolic=False, reduce=False, apply_deltas=False,use_zero_mass=False):
    
        #todo finish second polynomaisl which is U*mass_sum + sum of forests multiplied by the momentum across the cut edges
        delta_functions = internal_edge_set(self).X_delta
        
        first = [self.edge_complement(k) for k in self.k_forests(order=-1)]
        #needs to have at least three vertices for this to make sense - otherwise no way to disconect the vertices and still keep an edge
        second = [self.edge_complement(k) for k in self.k_forests(order=-2)]
        
        to_alpha_symbols = lambda l : [Symbol("alpha_"+str(i)) for i in l]
        
        if symbolic:
            a,b= [to_alpha_symbols(t) for t in first],[to_alpha_symbols(t) for t in second] 
            massum = _sum_([Symbol("alpha_"+str(i))*Symbol("m_"+str(i)) for i in range(len(self.edges))]) if not use_zero_mass else 0 
            if reduce: 
                a,b= [_product_(l) for l in a], [_product_(l) * Symbol("\hat{p}") for l in b]
                a,b = _sum_(a),_sum_(b) 

                b = a*massum + b
                if apply_deltas:
                    for d in delta_functions: a,b = d.apply(a),d.apply(b)  

                return a,b
            return a,b
        return first,second 

