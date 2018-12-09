#TODO
#1. add colour (integrate with drawing), residual operators, loop hash,
# RG can have global species colours, species labels, etc.
#having dones this it is possible to create a propagator matrix format to construct propagators the vector way

#implement edge sets by wrapping set and overloading /set notation

#temp: later pick your objects - just is just to prototype

#list of cruicial operations: covertex is the hard one (e, dir, colour)
from sympy import *

import operator
import functools
import itertools
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
        obj._species_vector = species_vector if species_vector != None else [0 for i in range(edge_count)]
        obj[:] = 0 
        obj._num_internal_edges = len(internal_edges)     
        obj.add_edges(edges)
        obj._external_vertices = external_vertices
        obj._internal_vertices = list(set(list(range(vertex_count))) - set(external_vertices))
        obj._kmatrix = obj.kirchhoff
        obj._betti_1 = 0  
        obj.is_tree = obj._num_internal_edges == vertex_count - 1
        if len(obj.spanning_trees)  > 0 and not obj.is_tree:
                #order = obj._kmatrix.degree() - 1
                p = obj.kirchhoff.coeff_monomial(Symbol("z"))
                #K gives the edges IN the tree and we want the number of cuts which gives the loop number
                obj._betti_1 = obj._num_internal_edges - len(p.as_ordered_terms()[0].free_symbols)
       
        return obj
     
    #add external edge and internal edges sets for convenience
  
    #def __eq__(self, b):return self.graph_hash() == b.graph_hash()
    #def __hash__(self):  return hash(self.graph_hash() ) 
    
    @property
    def species_vector(self):
        try: return self._species_vector
        except: return False
    @property
    def  bridge_count(self):
        """
        need to check if a proper defintion or there is not something cleaner, but if we find edges in all spanning trees
        then these are bridges
        """
        trees = self.spanning_trees
        mask = np.zeros((len(trees),self.shape[-1]),np.int)
        for i, r in enumerate(trees.as_matrix()):  mask[i][r] = 1
        mask_tot = mask.sum(axis=0)
        return len(np.where(mask_tot == len(trees))[0])

    @property
    def num_ex_edges(self): return len(np.where(self[-1]!=0)[0] )
    
    @property #the last row determines if the edge is connected to vinty i.e. external
    def num_internal_edges(self):len(np.where(self[-1] == 0)[0])
    
    def graph_hash(self, to_string=True):
        spec = np.array(self.species_vector) + 1
        iedges = self
        def pred(row):return len(np.where(row == 1)[0]) + len(np.where(row == -1)[0])*100
        sel = np.apply_along_axis(pred, 1, iedges)
        h = np.matmul(iedges[np.argsort(sel)], np.diag(spec))
        return h if not to_string else str(h.flatten())

    
    def loop_hash(self, to_string=True):
        eidx = np.where(self[-1]==0)[0]
        spec = np.array(self.species_vector)[eidx] + 1
        iedges = self.T[eidx].T[:-1]
        def pred(row):return len(np.where(row == 1)[0]) + len(np.where(row == -1)[0])*100
        sel = np.apply_along_axis(pred, 1, iedges)
        h = np.matmul(iedges[np.argsort(sel)], np.diag(spec))
        return h if not to_string else str(h.flatten())

    def __mul_hash__(self,use_full_vertex=True,ignore_residue=False):
        vals = []
        
        flat = lambda l : tuple(list(l.sum(axis=0).flatten()))
        
        all_stars = list(self.stars)# [s if ignore_residue else s[np.where(s[:,-1]!=-1)[0]] for s in ]
        for i,r in enumerate(all_stars):
            for j, v in enumerate(r[:-1]):   
                IN = v[:,-1]
                if np.sum(IN) > 0: 
                    #the decodated residual in the species IO basis is the key
                    if use_full_vertex: k = ( flat(all_stars[i]) , flat(all_stars[j]),tuple(IN))
                    else: k = (i,j,tuple(IN) )
                    if k not in vals: vals.append( k )
        #if the set is empty, just return myself and the void mull
        if len(vals) == 0: return set([(flat(all_stars[0]))])
        return set(vals)
    
    @property
    def stars(self):
        """v1: This constructs all information for vertices 
               The external vertex is labelled dim max, then we construct a tensor that contains IO for each species, and a slab for each adjacent vertex
               This allows a summation in z to give the full vertex information but also we can reconsturct the muls
        """
        species = self.species_vector if self.species_vector is not None else list([0 for i in range(self.shape[-1])])
        sn = len(np.unique(species))  
        num_vertices = self.shape[0]
        for i, v in enumerate(self[:-1]):
            tense = np.zeros((num_vertices,sn,2),np.int)
            for e,edge in enumerate(self.T):
                if edge[i] == 0: continue
                IO = 0 if edge[i] == -1 else 1
                vs = np.where(edge!=0)[0] #some edge data
                cov = vs[vs!= i][0] #the vertex that is not me
                tense[cov][species[e]][IO] += 1     
            yield tense
      
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

        external_edges = self[-1][np.nonzero(self[-1])[0]]
        residual_info = np.array(list( zip(  range(len(external_edges)), external_edges*-1, self.species_vector if self.species_vector is not None else [0 for i in range(len(external_edges))])))
        if sort: residual_info=residual_info[residual_info[:,-1].argsort()]
        if exiting == False:residual_info=residual_info[residual_info[:,1] != -1]
        if entering == False:residual_info=residual_info[residual_info[:,1] !=  1]
        return residual_info if not as_tensor else _residual_tensor_(residual_info)

    @property
    def loop_part(self):
        edges = np.array([self.directed_edge_as_vertices(e) for e in range(self.shape[-1])])
        iedges = np.where(self[-1]==0)[0] 
        spec = list(np.array(self.species_vector)[iedges])
        return incidence_matrix(edges=edges[iedges], species_vector=spec)

    @property
    def residue_part(self):
        vinf = self.shape[0] - 1
        def map_inf(e):
            #because this is a residue there can only be two vertices
            return [
                0 if e[0] != vinf else -1,
                0 if e[1] != vinf else -1,
            ]   
        edges = np.array([map_inf(self.directed_edge_as_vertices(e)) for e in range(self.shape[-1])])
        iedges = np.where(self[-1]!=0)[0] 
        #remove internal lvertices
        spec = list(np.array(self.species_vector)[iedges])
        return incidence_matrix(edges=edges[iedges], species_vector=spec)

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

    ###useful for trees where there is a distinguished root vertex
    def root_vertex(self):
        #incomding edges
        v_infty = self.shape[0] - 1
        def _check_covertex_for_in_edges_(v):
            
            edges = np.where(self[v]==1)[0]
            #all covertices of incoming edges
            return [v for v in np.unique(np.where(self.T[edges].T == -1)[0]) if v != v_infty] 
        for i, v in enumerate(self):
            l = _check_covertex_for_in_edges_(i)
            if len(l) == 0 and i != v_infty: return i
        
    #dont override these on a class that overrides ndarray - crazy!
    #def __mul__(self, B):
    #    return self.dprod(B)
    #def __rmul__(self, B):
    #    return self.dprod(B)
    
    #I may implement a graft shuffle product but as a generator it is not well defined 
    #this is because there are inaccessible graphs e.g. those that require trasmutation
    #if we can not rely on it as a complete product then it is only a nice-to-have
    
    @staticmethod
    def shuffle(A,B,allow_split=False):
        """
        return tuples of size species count (n,m,..) for a mul: how many chords to join
        """
        def _complete_marriage_tensor_(A,B):
            rt = A.residual(as_tensor=True)
            lt = B.residual(as_tensor=True)
            return np.minimum(rt[:,-1], lt[:,0])

        sets = [list(range(s+1)) for s in _complete_marriage_tensor_(A,B)]

        return list(itertools.product(*sets))

    
    def dprod(self, B, species_chord_set,outd=None, contract_tree_level=False):
        #get my ins
        #todo figure out a maximal chord set if nothing passed in
        inputs = self.residual(exiting=False)
        outputs = B.residual(entering=False)
        
        ind= {}
        for i in inputs:
            if i[-1] not in ind: ind[i[-1]] = [] #check species type
            ind[i[-1]].append(i[0]) #add eid to correct species

        #this allows us to choose specific edges for the chord set of the left hand side and pass it in. 
        #this is useful for example to combine tree level radicals in specific ways
        if outd == None:
            outd = {}
            for i in outputs:
                if i[-1] not in outd: outd[i[-1]] = [] #check species type
                outd[i[-1]].append(i[0]) #add eid to correct species

        #return ind,outd
        bs = []
        for i,bonds in enumerate(species_chord_set): #this is a tuple e.g. (2,0,1) for species A,B,C in the theory
            for b in range(bonds): bs.append([ind[i].pop(), outd[i].pop()])
        #print(bs)
        
        all_my_edges = list(self.vedges())
        num_my_vertices =  self.shape[0]-1
        all_their_edges= list(B.vedges(voffset=num_my_vertices))
        my_edge_species = self.species_vector
        their_edge_species = B.species_vector
        
        #print(their_edge_species,"before")
        new_edges = []
        new_edges_species = []
        for b in bs:
            le = all_my_edges[b[0]]
            re = all_their_edges[b[1]]
            new_edges.append([re[0], le[-1]]) #from right to left
            new_edges_species.append(my_edge_species[b[0]])
        
        #filter by what we have merged
        #print(all_their_edges)
        
        all_my_edges = [all_my_edges[e] for e in range(len(all_my_edges)) if e not in list(np.array(bs)[:,0])]
        my_edge_species = [my_edge_species[e] for e in range(len(my_edge_species)) if e not in list(np.array(bs)[:,0])]
        
        all_their_edges = [all_their_edges[e] for e in range(len(all_their_edges)) if e not in list(np.array(bs)[:,1])]
        their_edge_species = [their_edge_species[e] for e in range(len(their_edge_species)) if e not in list(np.array(bs)[:,1])]
        #print(all_their_edges)
        
        new_edges = all_my_edges + all_their_edges + new_edges
        new_edges_species = my_edge_species + their_edge_species + new_edges_species
  
        #print(their_edge_species,"after")
        #print("check speciation lengths", len(new_edges),len(new_edges_species))
        #print("check speciation lengths", len(all_my_edges),len(my_edge_species))
        #print("check speciation lengths", len(all_their_edges),len(their_edge_species))  
        res = incidence_matrix(edges=new_edges,species_vector=new_edges_species)
        
        #tree level we always contract the internal edges - only loops hold on to this information
        return res.residue_part if res.is_tree and contract_tree_level else res

    @property
    def cycle_basis(self): return cycle_finder(self.edges).cycle_basis
    
    #this is a bit problematic, I have not decided on the best conv. unclude or exclude external by default?

    @property
    def edges(self): return [self.directed_edge_as_vertices(e) for e in range(self._num_internal_edges)]
        
    @property
    def internal_vertices (self):return self._internal_vertices
    
    @property
    def external_vertices (self):return self._external_vertices
    
    @property
    def internal_edges(self):return list(np.where(self[-1])[0])
        
    
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
        ex_edge_labels = [ "q"+str(i) for i in range(len(ex_edges))]
        labels = np.array(  [ "p"+str(i) for i in range(self.shape[-1]) ]   )
        labels[ex_edges] =  ex_edge_labels
        return list(labels)
       
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
    
    def __sym_i_cofactor__(self, mat,r=-1,c=-1):
        """remove ith row and ith column"""
        temp = Matrix(mat)
    
        #internal_edges = self.internal_edges
        #return temp[internal_edges,internal_edges]
    
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
    def reduced_lap(self):
        #get the laplacian
        el = self.edge_label_matrix
        
        L = el * el.T
        #scale the terms because of the way we did it
        for s in self._internal_symbols_(): L=L.subs({s**2:s})
        #treat all external edges as z so we can get a polynomial in z for forest order
        for s in self._external_symbols_(): L=L.subs({s**2:Symbol("z")}) 
        #take the reduced laplacian
        return  self.__sym_i_cofactor__(L)  
    
    @property
    def kirchhoff(self):
        #if we only have one vertex then there is nothing to do
        if self.shape[0] == 2: return  Poly(0, Symbol("z"))

        reduced = self.reduced_lap
        #solve determinant
        P= reduced.det(method='berkowitz')
        #make it pretty
        res = collect(expand(P),Symbol("z"))
        res = Poly(res, Symbol("z"))
        return res
    
    @property
    def first_betti_number(self): return self._betti_1
    
    @property
    def spanning_trees(self): return pd.DataFrame(self.k_forests(1))
    
    def k_forests(self, order): 
        #hack - i messed up my symbols
        if self.kirchhoff != 0:
            p = self.kirchhoff.coeff_monomial(Symbol("z")**order)
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
        
        first = [self.edge_complement(k) for k in self.k_forests(order=1)]
        #needs to have at least three vertices for this to make sense - otherwise no way to disconect the vertices and still keep an edge
        second = [self.edge_complement(k) for k in self.k_forests(order=2)]
        
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

