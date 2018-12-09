from . import incidence_matrix
import functools, itertools, operator
import numpy as np

def coloured_edges(self, direction = 1,stack_colours=False, as_c_residual=True):
    spec = np.array(self.species_vector)
    edges = np.where(self[-1]==direction)[0]
    vids = np.where(self.T[edges].T[:-1] != 0)[0]
    
    #if it is as_c_residual and we have a large radical, we only consider head and tail edges either index 0 or V-2 by convention
    if as_c_residual and self.shape[0] > 3: 
        edges = [e for i, e in enumerate(edges) if vids[i] in [0,self.shape[0]-2]] 
        #recompute
        vids = np.where(self.T[edges].T[:-1] != 0)[0]
  
    if stack_colours:
        d = {}
        for i,e in enumerate(edges):
            s = spec[e]
            if s not in d: d[s] = []
            d[s].append(e)
        return d
    return list(zip(edges, vids, spec[edges]))

def binnings(n, vs): return set(list(itertools.combinations(vs,n)))
    
def expand(be,binf):
    items = []
    for v, c in enumerate(be):
        #for a particular vertex, maybe there are some data - get the first c rows
        ar = binf[np.where(binf[:,1]== v)[0]][:c,0]
        #if there are rows, yield the edge ids for them
        if len(ar) > 0:items +=  list(ar)
    return items
    
def radical_shuffle(A,B,maximal=False):
    ac = coloured_edges(A, -1,True)
    binf = np.array(coloured_edges(B, 1))
    #get spec capacities on A LHS
    #print(ac)
    #print(binf)
    all_edges = {}
    for s in ac.keys():
        sc = len(ac[s])
        #print(s,sc)  
        finfo = binf[np.where(binf[:,-1]==s)[0]]
        #print(binf)
        #for each species generate the labelled objects on v for each species and choose combinations from the vset with sets
        vset = list(finfo[:,1])
        bins = binnings(sc,vset)
        #print(s,sc,bins)
        #map these to edge choices - be careful to use the occurances of each vertex in the tuple and index the edges that way
        bincounts = [list(np.bincount(b)) for b in bins]
        #expand on the edges
        all_edges[s] = [tuple(expand(e,finfo)) for e in bincounts]
    
    #cartesian product over all species edges
    all_edges = [ set(functools.reduce(operator.add, t , ())) for  t in list(itertools.product(*all_edges.values()))]
 
    if maximal: return all_edges
    
    all_sets = []
    for t in all_edges:
        #t = list(t)
        for i in range(1,1+len(t)):
            s= list(itertools.combinations(t,i))
            all_sets += s

    return set(all_sets)

#what is the chord set
def chord_set_and_edge_stack(self, r, ):
    sbasis=np.zeros(2,np.int)
    spec = np.array(self.species_vector)
    edges = list(r)
    spec = spec[edges]
    outd = {}
    for s in range(len(sbasis)): outd[s] = []
    for i,e in enumerate(edges):
        s = spec[i]
        sbasis[s]+=1
        outd[s].append(e)
    return sbasis,outd

def shuffle_product(A,B):
    #try:
    res = radical_shuffle(A,B)
    for r in res:
        chords, outd = chord_set_and_edge_stack(B,r)
        yield A.dprod(B,species_chord_set=chords,outd=outd)
    #except:print("shuffle failed")
        
        
        
#first pass
def binary_combinations(prim, left_right_restriction=None,skip_loops=True):
    idx = list(range(len(prim)))
    mulfodder= [tup for tup in list(itertools.permutations(idx,2))]
    def make(comb): return tuple(prim[i] for i in comb)
    for m in mulfodder:
        A,B = make(m)
        if skip_loops and (A.first_betti_number >= 1 or B.first_betti_number >= 1): continue
        if left_right_restriction is not None:
            #if there are more vertices than allowed skip - i am doing greater than to be more exhasutive - could be not equal to on RHS of check
            if A.shape[0] != left_right_restriction[0] or B.shape[0] > left_right_restriction[1]: continue
                
        for r in shuffle_product(A,B):
            yield r
            
            
def iso_reduce(l,add_new_residues=True, max_residue_legs=3, max_grahp_legs=4):
    u = {}
    for r in l: 
        rp = r.residue_part
        #only add if the residue is less than the dim any
        if rp.shape[-1] <= max_grahp_legs:
            u[r.graph_hash()] = r
        if add_new_residues:
            if rp.shape[-1] <= max_residue_legs and  rp.shape[-1] > 2: #should not have or add propagator types
                u[rp.graph_hash()] = rp
            
    return list(u.values())


    
def sep_loops(l):
    d ={0:[], 1:[], "primitives": [], "loops":[]}
    for r in l:
        if r.is_tree: 
            if r.shape[0] <= 2:
                d["primitives"].append(r)
            else:   d[0].append(r)
        if r.first_betti_number==1 and r.bridge_count == 0: 
            d[1].append(r)
    
    loops = {}
    for c in d[1]:   
        lh = c.loop_hash()
        if lh not in loops: loops[lh] = { "graph" : c.loop_part, "counter" : 0 }
    for v in loops.values():
        d["loops"].append(v["graph"])
    
    return d
    

    
#maybe find somew
#combinations over the primitives choose 4           
def _shuffle_(A,B,contract_tree_level=False):
    for s in incidence_matrix.shuffle(A,B):
            
        #if A or B are tree level, then only use the residual in the operation - this is true but need a clean way to describe it
        #internal edges are contracted always unless there is a loop
        
        #TODO this should be a sum over arbitrary species basis
        if s[0] + s[1] != 0: #the null case     
             yield A.dprod(B,species_chord_set=s, contract_tree_level=contract_tree_level)
                
def brute_mul(A,B,C,D,l=1):
    genset = list(_shuffle_(A,B))
    news = []
    for g in genset:news += list(_shuffle_(g,C))
    genset+= news
    for g in genset:news += list(_shuffle_(g,D))
    genset+= news
    return  [g for g in genset if g.first_betti_number == l and g.bridge_count == 0]# and



def brute_mul_hack(L,l=1):
    A,B = L[0], L[1]
    genset = list(_shuffle_(A,B))
    news = []
    if len(L) > 2:
        C = L[2]
        for g in genset:news += list(_shuffle_(g,C))
        genset+= news
        if len(L) > 3:
            D = L[3]
            for g in genset:news += list(_shuffle_(g,D))
            genset+= news
            
    return  [g for g in genset if g.first_betti_number == l and g.bridge_count == 0]# and