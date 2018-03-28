#TODO
#INTERACTION SYSTEM SHOULD BE FUNCTIONAL AND MULT-THEORY CENTRIC
#THE RESOLUTION IS WHEN ALL BASE FIELDS ARE KNOWN - ONE CAN CREATE A BASIS - THERE MAY BE RESIDUAL PROLONGATIONS - WE SHOULD CHECK WHICH WERE IN THE ORIGINAL
#PROBABLY ALL COMMON FIELDS ARE THE ORIGINAL _ BECAUSE I DO NOT THINK WE STOP UNTIL WE FILL IN - BUT NOT REQUIRED
#SHOULD BE EASIER TO SEE RESOLUTION PATHS; I think 0 happens when one of the original vertices is bypassed when resolving the base fields.
#


#also a class for diagrams that takes interactions or constructions
#constructions are nested interactions with some more information -a construction can have an effecetive interaction for validation
#constructions can have k-order loops from 0
import pandas as pd
import operator
import functools
import itertools
import warnings
from scipy.special import factorial
warnings.filterwarnings("ignore")
#later come back and test all the ops without this and clean up

import numpy as np
import sympy
from sympy import sympify
from sympy import Symbol
from sympy import Matrix
from sympy import init_printing
from sympy import latex
sympy.init_printing()

from IPython.display import Latex#
#special symbols

from . import diagrams

dim = Symbol("d")
T = Symbol("T")
L = Symbol("L")
INV_MEASURE = 1/ ( T * (L**dim))


def __display_all_(interactions, show_mat=False):
    return [d.display(show_mat) for d in interactions]

# class coupling(list):
#     def __init__(self, *args):
#         list.__init__(self, *args)
        
class _interaction_identity(object):
    def __and__(self,other):  return other
    def __add__(self,other): return other
    #def __sub__(self,other): return other
    
interaction_identity = _interaction_identity()

#add sympy and see what we can do with it
#must generate symbols for the fields
#assumed all fields can be equated to the inverse time-space measure because [L][D]=1
class interaction(object):
    #exnend a grammar for the matrix ops so they are hidden away
    ##todo - global system to map field_index to field_letter or new symbol as a configuration option
    def __init__(self, mat, coupling_symbol=None):
        self._mat = mat
        self._tensor = np.array(mat)
        self.shifted=True
        self.symmetry_factor = 1#not sure if this is the one for the externals??
            
        
        self.in_fields = []
        self.out_fields = []
        out_ = lambda power,index : Symbol('\\tilde{\phi}_'+str(index))**power
        in_ = lambda power,index : Symbol("\phi_"+str(index) )** power
        ar = np.array(mat)
        
        self.num_species = ar.shape[0]
        
        for s in range(ar.shape[0]):
            self.in_fields.append(in_( ar[(s,0)] , s))
            self.out_fields.append(out_( ar[(s,1)] , s))
            
        #display(Latex(self.__repr__()))

        
    def __common_fields__(self, other):#not sure if this should take cardinality into account - no harm if we are doing gradual reduction
        #easy way to do this would be to subtract them from me. we have count non negative fields in common + the non zero sum
        #SELF - OTHER + 1 summed iun numpy would count preceisely

        ar = (np.array(self._mat) - np.array(other._mat) + 1)
        ar[ar<0]=0
        return ar.sum()
    
    @property
    def tensor(self): 
        
        return self._tensor
    
    @property
    def is_trivial(self): return np.array(self._mat).sum() == 1 and np.count_nonzero(self._mat)== 1
    
    @property
    def is_constant(self):  return np.count_nonzero(self._mat)== 0
    
    @property
    def is_positive(self):  
        ar = np.array(self._mat)
        return len(ar[ar<0]) == 0
    
    def __hash__(self):  return hash(str(self._mat))

    def __eq__(self,other): return np.array_equal(self._mat, other._mat)
    
    def __and__(self,other):
        if isinstance(other,_interaction_identity):return self
        c = compound_interaction(self,other)
        print("symmetry factor:",c.symmetry_factor)
        return c.effective_interaction
        
    def __add__(self,other):
        a = np.array(self._mat) + np.array(other._mat)
        #a[a < 0] = 0
        return interaction(a)
    
    def __sub__(self,other):
        a = np.array(self._mat)- np.array(other._mat)
        #a[a < 0] = 0
        return interaction(a)
    
    @property
    def symmetries(self):
        left_external = np.array(self._mat)[:,1]
        return factorial(left_external, exact=True).prod()
    
    @property
    def diagram(self):  return diagrams.diagram(self)
    
    def shift(self):  pass #returns new set of fields
    
    def to_generator(self):   pass 
    
    def display(self, show_mat=False, printer='mathjax'):
        init_printing(use_latex=printer)
        if show_mat:  return sympy.Matrix(self._mat)
        else:  return functools.reduce(operator.mul, self.in_fields + self.out_fields, 1)
        
    def __repr__(self):  return self._repr_latex_()
    
    def _repr_latex_(self):
        init_printing(use_latex='mathjax')
        
        return latex(self.display(False),  mode='inline')
    
#     def _repr_html_(self):
#         init_printing(use_latex='mathjax')
#         latex(self.display(True))
    
    def __truediv__(self,b): return self.__div__(b)
    def __floordiv__ (self,b):  return self.__div__(b)
    def __rdiv__(self,b):   return self.__div__(b)
    def __rtrue_div__(self,b): return self.__div__(b)
    def __div__(self,b):  
        other = np.array(b._mat)
        mine = np.array(self._mat)
        art = mine/other
        #print(art)
        art[np.array(other==0)]=0
        remains = art[art != 0]
        return int(remains.min()) if len(remains) > 0 else 0
        
    def __mod__(self,b):  
        coef = self/b
        other = np.array(b._mat) * coef
        return self - interaction(other.astype(int))
           
    #assert - applying unresolved complexes against resolved interactions
    #assume therefore that self has inverse measure weight and also we are producing enough alternate factors
    
    def reduce(self, coupled, new_couplings = {}):
        interaction._reduce(self, coupled, new_couplings)
        
        #check if all base fields exists in the system and return true or false for resolved
        
    
    #only coupled can take part at all!!!!! until we find all trivial interactions and then find the final lambdas
    #process is (ANSATZ | DEGENREATE COUPLINGS) | FIND ALL TRIVIAL INTERACTIONS i.e. FIELDS | RESOLVE FINAL 
    def _reduce(rem, coupled, new_couplings = {}):
        for I,v in coupled.items():#
            divides = rem / I
            if divides > 0:               
                rem_copy = rem % I # self has inverse measure weight
                
                #skip negative for now
                if not rem_copy.is_positive: continue
                
                #print("reduction", rem, "by", I, "is", rem_copy)
                if rem_copy not in new_couplings and rem_copy not in coupled  and not rem_copy.is_constant: #only if we have never seen it (could try to reduce again but staying atomic for now)
                    #this is a simple statement - wehn we devide in the known factor *weight, the remainder has a couplings of that same factor
                    # \phi[\phi\psi]   =>    \phi[\lambda], where \lambda is whatever we factored out
                    res = (interaction_system.__invert_coupling__(v)   ** divides)
                    if res != 0:
                        new_couplings[rem_copy] = res.simplify()  
                        interaction._reduce(rem, dict(new_couplings), new_couplings) # this makes no sense because rem is a coupled field not an uncoupled field
                        
    
class interaction_system:
    def __init__(self, interactions):
        self._interactions = interactions
        self.__reset__()
        
    def __reset__(self):
        self.couplings = dict(zip(self._interactions, [0 for i in range(len(self._interactions))]))
        self.num_species = self._interactions[0].num_species #check consistent basis TODO
        self.theories = []
        
    def __display_couplings__(self):
        return [[a.display(), b] for a,b in self.couplings.items()]
    
    def __display_interaction_dims__(self):
        return [[a.display(), b] for a,b in self.inverse_couplings.items()]
    
    #this doesn matter if its the wieght or coupling COUP*WEIGHT=INV_MEASURE so they are adjoint and same diff
    def __invert_coupling__(lamb):
        return 0 if lamb == 0 else INV_MEASURE/lamb
        
    def invert_couplings(chi):
        d = {}
        for a,b in chi.items():  d[a] = interaction_system.__invert_coupling__(b)
        return d
        
    def get_basis(self,scaled=None):
        items = []
        for idx, k in enumerate(range(self.num_species)):
            L = [None,None]
            ar = np.zeros((self.num_species,2),np.int)
            ar[k, 0] = 1
            L[0] = self.inverse_couplings[interaction(ar)]
            if scaled is not None:  L[0] = L[0]**scaled[k, 0]
            ar = np.zeros((self.num_species,2),np.int)
            ar[k, 1] = 1
            L[1] = self.inverse_couplings[interaction(ar)]
            if scaled is not None:  L[1] = L[1]**scaled[k, 1]
            items.append(L)
        return Matrix(items)

    @property
    def inverse_couplings(self): return interaction_system.invert_couplings(self.couplings)
    
    @property
    def coupled_interactions(self):
        d = {}
        for a in  self._coupled_interaction_filter(self.couplings):   d[a] = self.couplings[a]
        return d
    
    @property
    def uncoupled_interactions(self):
        d = {}
        for a in  self._uncoupled_interaction_filter(self.couplings):   d[a] = self.couplings[a]
        return d
                
    def _coupled_interaction_filter(self, K):
        for k in K.keys():
            if K[k] != 0:   yield k
                
    def _uncoupled_interaction_filter(self, K):
        for k in K.keys():
            if K[k] == 0:   yield k
                
    def display(self, show_mat=False):
        return [d.display(show_mat) for d,c in self.couplings.items()]
    
    
    #returns as last arg the value of the fields
    def is_biproduct_of_known_interactions(I,CIs):
        L=list(itertools.combinations(list(CIs.keys()), 2))
        for t in L:
            P = t[0].__add__(t[1])
            if P == I:
                v1 = interaction_system.__invert_coupling__(CIs[t[0]])
                v2 = interaction_system.__invert_coupling__(CIs[t[1]])

                return (t[0], t[1], v1 * v2)
        return None
    
    
    def add_trivial_adjoint_of_known_interactions(self, I):
        result = None
        new_couplings = {}
        for c in self.couplings:
            res = I % c
            weight = self.inverse_couplings[I]
            if res not in list(self.couplings.keys()):
                if res.is_trivial and res.is_positive:
                    #having devided my self by something we know, i find the residual is a simple field that is not known
                    #i have a field dimension and now i will div by the weight if this other field to find the residual weight
                    #print("trivial:", res.is_trivial, "postive", res.is_positive, "const", res.is_constant, "new:", res ) 
                    result =  weight/ self.inverse_couplings[c]
                    #i save it as a new coupling
                    new_couplings[res] = interaction_system.__invert_coupling__(result)
        
        self.couplings.update(new_couplings)
     
    def update_theories(self, new_couplings, max_L=None):
        self.__reset__()
        self.couplings.update(new_couplings)
        L = list(self.uncoupled_interactions.keys())
        max_L = len(L) if max_L == None else max_L

        self.theories = []
        #TODO - count equations and unknowns - how many ansatx-lmbdas do we need, max combinations there
        #ACTION TERMS: equations, unknowns= Species*2 i.e. 4 - |NormalAnsataz| therefore TERMS - 4 + X 
        #One coupling unkown for each equation + costant 4 more. 
        for i in range(1,max_L):
            #perms could be used for all orders and test for equiv class
            for tup in list(itertools.combinations(L,i)):  
                P= list(tup)
                nt = new_couplings.copy()
                for p in P: nt[p] = 1
                #create a new system - bad design - fix later
                IS_temp = interaction_system(self._interactions)
                #use this template of ansatz
                IS_temp.couplings.update(nt)
                resolved = IS_temp.reduce(nt)    
                if resolved:
                    self.theories.append(nt)
    
    def __get_dimensionless_fields_for_theories__(self):
        fixed = [[] for i in range(len(self.theories))]
        for i,t in enumerate(self.theories):
            for f in t.keys():
                if t[f] == 1:#dimensionless
                    fixed[i].append(f.display())
        return fixed
    
    @property
    def __equiv_classes__(self,include_base_fields_ony=True):
        mat = self.__get_resolved_fields_for_theories__(base_fields_only=include_base_fields_ony)
        equiv_classes = {}
        keys = list(mat.keys())

        def is_equiv(a,b):
            for k in keys:
                if mat[k][a] != mat[k][b]:return False
            return True
        
        l = len(self.theories)
        
        for i in range(l):
            equiv_classes[i] = [i] #im my own equiv class
            for t_ex in equiv_classes.keys():
                #unless im like something upstream, then map to that
                if is_equiv(i,t_ex) and i != t_ex:  
                    equiv_classes[i].append(t_ex) 
                    equiv_classes[t_ex] = equiv_classes[i]

        #go through each theory and put in an equiv class
        return equiv_classes

    def eval_interaction_for_theory(self, V, idx,should_reduce=True):
        IS = interaction_system(self._interactions)
        IS.couplings.update(self.theories[idx])
        resolved = IS.reduce({})#if we have the base fields i.e. basis - all theories with the same basis are equivalent
        if resolved:  
            res = IS.get_basis(np.array(V._mat))
            if should_reduce: return functools.reduce(operator.mul, res, 1)
    
        return False
                
    
    def __get_resolved_fields_for_theories__(self, as_matrix=False, as_powers_of_L=False, base_fields_only=False, show_relevance_at=None, resolve_all=True):
        ALL = {}

        l = len(self.theories)
        bases = []
        
        theory_couplings = []
        
        for index, t in enumerate(self.theories):
            IS = interaction_system(self._interactions)
            IS.couplings.update(t)
            resolved = IS.reduce({})#if we have the base fields i.e. basis - all theories with the same basis are equivalent
            
            theory_couplings.append(IS.couplings)
            
            if resolved:
                for c in IS.couplings.keys():
                    #only base fields to be shown
                    if base_fields_only and not c.is_trivial: continue
                        
                    if c not in ALL:
                        #place holder for each ordered theory
                        ALL[c] = [None for i in range(l)]
                        
 
                    ALL[c][index] = IS.couplings[c]
                    
        def get_L_power(s,D=4):
            try:#using globals in this file
                if s == None: return None
                LL = s.subs({T : L**2})
                return LL.as_powers_dict()[L]
            except: return 0 #dangerous PLUS need to find out when things are not proper symbols

        def get_rel(s):
            try:#using globals in this file
                if s == None: return None
                LL = s.subs({T : L**2, dim:show_relevance_at})
                P= LL.as_powers_dict()[L]
                return P
            except: return -1 #dangerous PLUS need to
            
        if resolve_all:
            print("resolving for all known fields")
        
            for coup in list(ALL.keys()):
                ALL[coup] =  [interaction_system.__invert_coupling__(self.eval_interaction_for_theory(coup, ind)) for ind in range(len(self.theories))]
      
        if as_powers_of_L or show_relevance_at != None:
            for key in ALL:
                ls = ALL[key]
                if show_relevance_at != None: ls = [get_rel(i) for i in ls]
                else: ls = [get_L_power(i) for i in ls]
                ALL[key] = ls
                

            

        return ALL if not as_matrix else Matrix([ [k.display() ] +  ALL[k]  for k in list(ALL.keys())])
            
    #ideas is to keep trying to find new factors until we fail out
    #we may need to add redundant i.e dimensionless couplings
    #during the process, we should find all trivial interaction couplings which terminates
    def reduce(self, Chi, order=1, finalize_by_products=True):
        #if len(Chi) == 0: raise ("Must set Chi")
        self.couplings.update(Chi)   
        found_new = True
        counter_safety = 0    
        new_couplings = {} 
        
        for I,v in self.coupled_interactions.items(): 
            I.reduce(self.couplings, new_couplings)
            self.couplings.update(new_couplings)
    
        #repeat while finding new things
        
        #update the Js that are products of things we know
        for k in self.uncoupled_interactions:
            res = interaction_system.is_biproduct_of_known_interactions(k, self.coupled_interactions)
            if res != None:#if we find it, we can update the non set guy with the tuple tail which is the product
                self.couplings[k] = INV_MEASURE/res[-1]
           
        for k in self.coupled_interactions:
            self.add_trivial_adjoint_of_known_interactions(k)  
            
        #construct the base field for all species and check that we have a resolved couplings for all
        for k in range(self.num_species):
            ar = np.zeros((self.num_species,2),np.int)
            ar[k, 0] = 1
            if interaction(ar) not in self.couplings:  return False
            if self.couplings[interaction(ar)] == 0 : return False #cannot be zero - recall why sometimes 0
            ar = np.zeros((self.num_species,2),np.int)
            ar[k, 1] = 1
            if interaction(ar) not in self.couplings:  return False
            if self.couplings[interaction(ar)] == 0 : return False #cannot be zero - recall why sometimes 0
            
            
        #if finalize_by_products:
        #    for k in self.uncoupled_interactions:
        #        res = self.get_basis(np.array(k._mat))
        #        self.couplings[k] = INV_MEASURE/res
                
        #at the very end we may have got all our base fields but there may be some by-products for which we did not update couplings
        #we can optinally do this now
                
                
        return True
            
        
        
#manage heirarchy of nodes
#validator is global and if we can validate it, we can add it to the pile 
class compound_interaction(object):
    def __init__(self,nodea,nodeb):
        """
        nodes can be either interactions with a ._mat or they can be wrappers i.e. compound_interaction
        """
        self.left,self.right = nodea,nodeb
        self._symmetry_factor = 1
        self._child_syms = 1
        self.loops = 0
        if isinstance(nodeb, compound_interaction):  self.right = nodeb.effective_interaction
        if isinstance(nodea, compound_interaction):  self.left = nodea.effective_interaction
        self._child_syms = nodea.symmetry_factor * nodeb.symmetry_factor
        self._effective = self.__merge__(self.left,self.right)#
        self._external_sym_factors = 1 # this one is determine by each species external choices
        #foreach OUT species add power!
        left_external = self._effective[:,1]
        self._external_sym_factors = factorial(left_external, exact=True).prod()
        
    def __merge__(self, interaction_l, interaction_r):
        left = interaction_l._mat
        right = interaction_r._mat
        self._symmetry_factor = 1
        residual, self._symmetry_factor, pairings = compound_interaction.pair_matrices(left, right)
        if pairings.sum() > 2: self.loops += 1#i think?
        return residual
        
    def pair_matrices(l,r):
        l,r = np.array(l),np.array(r)
        l_,r_ = l[:,0],r[:,1] #left internal, #right internal
        s,e = l[:,1], r[:,0] #left external, #right external
        binding = np.array([l_,r_]).T
        pairings = np.array([binding.min(1),binding.min(1)]).T 
        binding, pairings       
        residual = binding - pairings #what is left over after binding to add to externals
        has_pairings = (pairings>0).astype(int) #mask
        multiplicity = l_ * has_pairings #left internal used for symmetry factor (for all species)
        res = np.stack([e,s]).T + residual
        return res, multiplicity.sum(), pairings
        
    @property
    def symmetry_factor(self): return self._symmetry_factor * self._child_syms * self._external_sym_factors
    
    @property
    def effective_interaction(self):  return interaction(self._effective)
    
    def _repr_html_(self):  
        r = 4 if self.loops == 1 else 2#hard code to make it loopy for now
        return diagrams.diagram(self.effective_interaction, props={"loop_radius" : r})._repr_html_()
    
    def display_construction(self): pass
    
    def combinations(primitives, loop_orders = [0], max_k=3):
        if not isinstance(loop_orders, list):loop_orders = [loop_orders]
        l = []   
        for i in range(2,max_k+1):
            for tup in list(itertools.permutations(primitives,i)):
                res = compound_interaction._combine(list(tup))
                if res.effective_interaction in primitives and res.loops in loop_orders: l.append(res)
        #todo: define uniquness for compound_interaction set
        return list(set(l))
        
    def _combine(set_of_interactions):
        def _merge(a,b): 
            if isinstance(a, _interaction_identity):return b
            if isinstance(b, _interaction_identity):return a
            return compound_interaction(a,b)
        return functools.reduce(_merge, set_of_interactions, interaction_identity)
    
    
    #this entire method of combinations needs to improve 
    #- we need to generalise between interactrions and compund interactions so they are more like the same thing with the same operations
    #and they should in general be renormalised. The way I have done it here i have sort of tacked it on
    #the way forward is obvious. just remember symmetry factors and loops and add this bominer stuff to the boolean operators
    #uniqueness then defined for the matrix plus the symmetry factors ?? hmmmm - maybe that is a good reason to keep them as seperate things    
    
    
    
    
    
    
    