import numpy as np
import pandas as pd
import itertools
#from scipy.special import factorial
#import sympy
from sympy import Symbol, symbols
from IPython.display import display, Math, Latex
from functools import reduce
from operator import add,mul

class cpropagator(object):
    
    def __init__(self, pid=0, power=1,causal_type=1, D=0, show_subscript=True, use_compact_symbols=False ):
        w,k,  = symbols("\omega k")
        qual = str(pid) if show_subscript else ""
        m = Symbol("m"+qual)
        D = Symbol("D"+qual)
        A = Symbol("A"+qual)
        self.Pinv = causal_type*I*w+D*k**2+m if use_compact_symbols==False else -I*w+A
        self._pole = solve(self.Pinv, Symbol("\omega"))[0]
        self.Pinv = self.Pinv**power
        self._pole_order = power
        self._causal_type = causal_type
    def display(self):return self.Pinv**-1   
    
    #todo - define equality on Pid and causality so that we can make basis sets
    @property
    def value(self): return self.Pinv**-1
    @property
    def pole(self): return self._pole
    @property
    def pole_order(self): return self._pole_order
    @property
    def causal_type(self): return self._causal_type
    def __repr__(self):  return self._repr_latex_()  
    def _repr_latex_(self):   return latex(self.display(),  mode='inline')
    
    def _pset__():
        pass #some helper to make sets but compact into powers of propagators - that way we can ask for combinations of higher orders but limit types with k
    
    def __display__(pset):
        return [p.display() for p in pset]
    
    def generate_propagator_basis(n,power=1,combination_order=1,filter_loop_degenerates=True,for_display=False):
        """alot going on here - we can choose n distinct or n of k types of propagators 
        - we can combine them into products of propagators which are compacted into powers leading to multiple poles"""
        Ps = []
        for i in range(n):
            for p in range(0,power):
                for causality in [-1,1]: Ps.append(cpropagator(i, power=p+1, causal_type=causality))
         
        L =  Ps#list(itertools.chain(*Ps))
        if len(L) < combination_order:  combination_order = len(L)
        if combination_order == 1:  return [cpropagator.__make_list__(L,for_display)] #make it a list of lists for consistency
        else:
            res = []
            for i in range(2,combination_order+1):
                combo = list(itertools.combinations(L,i))#[t for t in  ]
                #i think this is where i check they are non degenerate i.e. cannot all have the same causal type if they are to be in a loop by momentum conservation
                combo = [cpropagator.__make_list__(c,for_display) for c in combo  ]#if not cpropagator.__is_degenerate_loop__(c)
                
                res.append(combo)
            return res
    
    def __make_list__(pset,display):
        return list(pset) if not display else cpropagator.__display__(pset) 
        
    def __is_degenerate_loop__(pset):
        return len(np.unique([c.causal_type for c in pset])) == 1      
    
    def __eval_combinations__(combinations):
        res_set = []
        for order in combinations:
            res = []
            for pset in order:  
                #we add the sum of residues knowing terms in the definition cancels with integration measure
                res.append(reduce(operator.add, [r["integration"] for r in cpropagator.residues(pset)], 0))
            res_set.append(res)
        return res_set
       
    def residues(pset,return_smallest=True):
        W = Symbol("\omega")
        def _diff_(s,o):return s
        def sum_poles(st): return np.sum([p["order"] for p in st])
            
        sets = {1:[], -1:[]}
        residual_base = reduce(operator.mul, [p.value for p in pset], 1)
        for p in pset:
            d = {"pole":p.pole}
            d["order"] = p.pole_order
            d["residual"] = simplify(residual_base / p.value)
            #generally we take the mth derivitive of g and evaluate there, where m is pole order
            if p.pole_order > 1:d["residual"]= diff(d["residual"],W,n=p.pole_order)
            d["causal_type"] = p.causal_type
            d["integration"] = d["residual"].subs(W,p.pole)
            sets[d["causal_type"]].append(d)
        #this could become a more complex choice; maybe if the smaller set gives rise to imaginary derivs/residuals I dont want them?
        #unless by some cheeky invariance argument i can ignore the causality here - but I doubt it.
        #if they are the same size, take the one with the fewest poles - probably what i should have done anyway
        #if len(sets[1]) == len(sets[-1]): return sets[1] if sum_poles(sets[1]) < sum_poles(sets[-1]) else sets[-1]
        return sets[1] if len(sets[1]) < len(sets[-1]) else sets[-1]
    
    