#wrap an incidence matrix and some species info
#make a propagator product
#apply omega integration
#simplify k polynomials
#graph polynomials
#etc.

#temp: later pick your objects - just is just to prototype
from sympy import *

import operator
import functools
import math
import numpy as np

from . graph import scalar_propagator,pole_descriptor

from . import _sum_, _product_ 

class integral(object):
    """
    wrapper around integrand function used for reducing down the omega integrations followed by ks
    """
    def __init__(self, f): 
        if f != 0:
            self.value = f
            is_loop_freq = lambda t : "omega'" not in str(t) and "omega" in str(t)
            #get all the bits i.e. we parse each time the function is passed in
            self.int_vars = [t for t in f.free_symbols if is_loop_freq(t)]

       
    ##all the causal work is done here dont use parser of scalar propagator which so far only works for init of elementary propagators
    def poles(self,v): 
        bottom = self.value.as_numer_denom()[-1]
        for f in bottom.as_ordered_factors():
            keys, values = zip(*f.as_powers_dict().items())
            #assert only one
            exp,order = simplify(keys[0]),values[0]
            mass = solve(exp,v)
            if len(mass) > 0:
                #the pole descriptor needs the f to divide in contour integral, the mass, the pole order and the imaginary half-plane
                yield pole_descriptor(f**-1, mass[0], order, exp.coeff(I*v))
        
    
    def __call__(self,eval_point=None): return self.value if eval_point == None else self.value.subs(eval_point)
    def _sympy_(self): return self.value
    def __repr__(self):  return self._repr_latex_()
    def _repr_latex_(self):  return latex(self.value,mode='inline')
    
    def residues(self, v):
        """
        - everything below here should be pure sympy 
        - for self and p to work they need to be sympified
        - also all ops should preserve propagator factor structure
        - in the division, (self.p), g should still be type integrand? if not, diff some other way
        """
        for pd in self.poles(v):
            #(n-1)th derivitive of G/P evaluated for the pole of order n
            if pd.plane == -1: #'neg freq'
                #print("adding residue for", v)
                g = simplify(self/pd.propagator_factor)
                #print(g)
                if pd.order > 1:g = diff(g,v,n=pd.order-1)
                #print("eval", str(g), "at", pd.mass)
                val = pd.factor * g.subs({v: pd.mass})
                #print("val", val)
                yield _product_([p.simplify() for p in val.as_ordered_factors()])
            
    def integrate(self,int_vars=None, stages=[]):
        int_vars = self.int_vars if int_vars == None else int_vars
        
        for v in int_vars:
            #print("integrating ", v)
            res = simplify([r for r in self.residues(v)])
            stages.append(res)
            self.value = simplify(_sum_(res))
        return self
    
