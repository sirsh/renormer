"""

work in progress propagator handling - still setting scope and learning sympy interface
There is a basus for creating propagators or one can use the constructor. For example the vector basis
contains coefficients for say two loops, two externals and then species and power details


mat = [[1,-1,0,1,1,2],
       [0,-1,0,0,0,1]]
factors = scalar_propagator.from_basis_vector(mat,2) # the 2 is the loop order needed for this quick tensor constructor convention

# the expression can also be parsed but maybe there is a better way
here there are two propagators which should be the ordered factors in the expression
P1 = list(scalar_propagator.parse( factors ))[0]
P2 = list(scalar_propagator.parse( factors ))[1]

#I can ask for poles in short hand for the omega loop variables or supply the Symbol omega_1 etc.

P1.pole(1)

"""


from sympy import *

import operator
import functools
import math
import pandas as pd
import numpy as np

def _sum_(terms): return functools.reduce(operator.add, terms , 0)
def omega_var(i,is_external=False):return I*Symbol("omega_"+str(i)) if not is_external else I*Symbol("\omega'_"+str(i))
def k_var(i,is_external=False):return Symbol("k"+str(i)) if not is_external else Symbol("k'_"+str(i))
def __prod__(l): return functools.reduce(operator.mul, l,1)


class pole_descriptor(object):
    def __init__(self, propagator, mass, order,plane):
        self.mass = mass
        self.order = order  
        self.prop = propagator
        self.plane = plane
        # I have messed up my I orders here because I need to check with the deriv is doing
        #there is a trailing I hypothesis that is however many derivitive you have, e.g. order-1 giving i's, you also have order of i's of the same sign leaving you with i
        #no matter we do i think each residue should have a factor of -i which then cancels with the one in the formula
        factor =  (  ( plane*-1*I**(order) ) * ( (-1*I) **(order-1) )  ) if order > 1 else Symbol(str(1))
        self.factor = factor/np.math.factorial(order-1)
    
    def __repr__(self):  return self._repr_latex_()
    def _repr_latex_(self): return latex(self.mass,  mode='inline')
    #if evalauted, this returns the propagator itself - maybe this is confusing but...
    def __call__(self):return self.value
    def _sympy_(self): return self.value
    
    @property
    def propagator_factor(self):return self.prop
    
#this guy is solid but to signs in omegas. The cycle basis MUST produce the correct edge bases and also signs used for omegas
class scalar_propagator(object):
    def __init__(self, loop_basis, ex_edge_basis, index=0, species="a",  power=1):
        species_variable = lambda name : Symbol(name+"_"+species)
        self._index = index
        self._species = species
        self._power = power
        k_vars = [ abs(v) * k_var(i) for i,v in enumerate(loop_basis)] + [ abs(v) * k_var(i,True) for i,v in enumerate(ex_edge_basis)]
        o_vars = [ v * omega_var(i) for i,v in enumerate(loop_basis)]  + [ v * omega_var(i,True) for i,v in enumerate(ex_edge_basis)]

        self._core =  _sum_(o_vars)+species_variable("D")*(_sum_(k_vars)**2) +species_variable("m")**2
                    
        self._den = self._core**power
                
        self._value = self._den**(-1) 
        
        #self._associated_delta_functions = [signed_delta("omega",index)._rep,signed_delta("k",index)._rep]

    def pole(self, omega):
        #convenience option of testing - extracts the loop variable by index
        if isinstance(omega,int):omega = Symbol("omega_"+str(omega))
        
        res = solve(self._core,omega)
        #zero or one solutions for each factor in the 
        return None if len(res) == 0 else pole_descriptor(self, res[0], self._power,plane=self._core.coeff(I*omega))
    
    def from_basis_vector(mat, loops, lower_species=True):
        """
        temporary helper to map a matrix to an integral - assumes you ass L and there are less than 7 species and uses
        [loop_basis|edge_basis|power|species_id]
        example mat for two loops
    
        mat = [
            [0,1,0,1,1,1],
            [0,-1,0,0,2,2],
        ]
        """
        species = ["A", "B", "C", "D", "E", "F", "G", "H"]
        if lower_species: species = [s.lower() for s in species]
        n = len(mat)
        props = []
        for i, v in enumerate(mat):           
            props.append(scalar_propagator(loop_basis=v[:loops],
                                           ex_edge_basis=v[loops:-2],
                                           index = i, 
                                           power = v[-1],
                                           species = species[v[-2]]).value)
         
        return __prod__(props)
           
    def _parse_(exp,i=0):
        """A horrible method to parse the sympy object. Must improve this"""

        keys, values = zip(*exp.as_powers_dict().items())
        #assert only one
        exp,power = simplify(keys[0]),values[0]
        
        s = None
        ex,ex_signs,loop,loop_signs = [],[],[],[]
        sign_dict = exp.as_coefficients_dict()
        #print(exp.free_symbols, "used in parse")
        for t in exp.free_symbols:
            if str(t)[:1] == "m":s = str(t)[-1].lower()
            elif  "omega'" in str(t): 
                ex.append(int(str(t)[-1]))
                ex_signs.append(sign_dict[I*t])
            elif  "omega" in str(t): 
                loop.append(int(str(t)[-1]))
                loop_signs.append(sign_dict[I*t])   

        if len(loop)==0:return None # if this happens we may be dropping a factor
        
        loop_base = list(np.zeros(np.max(loop)+1,np.int))
        ex_base = list(np.zeros(np.max(ex)+1,np.int)) if len(ex) > 0 else []
        for i,e in enumerate(ex): ex_base[e] = ex_signs[i]
        for i,e in enumerate(loop): loop_base[e] = loop_signs[i]

        return scalar_propagator(loop_base,ex_base,index=i,species=s,power=abs(power))

    def parse(exp): 
        for i, f in enumerate(exp.as_ordered_factors()):
            sp= scalar_propagator._parse_(f,i)
            if sp != None: yield sp
            
    @property
    def value(self):
        
        ret = self._value #* functools.reduce(operator.mul, self._associated_delta_functions , 1)
        return ret
        return ret.simplify()
        
    def __repr__(self):  return self._repr_latex_()
    def _repr_latex_(self): return latex(self.value,  mode='inline')
    def __call__(self, i):pass
    
    def subs(self, d):
        #as below for deltas, we need to keep track of terms seperately to avoid rediscovery, which needs cleaning up
        self._den = self._den.subs(d)
        self._value = self._value.subs(d)
        return self
    
    def apply_delta(self, d,i):
        #hacky - I have not thought out how these workflows should be done yet but I will
        self._den = self._den.subs(d.substitution(i))
        self._value = self._value.subs(d.substitution(i))
        
        return self
    
    def apply_delta_functions(self):
        #hack
        self._associated_delta_functions = []
        if self._sense == 1: 
            self._value=self._value.subs({"omega_"+str(self._index): -1*Symbol("omega_"+str(self._index)) })
            #hack see above
            self._den = self._den.subs({"omega_"+str(self._index): -1*Symbol("omega_"+str(self._index)) })
            
        return self
    
    def Ltransform(self):
        self._value = self._value #* exp(-I*Symbol("omega_"+str(self._index))*Symbol("t_"+str(self._index)))
        return self
    
    def Ftransform(self):
        self._value = self._value * exp(-I*Symbol("omega_"+str(self._index))*Symbol("t_"+str(self._index)))
        return self