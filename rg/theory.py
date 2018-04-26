import pandas as pd
import operator
import functools
import itertools
import warnings
from scipy.special import factorial
import numpy as np
import sympy
from sympy import sympify
from sympy import Symbol
from sympy import Matrix
from sympy import solve
from sympy import init_printing

from sympy import latex
sympy.init_printing()

from . interaction import  interaction as J, dim, T, L, INV_MEASURE, IN_FIELD, OUT_FIELD

class matrix_accessor(object):
    def __init__(self,theories): self.ts = theories
    @property    
    def general_form(self): return Matrix([ t.general_form for t in self.ts ])
    def criterion(self, d=4): return Matrix([ t.criterion(d) for t in self.ts ])
    def dimensionless_terms(self): return Matrix([  [c.display() for c in t.dimensionless ] for t in self.ts if t.is_valid]) 
#TODO: ADD this one which presents all the couplings side by side for all the theories
#dfs = []
#for k in cs:
#    M = k.interpret_couplings(interactions)
#    df = pd.DataFrame(M.tolist(), columns=["vertex", "coupling"]).set_index("vertex")
#    dfs.append(df)
#dfs = pd.concat(dfs,axis=1).reset_index()
#return Matrix(dfs.as_matrix())
    
class ftheory(dict):
    def __init__(self, *args, **kwargs):
        self.size = None
        self.update(*args, **kwargs)
        self.is_valid = False
        self.dimensionless = []
        if "dimensionless" in kwargs:
            self.dimensionless = kwargs["dimensionless"]
            for k in kwargs["dimensionless"]:  self[k] = 1        
        self.fields = None
        try:
            self.fields = self.solve()
            if len(self.fields) == np.array(self.size).prod():self.is_valid = True
        except Exception as ex: 
            #print("failed to solve", repr(ex))
            pass
        
    @property
    def general_form(self):
        field_weights = [self.fields[k] for k in self.__base_fields__]
        field_weights = [ v.subs({T : L**2}).as_powers_dict()[L] for v in field_weights]
        return field_weights
    
    def criterion(self,d=dim):
        field_weights = [self.fields[k] for k in self.__base_fields__]
        field_weights = np.array([ v.subs({T : L**2, dim:d}).as_powers_dict()[L] for v in field_weights]) * - 1 #invert for the coupling form
        inv_mes_pow= [INV_MEASURE.subs({T : L**2, dim : d}).as_powers_dict()[L]]
        return list(field_weights) + inv_mes_pow
    
    def __getitem__(self, key):
        val = dict.__getitem__(self, key)
        return val

    def __setitem__(self, key, val):
        if isinstance(key,J):
            size = key.tensor.shape
            if self.size == None: self.size = size
            assert size == self.size, "Illegal action: the interactions must all have the same rank. For example a 2*2 to matrix for 2 species with creation-annhilation fields"
            dict.__setitem__(self, key, val)

    def display(self):
        return Matrix( [[ k.display(),v] for k, v in self.items()] )
    
    @property
    def __dimensions__(self):
        d = {}
        for k,v in self.items():
            d[k] = INV_MEASURE/v
        return d
    
    @property
    def __system__(self):
        d = []
        for k,v in self.__dimensions__.items():   d.append(k.display()-v)
        return d
    
    @property
    def dimensions(self):
        d = self.__dimensions__
        return Matrix( [[ k.display(),v] for k, v in d.items()] )
    
    @property
    def __base_fields__(self):
        temp, l = J(np.ones(self.size,np.int) ), []
        for species in range(temp.tensor.shape[0]):
            l.append(temp.in_fields[species])
            l.append(temp.out_fields[species])
        l = np.array(l)
        return Matrix(l.reshape(temp.tensor.shape))
    
    def update(self, *args, **kwargs):
        for k, v in dict(*args, **kwargs).items():  
            if isinstance(k,J):  self[k] = v
           
    def __repr__(self):  return self._repr_latex_()
    
    def _repr_latex_(self):
        init_printing(use_latex='mathjax')
        return latex(self.display(),  mode='inline')
    
    def solve(self):
        return dict(zip(list(self.__base_fields__), solve(self.__system__, list(self.__base_fields__))[0]))
    
    def interpret_couplings(self, interactions, as_matrix=True, l_power_dim=None):
        res = {}
        for k,v in self.interpret_dimensions(interactions,as_matrix=False).items(): res [k] = INV_MEASURE/v
        
        if l_power_dim != None:
            for k,v in res.items():
                res[k] = v.subs({T : L**2, dim : l_power_dim})
           
        if as_matrix:return Matrix([ [k.display(),v] for k,v in res.items()  ])
        return res

    def interpret_dimensions(self, interactions, as_matrix=True):
        res = {}
        if not isinstance(interactions, list): interactions = [interactions]
        field_keys = list(self.__base_fields__)
        field_weights = [self.fields[k] for k in field_keys]
        for k in interactions:
            ordered_weights = list(k.tensor.flatten())    
            #useful test of ordering is to use field_keys instead of field_weights below
            res[k] = functools.reduce(operator.mul, [field_weights[i]**ordered_weights[i] for i in range(len(field_keys))])
        if as_matrix:return Matrix([ [k.display(),v] for k,v in res.items()  ])
        return res

    def uncoupled(interactions, couplings):
        unc = [c for c in interactions if c not in couplings]
        return unc

    def combinations(interactions,couplings, comb_order=3):
        c = ftheory.uncoupled(interactions, couplings)
        L=list(itertools.combinations(list(c), comb_order-1))
        return [list(c) for c in L]

    def theories(interactions,couplings, comb_order=3, valid_only=True):
        cs = ftheory.combinations(interactions, couplings, comb_order)
        ts = [ftheory(couplings.copy(), dimensionless=c) for c in cs]
        if valid_only: ts = [t for t in ts if t.is_valid]
        return ts
    
    def matrices(theories): return matrix_accessor(theories)
    
    