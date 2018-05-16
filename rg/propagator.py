from sympy import Symbol, symbols, Poly, I,solve, simplify,latex,expand,Function,collect
from IPython.display import display, Math, Latex
import itertools
from functools import reduce
from operator import add,mul
import numpy as np
from sympy import Function

#helpers
def __alpha_params__(ps):
    for b in ps.atoms():
        if b.is_constant():continue
        if b.name[:5] == "alpha":yield b
            
def determine_elimination_variable(p):
    #actually i will implement a comparator: we choose one that has only one entry with the most complex coefficient - here i just assumed we have 1 or 2
    coeffs = p.all_coeffs()
    coeff = coeffs[0]
    
    #this is the condition where we have already got a polynomail that has no linear k term 
    if len(coeffs) == 3 and coeffs[1] == 0:return None
    
    d = {}
    set_rank = {}
    rank = 0
    bad_set = None
    for a in __alpha_params__(coeff):
        co = coeff.coeff(a)
        if co not in d: d[co] = []
        d[co].append(a)
        set_rank[co] = len(co.atoms())
    
    if len(d) <= 1: return None
        
    for k in set_rank:
        if set_rank[k]>rank:rank,bad_set = set_rank[k],k
    
    return d[bad_set][0]


class kpropagator(object):
    def __init__(self,P):
        self._P = P
        self._num, self._denom = P.as_numer_denom()
        self._terms = self._denom.as_ordered_factors()
        self.used_chung_wu = False
        self._p = None
          
    #static measure for this type of integration
    measure = (Symbol("d")**Symbol("d")*Symbol("k_l"), (2*Symbol("\pi"))**(-1*Symbol("d")))
        
    @property
    def numerator(self): return self._num
    @property
    def denominator(self): return self._denom
    @property
    def __iter__(self): 
        for t in _terms: 
            if not c.is_constant: yield t
                
    @property
    def constant_factors(self):return [t for t in self._terms if t.is_constant()]
    
    @property
    def parametric_integral(self):return self._num/kpropagator.__alpha_parameterize_denom__(self._terms)
    
    @property
    def reduced_parametric_integral(self):
        r, self.used_chung_wu,self._p = kpropagator.__lin_reduce__(self.parametric_integral)
        return r
         
    def __alpha_parameterize_denom__(terms):
        al = lambda i : Symbol("alpha_"+str(i+1))
        if len(terms) == 1:return terms[0]       
        constants = 1
        _nus, _terms = [],[]
        for t in terms:
            if t.is_constant():constants*= t
            else:
                _nus = _nus + list(t.as_powers_dict().values())
                _terms = _terms + list(t.as_powers_dict().keys())
        nu =   reduce(add, _nus,0)
        body = reduce(add,[al(i)*_terms[i] for i in range(len(_terms))], 0)#expand in the alpha params etc. not always necessary
        
        return constants* body **nu
    
    def __assert_k2_coeff__(exp):
        _k  =Symbol("k_l")
        coeff = Poly(exp,_k).all_coeffs()
        assert len(coeff) == 3, "there is a problem. the polynomial in the propagator is not of order 2 in loop momentum"
        assert coeff[0] == 1, "the properly reduced integral should have a unity coefficient for k^2"
    
    def __lin_reduce__(par_integral):
        used_chung_wu = False
        _k  =Symbol("k_l")

        _num, _denom = par_integral.as_numer_denom()

        constant = 1
        factors = _denom.as_ordered_factors()
        if len(factors) > 1:
            _denom = factors[1]
            constant = factors[0]
            assert factors[0].is_constant(), "assumed first factor is a constant prefactor but might need to do something smarter"

        #assert has single term with power
        tup = list(_denom.as_powers_dict().items())[-1]#the single pair of; term and power but there may be an integer constant prefactor
        _p = expand(tup[0]) #expand inside so that we have just a polynommial in k
        _power = tup[1]
        #_p = _p.subs(_k, -1*_k)
        _p = collect(_p, Symbol("k_l")) #organise the poly into powers of k
        p = Poly(_p, _k)
        coeff = p.all_coeffs()
        assert len(coeff) == 3, "there is a problem. the polynomial in the propagator is not of order 2 in loop momentum"

        _denom.as_powers_dict().items()

        el = determine_elimination_variable(p)
        aps = list(__alpha_params__(_p))
        if el != None:
            print("Apply checng-wu, setting", el, "to 1")
            _p = _p.subs(el, 1)#cheng-wu says set one of the variables = 1
            coeff = Poly(_p, _k).all_coeffs()#refresh coefficients
            used_chung_wu = True
        elif len(aps) == 2:
            #TODO - there is a smarter choice to make here - we should coose the param that simplifies the expression the most
            _p = simplify(_p.subs(aps[0], 1-aps[1])) #special case simplification assuming we dont want to use them for linearisation
            coeff = Poly(_p, _k).all_coeffs()#refresh coefficients

        _p =  simplify(_p/coeff[0])
        p = Poly(_p, _k)

        coeff = p.all_coeffs() 

        if len(coeff) == 3 and coeff[1] == 0:  pass
        else:
            sub = Symbol("k_l") - (coeff[1])/2
            _p = expand(_p.subs(_k, sub))
            p = Poly(_p)

        kpropagator.__assert_k2_coeff__(_p)

        return _num/(constant*(_p**_power)), used_chung_wu, p
    
    def __reduced_terms__(exp):
        _k  =Symbol("k_l")
        _num, _denom = exp.as_numer_denom()
        constant = 1
        factors = _denom.as_ordered_factors()
        if len(factors) > 1:
            _denom = factors[1]
            constant = factors[0]
            assert factors[0].is_constant(), "assumed first factor is a constant prefactor but might need to do something smarter"
            
        _r = len(Poly(_num, _k).all_coeffs())-1    
        tup = list(_denom.as_powers_dict().items())[0]#the single pair of; term and power
        _nu = tup[1]#the power
        #return tup[0]
        coeffs = Poly(tup[0], _k).all_coeffs()
        _M = coeffs[-1]/coeffs[0]#the factor body as a polynomial and the mass from the  k^0. I normalise just in case
        return {"r": _r, "nu": _nu, "M": _M, "prefactor": constant}

    
    #write gamma functions as gamma(1) goes to n! with gamma(1) = 1
    
    def _beta_function_(a,b):
        _gamma = lambda arg:Function("\Gamma")(arg)
        return _gamma(a)*_gamma(b)/ (_gamma(a+b))
    
    def _form1(r,nu,M,prefactor):
        _gamma = lambda arg:Function("\Gamma")(arg)

        d = Symbol("d")/2
        a = r + d
        measure_pi = (2*Symbol("\pi"))**Symbol("d")
        measure_mod_with_sphere = (4*Symbol("\pi"))**d
        
        return   1/prefactor*  1/( measure_mod_with_sphere )* _gamma(a) *_gamma(nu-a)/(_gamma(d)*_gamma(nu)) * M

    #def _gamma_integral_(r, nu,M,prefactor):
    #    _gamma = lambda arg:Function("\Gamma")(arg)
    #    
    #    _d = Symbol("d")/2
    #    _a = _d - nu
    #    _sphere_vol = 2*Symbol("\pi")**(_d)/_gamma(_d)
    #    k2 = Symbol("k")**2
    #    measure_pi = (2*Symbol("\pi"))**_d
    #    mp = _d-1
    ##    #return  (1/prefactor)* _sphere_vol *  M**mp * _gamma(1-mp)*  _gamma(_d) /(_gamma(nu))
    #   return (1/prefactor) * M**alpha * 1/(4*Symbol("\pi") * _gamma(_d) ) * _beta_function_(r)

    @property 
    def reduced_terms(self):return  kpropagator.__reduced_terms__(self.reduced_parametric_integral)
        

    def _alpha_integration_(reduced_terms, elim=[], used_chunk_wu=False ):
        """
        for now im assuming one integration variable in alpha params
        """
        alpha_term = [a for a in reduced_terms["M"].atoms() if not a.is_constant() and a.name[:5]== "alpha"][0]

        from sympy import integrate
        power_of_quasi_mass = reduced_terms["r"] + Symbol("d")/2 - reduced_terms["nu"]
        kernel = reduced_terms["M"]
        factor = kernel.as_numer_denom()[1]

        for j in elim: 
            kernel = kernel.subs(j, 0)
            factor = factor.subs(j, 0)

        kernel = kernel* factor

        return  1/factor**(Symbol("d")/2) *  1/ (power_of_quasi_mass+1) *( (integrate(kernel, (alpha_term, 0, 1))))**(power_of_quasi_mass+1)


    def gamma_integral(self, eval_alpha=False,elim=[]): 
        """
        I have used one general formula but I need to test that is the appropriate one to use in all cases
        """
        terms = kpropagator.__reduced_terms__(self.reduced_parametric_integral)
        power_of_quasi_mass = terms["r"] + Symbol("d")/2 - terms["nu"]
        quasi_mass = terms["M"]**(power_of_quasi_mass)#this is the alpha power term
        if eval_alpha:
            quasi_mass = kpropagator._alpha_integration_(terms, elim=elim, used_chunk_wu=self.used_chung_wu)
        return kpropagator._form1(M=quasi_mass , r=terms["r"], nu= terms["nu"], prefactor=terms["prefactor"])
        
    
    def __repr__(self):  return self._repr_latex_()  
    def _repr_latex_(self):   return latex(self.parametric_integral,  mode='inline')
    
    

#example
#propagators = [
#        cpropagator(0, momenta= ["p", "l"], causal_type=-1),
#        cpropagator(0, momenta= ["l"], causal_type=1),
#]
#[p.display() for p in propagators]
#cpropagator.residues(propagators)[0]["integration"]

class cpropagator(object):
    def __init__(self, pid=0, power=1,causal_type=1, momenta =None,extra_mom_term=None):    
        if momenta == None: momenta = ["l"]
        if not isinstance(momenta,list): momenta = [momenta]
        ks = reduce( add, [Symbol("k_"+t) for t in momenta],0)
        #trial tool - to construct intermediates
        
        omegas = reduce( add, [Symbol("\omega_"+t) for t in momenta],0)       
        self.denom = causal_type * I* omegas + Symbol("D_"+str(pid))*ks**2 + Symbol("m_"+str(pid))
        if extra_mom_term != None: self.denom = self.denom+ extra_mom_term
        self.numer = 1
        self._pole = solve(self.denom, Symbol("\omega_l"))[0]
        self._oform = self.numer/self.denom
        #maybe do not rearrange as such just capture the effect. g is based on dividing out f' and then multipling by the coefficient factor. we read the pole "first" so to speak before rearranging
        #so we evaluate our correct g' at the pole. that means there are just two things to test. the right pole (easy) and the right g - just be careful.
        #the key observation is that in the formula, we never need to think about f, all we need is the correct g and the pole - see residual formula 
        #pole is invariant to any transformation of the function - but it so happens that the transformation of the function is always +-\omega to the pole order
        
        self.causal_prefactor = (Poly(self.denom, Symbol("\omega_l")).all_coeffs()[0]) **power
        
        self._pole_order = power
        self._causal_type = causal_type
        
    measure =  (Symbol("d\omega_l"), (2*Symbol("\pi")) **-1)
        
    def rearrange(self):
        """multiply top and bottom by I or -I and return the equivalent form and the pole"""
        #Needs to be of the form A(w-a) and f=(w-a)
        n,d = self.numer, self.denom  
        multiplier = I if d.as_coefficients_dict()[y] > 0 else -I
        self.denom = (simplify(multiplier*d))

        self.numer = multiplier*n
        self._pole  = solve(d, Symbol("\omega_l"))[0]
        self._f = Symbol("\omega_l") - self._pole
    
    def display(self):return self.value
   
    
    #todo - define equality on Pid and causality so that we can make basis sets
    @property
    def value(self): return self.numer/(self.denom**self._pole_order)
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
                res.append(reduce(add, [r["value"] for r in cpropagator.residues(pset)], 0))
            res_set.append(res)
        return res_set
    
    def integrate(pset):return simplify(I*reduce(add, [r["value"] for r in cpropagator.residues(pset)], 0))
        
       
    def residues(pset,return_smallest=True):
        from scipy.special import factorial
        W = Symbol("\omega_l")
        def _diff_(s,o):return s
        def sum_poles(st): return np.sum([p["order"] for p in st])
            
        sets = {1:[], -1:[]}
        residual_base = reduce(mul, [p.value for p in pset], 1)
        for p in pset:
            d = {"pole":p.pole} #this is the a in the formula
            d["order"] = p.pole_order
            d["residual"] = simplify( p.causal_prefactor * residual_base / p.value) #p.f not p.value
            #generally we take the mth derivitive of g and evaluate there, where m is pole order
            if p.pole_order > 1:d["residual"]= diff(d["residual"],W,n=p.pole_order)
            d["causal_type"] = p.causal_type
            # I need to check the algebra: actually I need to get the proper prefactor because in the case of g' i would need to carry the is around
            #need first to put f in the form (z-a) and then I need to reverse and put the result in the form having i\omega
            #make sure make g everything except for (z-a) so that is (prefactor) {(z-a) removed} (P2)(P3)...(Pn)
            prefactor =  int(1/factorial(p.pole_order-1,True))
            d["value"] =   prefactor * d["residual"].subs(W,p.pole)
            sets[d["causal_type"]].append(d)
        #this could become a more complex choice; maybe if the smaller set gives rise to imaginary derivs/residuals I dont want them?
        #unless by some cheeky invariance argument i can ignore the causality here - but I doubt it.
        #if they are the same size, take the one with the fewest poles - probably what i should have done anyway
        #if len(sets[1]) == len(sets[-1]): return sets[1] if sum_poles(sets[1]) < sum_poles(sets[-1]) else sets[-1]
        
        return sets[1] if len(sets[1]) < len(sets[-1]) else sets[-1]

    def from_edges(composite_graph_edges,k_node_index,enable_hack=False):
        propagators = []
        for index, d in enumerate(composite_graph_edges):
            species, causal_type, momenta =d["edge_species"], -1 if index == k_node_index else 1, ["l", "p"] if index == k_node_index else ["l"]
            if len(composite_graph_edges) > 2 and enable_hack: momenta = ["l"] #this is a hack because i dont yet know how we do it (i know how it is done in the lit)
            propagators.append(cpropagator(species, momenta= momenta, causal_type=causal_type))
        return propagators


##todo - high level wrapper around the other types
class _integral(object):   
    def __init__(self, integrand, props = {}):
        self._integrand = integrand  #- i think these are in omega at this stage and we need to treat them using cpropagator
        self._numerical_factors = []
        self._kpropagators = []     
        self._simple_form = fintegral.to_projective_space(self._kpropagators)
        #integrate out omega
        #create a collection of k terms
        #if there is only one then we can integrate - otherwise move to projective space
  
    @property
    def period(self): pass
    @property
    def uv_convergent_terms(self):  pass
  
    def to_projective_space(pset):
        #if already simple return the thing as is with now powers of alpha
        return pset 
    
    #handles
    def __iter__(self): 
        for k in self._kprops: yield k     
            
    def __call__(self):
        return self._simple_form()
        
    def __repr__(self):  return self._repr_latex_()  
    
    def _repr_latex_(self): return latex(self._integrand,  mode='inline')