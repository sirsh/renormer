import numpy as np
import pandas as pd
from IPython.display import SVG
from IPython.display import HTML
from .interaction import interaction
#todo -create a function to control the curvey lines so that we can do things with them - control the amplitudes - down to 0, colours, angles (transformed around diagram center) - long term. keep ugly for now
#make sure all the 2/3 species that we deal with can be rendered using two styles (red wavy, green straight)

species_styles = [
    { "stroke": "red",  "stroke-dasharray": None, "line" : "straight" },
    { "stroke": "black",  "stroke-dasharray":None, "line" : "wavy" }    
]

import numpy as np
import itertools
import math



class diagram_complex(HTML):
    def __init__(self, diagrams, layout_props={}):
        pass
        
    def _repr_html_(self): 
        #create a layout for all the diagrams and display them in a table
        return test

class diagram(HTML):
    def __init__(self, species_matrix, props={}, loop_radius=-1):
        if isinstance(species_matrix, interaction):species_matrix = species_matrix._mat
        arr = np.array(species_matrix)
        self.center = (80,50)
        self.legs = []
        self.loop_radius = 2
        if "loop_radius" in props: self.loop_radius = props["loop_radius"]
        
        ins_, outs_ = list(arr[:,0]), list(arr[:,1])
        ispec,ospec = [], []
        for spec, power in enumerate(ins_):  ispec += [spec for i in range(power)]
        for spec, power in enumerate(outs_):  ospec += [spec for i in range(power)]
            
        ispec = diagram.sort_species(ispec)
        iangles = diagram.angles(ispec)
        ofpsec = diagram.sort_species(ospec)
        oangles = diagram.angles(ospec,is_out=True)
        
        for c,i in enumerate(ispec):
            a = iangles[c]
            self.attach_leg(_leg(species_styles[i]),a)
            
        for c,i in enumerate(ospec):
            a = oangles[c]
            self.attach_leg(_leg(species_styles[i]),a)
  
    def inter_distance(l):
        score = 0
        ar = np.array(l)
        for u in np.unique(l):
            L = list(np.nonzero(ar==u)[0])
            L = list(itertools.combinations(L,2))
            Ds = [abs(t[0]-t[1]) for t in L]
            score +=np.array(Ds).sum()
        return score

    def sort_species(l):
        P = list(itertools.permutations(l))  
        return list(sorted(P, key=diagram.inter_distance)[0])

    def angles(l,is_out=False):
        offset = 0 if is_out == False else -180
        divs = 180 / (len(l) + 1)
        return [((i+1)*divs)-90+offset for i in range (len(l)) ]
    
    def attach_leg(self, leg, theta): 
        self.legs.append(leg.set_location(self.center[0],self.center[1], theta))
        
        return self
        
    def _leg_repr_(self):
        ls = ""
        for l in self.legs: ls += " {} ".format(str(l))
        return ls
    
    def __repr__(self):
        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="0" y="0" width="240" height="120">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round">
        {2}
      </g>
      <circle cx="{0}" cy="{1}" r="{3}" stroke="black" stroke-width="1" fill="black" />
    </svg>""".format(self.center[0], self.center[1], self._leg_repr_(),self.loop_radius)
    
    def _repr_html_(self): 
        return str(self)
        
    
#colour
#waviness
#break style
#label
class _leg(dict):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        self.x = 0
        self.y= 0
        self.rotation = 0
        
    def set_colour(self, col):
        self["stroke"] = col
        return self
    
    def set_location(self, x,y, rotation):
        self.x = x
        self.y= y
        self.rotation = rotation
        return self
        
    def __transform(self): 
        return " translate({0},{1}) rotate({2} 0, 0 )".format(self.x, self.y, self.rotation)
    
    def __repr__(self):
        if "line" in self and self["line"] == "wavy":
            return """<path d="M 0,0 C 2,-3.142 3,-5 5,-5 S 8,-3.142 10,0 S 13,5 15,5 S 18,3.142 20,0 C 22,-3.142 23,-5 25,-5 S 28,-3.142 30,0 S 33,5 35,5 S 38,3.142 40,0" 
                       stroke="{}" stroke-dasharray="{}" stroke-width:2" transform="{}" />""".format(self["stroke"],self["stroke-dasharray"], self.__transform())
        else: return """<line x1="0" y1="0" x2="40" y2="0" stroke="{}" stroke-dasharray="{}" stroke-width:2" transform="{}" />""".format(self["stroke"],self["stroke-dasharray"],self.__transform())
            
    

#mat = [[2,1],[1,1]]
#d = diagram(mat)#.\
#str(d)
#d
#presumably there are some 'equiavlent' permutations for the outs and we need to bias them based on the ins; todo!