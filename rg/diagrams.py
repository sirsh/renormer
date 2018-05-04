#todo:http://docs.sympy.org/0.7.0/modules/galgebra/latex_ex/latex_ex.html

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

class composition_diagram():
    def __init__(self, ci, compact=False, x_offset=0,y_offset=0):
        self.x = x_offset
        self.y = y_offset
        self.v,self.e = ci.graph_dataframes()
        self.default_offset = 30
        self.default_edge_x_off = 10     
        self.tensor = ci.tensor.copy()
        self.edges =ci.edges
        self.compact = compact 
        self.length = self.tensor.shape[0]
        self.__make_svg__()
    
    def partitions(self,  l, spacing=10,baseline=50):
        shift = ((l + 1) * spacing ) / 2
        for i in range(l):  yield int(((i+1)*spacing+baseline)-shift)

    @property
    def stubs(self):return pd.DataFrame(self.__stubs__).sort_values("vid")
    
    @property
    def __stubs__(self):  
        for i,k in self.e.iterrows():
            if k["internal"]:    
                edge = np.array([k["source"],k["target"]])#this can be confusing- there is directionality versus connectivity
                d = dict(k)
                d["id"], d["in"], d["vid"] = i,True, edge.max()
                yield d
                d = dict(k)
                d["id"], d["in"], d["vid"] = i,False, edge.min()
                yield d
            else:
                d = dict(k)
                source, direction = k["source"], False
                if source < 0: source, direction = k["target"], True
                d["id"], d["in"], d["vid"] = i,direction,source
                yield d

    def _assign_loop_coords(self, v,stubs):
        #pretend only 3 for now
        v["angle"] = 0
        v.loc[v["type"]=="source", "angle"] = 90
        v.loc[v["type"]=="sink", "angle"] = 270
        v["loop_coord"] = v["angle"].apply(composition_diagram._polarToCartesian)
        
        all_angles = v["angle"]
        #update the internal edges to say what "angle" they start and end with
        stubs["sangle"] = 0
        stubs["eangle"] = 0
        for index in v.index:
            stubs.loc[stubs["source"] == index, "sangle"] = all_angles[index]
            stubs.loc[stubs["target"] == index, "eangle"] = all_angles[index]
        
    @property
    def coords(self):
        v,e = self.v,self.e
        major_offset, minor_spacing,baseline = 50,18,50 #todo make these class props?
        centers = list(reversed([20+i*major_offset for i in range(len(v))]))#from left to right draw backwards
        v["x"], v["y"]  = centers, baseline    
        #add all the stubs
        newVs = []
        for i,k in v.iterrows():
            internals = k["degree_in"]
            externals = k["degree_out"]
            #todo add species to these, do we know if they are internal or not, can we draw edges - well 
            for p in self.partitions(internals, spacing=minor_spacing, baseline=baseline):
                newVs.append({ "x": centers[i]+minor_spacing, "y": p, "in": True, "vid":i})
            for p in self.partitions(externals, spacing=minor_spacing, baseline=baseline):
                newVs.append({ "x": centers[i]-minor_spacing, "y": p, "in": False, "vid":i})
        stubs = pd.DataFrame(newVs)
        
        #NEED TO BE CAREFUL HERE - I AM NOT SURE ABOUT THE ORDER OF THIS THING IN GENERAL - we should join the rank carefully
        stubs["rank"]= stubs.reset_index().groupby(["vid", "in"]).rank(method="min")["index"].astype(int)
        estubs = self.stubs
        estubs["rank"]= estubs.groupby(["vid", "in"]).rank(method="min")["id"].astype(int)
        stubs = pd.merge(estubs,stubs, on=["in", "vid", "rank"], how='inner') 
        assert len(estubs) == len(stubs), "there should be a coordinate for every stub entry which is not the case!"
        assert len(estubs) == len(stubs.groupby(["x", "y"])), "although there is a coordainte for each stub, they are not unique - they should be!"
        self._assign_loop_coords(v,stubs)
        return v, stubs, None

    @property
    def body(self):  return self._body
        
    def __repr__(self):
        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="500" height="120">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round"> {0}  </g> </svg>""".format(self.body, self.x, self.y)
    
    def _repr_html_(self): return str(self) if not self.compact else self.display_loop()
    
    def symmetric_offsets(self, l, offset=70): return [offset-20+20*i for i in range(l)]
    
    def draw_line(self, pts, colour, dash_array=None):
        dash_array = """stroke-dasharray='{}'""".format(dash_array) if dash_array != None else ""
        self._body+= """<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" stroke="{4}" stroke-width:2"  {5} />""".format(*pts, colour,dash_array )
        
    def _polarToCartesian(angleInDegrees, centerX=100, centerY=100, radius=50):
        angleInRadians = (angleInDegrees-90) *math.pi / 180.0;
        return {  "x": centerX + (radius * math.cos(angleInRadians)), "y": centerY + (radius * math.sin(angleInRadians))  };

    def _describeArc(startAngle, endAngle, x=100, y=100, radius=50):
        #dirty HAck because i dont want think about the general case yet
        if startAngle == 90 and endAngle == 0: startAngle, endAngle = endAngle, startAngle
        start = composition_diagram._polarToCartesian(endAngle, x, y, radius);
        end = composition_diagram._polarToCartesian(startAngle, x, y, radius );
        largeArcFlag = "0" if endAngle - startAngle <= 180 else "1"
        d = [  "M", start["x"], start["y"], "A", radius, radius, 0, largeArcFlag, 0, end["x"], end["y"]  ]
        d = " ".join([str(_d)+" " for _d  in d] )
        return d;       

    #todo - i should do the body in relative coords so that i can place lots of diagrams in the picture
    def __make_svg__(self):
        self._body = ""
        colours = ["green", "orange", "red"]
        vertices,stubs,links = self.coords
        vertices_coords = []
        #the original vertex cores are drawn on the horizontal
        for i,k in vertices.iterrows():
            vertices_coords.append([k["x"], k["y"]])
            self._body +="""<circle cx="{0}" cy="{1}" r="3" stroke="black" stroke-width="1" fill="black" /> """.\
            format(k["x"], k["y"])
        #we draw stubs from the original vertex which are connection points
        for i,k in stubs.iterrows():
            colour = colours[k["species"]]
            self._body +="""<circle cx="{0}" cy="{1}" r="2" stroke="{3}" stroke-width="1" {2} /> """.\
            format(k["x"], k["y"], "" if k["in"] else "fill='{}'".format(colour), colour)
        #draw edges - the edge id is the grouping - internal edges should be paired up neatly
        for i,grp in stubs.groupby("id"):
            l = len(grp)
            dash_array = "2,2" if l == 2 else None
            head = dict(grp.iloc[0])
            self.draw_line([head["x"], head["y"], *vertices_coords[head["vid"]]],colours[head["species"]], dash_array=dash_array)
            if l == 2: #internal edges
                tail = dict(grp.iloc[-1])
                self.draw_line([head["x"], head["y"], tail["x"], tail["y"]],colours[head["species"]], dash_array=dash_array)
                self.draw_line([tail["x"], tail["y"], *vertices_coords[tail["vid"]]],colours[tail["species"]], dash_array=dash_array)
        return self._body
  
    def display_loop(self):       
        colours = ["green", "orange", "red"]
        _body = ""
        iedges = self.coords[1]

        for k,v in self.coords[0].iterrows():
            _body +="""<circle cx="{0}" cy="{1}" r="3" stroke="black" stroke-width="1" fill="black" /> """.\
            format(v["loop_coord"]["x"], v["loop_coord"]["y"])
            _body+= diagram._residual_svg_(v["loop_coord"]["x"], v["loop_coord"]["y"], self.tensor[k],True)
            #append on residuals at each vertex    
        iedges = iedges[(iedges["internal"]==True)&(iedges["in"]==False)]       
        for k,e in iedges.iterrows():
            _body += """<path d="{}" stroke="{}"/>""".format(composition_diagram._describeArc(e["sangle"],e["eangle"]),colours[e["species"]])

        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="200" height="180">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round"> {0}  </g> </svg>""".format(_body, self.x, self.y)
  
    
#some way to control the scale of the sub diagrams e.g. the loops (and maybe other parameters being passed through as dicts
class diagram_set(object):
    def __init__(self, collection,offset=50,compact=False):
        DG = composition_diagram
        total_offsets = offset * len(collection) + offset
        entities = None
        if compact:  entities = [DG(c, y_offset=offset*(i)).display_loop() for i, c in enumerate(collection)]
        else:  entities = [DG(c, y_offset=offset*(i)).__repr__() for i, c in enumerate(collection)]     
        self.body =  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="0" y="0" width="500" height="{1}">{0}</svg>""".format(entities,total_offsets+10)
    def __repr__(self):  return  self.body
    def _repr_html_(self):   return str(self)     
    
    #todo - fluent - maybe columns for distinct class and then fill the representatives in rows
    def by_residual(): return self
    def by_residual_complement(): return self
        
        
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
    
    def _repr_html_(self):  return str(self)
    
    def _residual_svg_(x,y, residual_tensor, is_external=False, colour_map=["green", "orange", "red"]):
        species_dict = lambda i : { "stroke-dasharray": None, "line" : "straight", "stroke": colour_map[i] }
        ins_, outs_ = list(residual_tensor[:,0]), list(residual_tensor[:,1])
        if (len(ins_)+len(outs_)) == 0:return ""
        ispec,ospec,body = [], [], ""
        for spec, power in enumerate(ins_):  ispec += [spec for i in range(power)]
        for spec, power in enumerate(outs_):  ospec += [spec for i in range(power)]
        if is_external == False:#temp: I want to refactor this stuff but here i pretend as though there is an extra but write back an existing
            ins_.append(ins_[-1])
            outs_.append(outs_[-1])
        ispec = diagram.sort_species(ispec)
        iangles = diagram.angles(ispec)
        ofpsec = diagram.sort_species(ospec)
        oangles = diagram.angles(ospec,is_out=True)
        for c,i in enumerate(ispec): body+= _leg(species_dict(i)).set_location(x,y,iangles[c]).__repr__()
        for c,i in enumerate(ospec): body+= _leg(species_dict(i)).set_location(x,y,oangles[c]).__repr__()

        return body

        
    
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