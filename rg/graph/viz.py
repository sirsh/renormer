#temp: later pick your objects - just is just to prototype
from sympy import *

import operator
import functools
import math
import pandas as pd
import numpy as np
init_printing()

from . import colours

class arc(object):
    def __init__(self, angles, radius,center=None,size=None,pid=None,is_cut_line=False,colour=None):
        self.pid = pid
        self.radius = radius
        object_radius = 20 + radius
        self.center = (object_radius,object_radius) if center is None else center
        self.size=(2*object_radius,2*object_radius) if size is None else size
        self.body = ""
        self.header =  arc.get_header()
        self.colour = "black" if colour == None else colour
        
        #print(pid, ":angles for pid", angles)
        startAngle,endAngle = angles[0],angles[1]    
        mid_angle = (startAngle+endAngle)/2.
        self._is_cut_line = is_cut_line
        
        #self.body+= """<circle cx="{0}" cy="{1}" r="2" stroke="black" stroke-width="1" fill="red" /> """.format(*self.center)
            
        self.body+= self.__describeArc__(startAngle, endAngle)
        if not self._is_cut_line:
            self.body+= self.__describeArc__(startAngle, mid_angle, decorations="marker-end='url(#arrow)'")
           
    @staticmethod
    def get_header():
                 return """<defs><marker id="arrow"  markerUnits="strokeWidth" markerWidth='5' markerHeight='8' refX='0' refY='2' orient="auto">
                                      <path d="M0,0 V4 L2,2 Z" style="fill: #000000;" /> </marker></defs>"""
     
    
    def __polarToCartesian__(self, angleInDegrees, radius, center):
        centerX, centerY = center[0], center[1]
        angleInRadians = (angleInDegrees-90) *math.pi / 180.0;
        return int( centerX + (radius * math.cos(angleInRadians))), int( centerY + (radius * math.sin(angleInRadians))), angleInDegrees
          
    def __describeArc__(self, startAngle, endAngle, decorations=""):
        radius = self.radius
        wrap,sweepFlag = (0,0) if endAngle < startAngle  else (0,1)

        start = self.__polarToCartesian__(startAngle+wrap, radius, self.center) 
        end = self.__polarToCartesian__(endAngle+wrap, radius, self.center)
        largeArcFlag = "0" if endAngle - startAngle <= 180 else "1"
        #print(startAngle+wrap,endAngle+wrap, largeArcFlag, sweepFlag)
        d = [  "M", start[0], start[1], "A", radius, radius, 0, largeArcFlag, sweepFlag, end[0], end[1]  ]
        d = " ".join([str(_d)+" " for _d  in d] )   
        
        transfrm =  ""#  """transform="rotate(180, {}, {}) scale(1, -1) translate(0, -150)" """.format(*self.center) if self.change_sense else ""
        pidattr = "" if self.pid == None else """id='edge{}'""".format(self.pid)
        
        dasharray = "" if not self._is_cut_line else """stroke-dasharray='1,5' """
        
        return """<path {} d="{}" stroke="{}"  stroke-width="2px" {} {} {}/>""".format(pidattr, d,self.colour, decorations, transfrm, dasharray) 
    
    def __repr__(self):return self.header + self.body
    
    def _repr_html_(self): 
        return   """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round"> {0}  </g> </svg>""".format(self.__repr__(),  *self.center, *self.size)

#arc((130,250),radius=50,change_sense=False)

class ring_diagram(object):
    """
    This one uses a concentric circle arrangement starting from nodes connected to v_\{infty}
    This is a rough prototype to play around with the SVG
    Actually this needs a better integration with the incident matrix i.e. the incident matrix needs some concepts
    that are used in drawing such as drawing callbacks for arcs and edges mapped to a coodinate system
    then we should be able to enumerate through all arcs and edges in the diagram
    that would make it easy to show only one edge w.r.t the origin and be alot more functional
    important thing is to make sure the drawing is driven from the incidence matrix 
    this avoids the problem of connecting out to the outer rings which was just a test run
    the incident matrix then needs to be validated by independent cycles
    """
    def __init__(self, inc_matrix, options={}):
        self.x = 0
        self.y = 0
        self.margin = 10
        self.radius = 90
        center = self.radius+self.margin
        self._size=(200,200)
        self._center = (center,center)
        self._options = options
        self._inc = inc_matrix
        self._vertex_depths = [v[1] for v in list(self._inc.vertex_hierarchy())][1:]
        self._metadata = {}
        self.body = ""
        
        #for spanning tree - this is not general enough but quick fix
        self._dashed_edges = [] if not "cut_edges" in options else options["cut_edges"]
            
    @property
    def external_vertices(self):pass
     
    def __measure_inter_layer_chords__(self, existing_coords, proposed_coords):
        return 0
    
    def __polarToCartesian__(self, angleInDegrees, radius, center):
        centerX, centerY = center[0], center[1]
        angleInRadians = (angleInDegrees-90) *math.pi / 180.0;
        return int( centerX + (radius * math.cos(angleInRadians))), int( centerY + (radius * math.sin(angleInRadians))), angleInDegrees
    
    def __set_vertex_property__(self, i, key, value):
        if i not in self._metadata: self._metadata[i] = {}
        self._metadata[i][key] = value
        
    def __get_vertex_property__(self, i, key,value):
        if i in self._metadata and key in self._metadata[i]:return self._metadata[i][p]
    
    def __arrange__(self, angle, radius, vset, center=None, prev_layer_connections={}):
        center = center if center is not None else self._center
        if angle == 360: 
            self.__set_vertex_property__(vset[0], "point", center)
            return [center]#if the angle is 360 this can only be a dot in the center
        m,accepted = 1,None
        for offset in [angle]:
            proposed = [self.__polarToCartesian__( (i+1)*angle, radius,center) for i, v in enumerate(vset)]
            _m = self.__measure_inter_layer_chords__(prev_layer_connections, proposed)
            if m ==-1 or _m < m: 
                accepted, m = proposed, _m
                for i,v in enumerate(vset):
                    self.__set_vertex_property__(v, "point", [proposed[i][0],proposed[i][1]])
                    self.__set_vertex_property__(v, "angle", proposed[i][-1])
                    self.__set_vertex_property__(v, "radius", radius)
                            
        return accepted

    def edge_colour(self,pid=None):
        #this is temporary - we need an entire style format but not sure how i want to do it yet
        if pid != None:
            
            if self._inc.species_vector is not None and len(self._inc.species_vector) > pid and self._inc.species_vector[pid] < len(colours):
                return  colours[self._inc.species_vector[pid]]
                
        return "black"
    
    def __describeArc__(self, startAngle, endAngle, radius,pid=None):
        colour = self.edge_colour(pid)
        return arc((startAngle,endAngle),radius, self._center,self._size, pid=pid, is_cut_line=(pid in self._dashed_edges),colour=colour).body    
           
#     def try_get_connecting_edges(self, v):
#         #look up the metadata for v and get all edges that are adjacent to it that we know about
#         #think about us and vs here - i may be mixing them up
#         for u in self._inc.get_coincident_vertices(v):
#             if u in self._metadata and v in self._metadata: 
#                 _u,_v = u,v
#                 #if the max length is exceeded, try to get a nice arc otherwise....
#                 #max length is the length of the verticle chord connecting the bottom of one disk to the point on the outer disk          
#                 pid = self._inc.try_directed_edge_index(_v,_u)
#                 if pid == None: 
#                     _u,_v = _v,_u #flip for drawing
#                     pid = self._inc.try_directed_edge_index(_v,_u)
                    
#                 yield self.get_line([self._metadata[_u]["point"][0],
#                                      self._metadata[_u]["point"][1],
#                                      self._metadata[_v]["point"][0],
#                                      self._metadata[_v]["point"][1]],
#                                      pid=pid)
        
    def get_line(self, pts, dash_array=None,pid=None):
        
        colour = self.edge_colour(pid)
        
        #todo hash array passed in is redundant remove
        
        midpoint = [(pts[0] + pts[2])/2.,(pts[1] + pts[3])/2.]            
        #l = """<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" stroke="{4}" stroke-width:2px"  {5}   />""".format(*pts, colour,dash_array )
        pidattr = "" if pid == None else """id='edge{}'""".format(pid)
        dasharray =  "" if pid not in self._dashed_edges else """stroke-dasharray='1,7'"""
        hl =""
        l = """<path d="M{0} {1} L {2} {3}" stroke="{4}" stroke-width="2px" {5} {6}  /> """.format(*pts, colour,pidattr,dasharray )
        #print(l)
        if pid not in self._dashed_edges:#only display if not cut edge
            hl = """<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" stroke="{4}" stroke-width="2px"  {5} {6}  marker-end='url(#arrow)' />""".format(pts[0],
                                                                                                                                          pts[1], 
                                                                                                                                          midpoint[0], 
                                                                                                                                          midpoint[1], 
                                                                                                                                          colour,
                                                                                                                                          pidattr,
                                                                                                                                          dasharray )
        return l+hl
    
#     def __edge_segments__(self,a,b,angle, radius, vset):
#         #clean this up; assume that we have a vset but we dont know what the angles are
#         #we assume that they are connected in a ring (works for multiple edge)?
#         #the incidence matrix must start 1 and go -1 to check if the edge exists in that way as we have assumed a layout here
#         """
#         We always arrange vertices on disks separated by angles. Can only have vertices in [0,360]
#         We assume that the id is the arrangment id in order and that a<b for normal sense
#         if however b<a, then we have anticlockwise sense and we will change the angles for downstream reflection
#         """
#         #if there are multiple edges, then allow - otherwise restrict direction - wrapping because its a ring
#         a,b = a%len(vset),b%len(vset)
#         offset_angles = [self.__polarToCartesian__( (i+1)*angle, radius,self._center) for i, v in enumerate(vset)]
        
#         first = offset_angles[a][-1] % 360 # we should always start 
#         second = offset_angles[b][-1] 
        
#         all_con = list(self._inc.get_all_connections(a,b))
#         if len (all_con) > 0: print("mili edges detected")
        
#         pid = self._inc.try_directed_edge_index(vset[a],vset[b])
#         if pid == None:
#             a,b = b,a
#             first,second = second,first
#             pid = self._inc.try_directed_edge_index(vset[a],vset[b])

#         cap = lambda x : 0 if x == 360 else x     
#         if a< b: return (first,second), pid
#         return (first,cap(second)), pid

    def arrange_edges(self): 
        inc = self._inc
        MD = self._metadata
        def same_ring(a,b):
            if "radius" in MD[a] and "radius" in MD[b]:
                if MD[a]["radius"] == MD[b]["radius"]:return True
            return False      
        for p in MD.keys():
            for q in self._metadata.keys():
                if p >q:continue #diag
                pq_cons = list(inc.get_all_connections(p,q))
                for i, edge_meta in enumerate(pq_cons):
                    a,b = (p,q) if edge_meta[-1] == 1 else (q,p) #check convention
                    edge = [a,b]
                    edge_id = edge_meta[0]
                    
                    if i > 1 or not same_ring(a,b): #max ring links or on a different ring
                        self.body +=self.get_line([self._metadata[a]["point"][0],self._metadata[a]["point"][1],
                                       self._metadata[b]["point"][0],self._metadata[b]["point"][1]],
                                       pid=edge_id)
                    else: 
                        angles = [MD[a]["angle"],MD[b]["angle"]]
                        if i == 1:
                            angles = [angles[0]%360,angles[1]%360] 
                        #temp hack - not sure what the logic chould be but i know i should never have arcs bigger than 180
                        #for rings with 4 vertices, there might be a similar rule and i can generalise
                        if angles[0] - angles[1] >180 and angles[0] == 360:angles[0]=0
                        if angles[1] - angles[0] >180 and angles[1] == 360:angles[1]=0
                        
                        #print(edge_id,angles,i)
                        self.body += self.__describeArc__(*angles, MD[a]["radius"], pid=edge_id )

    def __display__(self):
        ##arrangement algorithm
        #each level of the hiearchy is given a radius >= 0 and maximally seperated items are added
        #there are symmetric locations and these are chosen to minimise arc lengths
        #I think maybe just try rotations of some integer factors of the freedom which is 360/|l| not including the one that gives the symmetry of course      
        marker_style = arc.get_header()
        self.body += marker_style
        levels = len(self._vertex_depths)
        
        radii_margins = [20*(levels-(i)) for i in reversed(range(len(self._vertex_depths))) ]
        angles = [(360./len(v)) for v in self._vertex_depths]

        for i, vset in enumerate(self._vertex_depths): 
            radius = self.radius-radii_margins[i]
            for coord in self.__arrange__(angles[i],radius , vset):
                self.body+= """<circle cx="{0}" cy="{1}" r="3" stroke="black" stroke-width="1" fill="black" /> """.format(*coord)
        
        self.arrange_edges()            
        if "show_labels" in self._options:
                for e in range(len(self._inc.edges)):
                    self.body+= """<text x="10%" style="fill: #000000;stroke: none; font-size: 14px;" ><textPath xlink:href="#edge{0}">e{0}</textPath></text>""".format(e)
                    
        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate(90,100,100)"> {0}  </g> </svg>""".format(self.body,  self.x, self.y, *self._size)

    #def __repr__(self):return self.__display__()           
    def __div2s__(self,n):
        while n - int(n) == 0:
            n/=2
            if n - int(n) ==0: yield n

                
    def _repr_html_(self): return self.__display__()
    
    @property
    def vertex_delta_functions(self):
        pass
    
    @property
    def loop_delta_functions(self):
        pass
                   
            
class show_spanning_trees(object):
    def __init__(self, inc_matrix, show_labels=False, max_width=None, min_split=5, alternate_background=True):
        
        total_edges = len(inc_matrix.edges)
        
        trees = inc_matrix.spanning_trees
        #todo limit the width so we do not have some crazy amount 
        l = int(np.ceil(np.sqrt(len(trees))))
        width,height = l,l
        if width < min_split: width,height=len(trees),1
        style = "" if alternate_background == True and height > 1 else "style='background:white'"
        #complement - cut the edges not in the tree
        edges_too_cut= [list(set(range(total_edges))-set(list(v))) for k, v in trees.iterrows()]
        rows = ""
        counter = 0
        for i in range(height):
            row = """<tr {}>""".format(style)
            for j in range(width):
                if counter == len(edges_too_cut):break
                o = options={"cut_edges": edges_too_cut[counter]}
                if show_labels: o["show_labels"] = True
                rd = ring_diagram(inc_matrix.copy(), options=o)
                row+="""<td>{}</td>""".format(rd.__display__() )
                counter+=1
            row+= "</tr>"
            rows+= row
                
        self.body = """<table >{}</table>""".format(rows)
        
    def _repr_html_(self): return self.body
    
   