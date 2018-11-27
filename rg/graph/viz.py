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
    def __init__(self, angles, radius,center=None,size=None,pid=None,is_cut_line=False,colour=None, is_directed=True):
        self.pid = pid
        self.radius = radius
        object_radius = 20 + radius
        self.center = (object_radius,object_radius) if center is None else center
        self.size=(2*object_radius,2*object_radius) if size is None else size
        self.body = ""
        self.header =  arc.get_header()
        self.colour = "black" if colour == None else colour
        self.is_directed = is_directed
        
        #print(pid, ":angles for pid", angles)
        startAngle,endAngle = angles[0],angles[1]    
        mid_angle = (startAngle+endAngle)/2.
        self._is_cut_line = is_cut_line
        
        #self.body+= """<circle cx="{0}" cy="{1}" r="2" stroke="black" stroke-width="1" fill="red" /> """.format(*self.center)
            
        self.body+= self.__describeArc__(startAngle, endAngle)
        if not self._is_cut_line and is_directed:
            arrow = "arrow" if self.colour not in ["grey", "orange"] else "arrow"+self.colour
         
            self.body+= self.__describeArc__(startAngle, mid_angle, decorations="marker-end='url(#{})'".format(arrow))
           
    @staticmethod
    def get_header():
        #temp - need to redo this like so much else
                 return """<defs><marker id="arrow"  markerUnits="strokeWidth" markerWidth='5' markerHeight='8' refX='0' refY='2' orient="auto" stroke="black" >
                                      <path d="M0,0 V4 L2,2 Z" style="fill: #000000;" /> </marker>
                                      <marker id="arrowgrey"  markerUnits="strokeWidth" markerWidth='5' markerHeight='8' refX='0' refY='2' orient="auto" stroke="grey" >
                                      <path d="M0,0 V4 L2,2 Z" style="fill: #000000;" /> </marker>
                                      <marker id="arroworange"  markerUnits="strokeWidth" markerWidth='5' markerHeight='8' refX='0' refY='2' orient="auto" stroke="orange" >
                                      <path d="M0,0 V4 L2,2 Z" style="fill: #000000;" /> </marker>
                                      </defs>"""
        
     
    
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
    
    @staticmethod
    def __describeArc_PTS__(start, end, radius=50, colour="black", pid=None,_is_cut_line=True):

        wrap,sweepFlag = (0,0) #if endAngle < startAngle  else (0,1)
        largeArcFlag = "0" #if endAngle - startAngle <= 180 else "1"
        #print(startAngle+wrap,endAngle+wrap, largeArcFlag, sweepFlag)
        d = [  "M", start[0], start[1], "A", radius, radius, 0, largeArcFlag, sweepFlag, end[0], end[1]  ]
        d = " ".join([str(_d)+" " for _d  in d] )   
        transfrm =  ""#  """transform="rotate(180, {}, {}) scale(1, -1) translate(0, -150)" """.format(*self.center) if self.change_sense else ""
        pidattr = "" if pid == None else """id='edge{}'""".format(pid)
        
        dasharray = "" if not _is_cut_line else """stroke-dasharray='1,5' """
        
        return """<path {} d="{}" stroke="{}"  stroke-width="2px" {}{}/>""".format(pidattr, d,colour, transfrm, dasharray) 
    
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
        self._is_directed = True if "is_directed" not in options else options["is_directed"]
        
        #for spanning tree - this is not general enough but quick fix
        self._dashed_edges = [] if not "cut_edges" in options else options["cut_edges"]
        
        if self._inc.shape[0] == 2:# if there is only one vertex then default to show ex edges
            self._options["show_ex_edges"] = True
            
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
            proposed_text = [self.__polarToCartesian__( (i+1)*angle, radius+10,center) for i, v in enumerate(vset)]
            _m = self.__measure_inter_layer_chords__(prev_layer_connections, proposed)
            if m ==-1 or _m < m: 
                accepted, m = proposed, _m
                for i,v in enumerate(vset):

                    self.__set_vertex_property__(v, "point", [proposed[i][0],proposed[i][1]])
                    self.__set_vertex_property__(v, "text_point", [proposed_text[i][0],proposed_text[i][1]])
                    
                    self.__set_vertex_property__(v, "angle", proposed[i][-1])
                    self.__set_vertex_property__(v, "radius", radius)
                            
        return accepted

    def edge_colour(self,pid=None):
        #this is temporary - we need an entire style format but not sure how i want to do it yet
        if pid != None:
            try:
                if self._inc.species_vector is not None and len(self._inc.species_vector) > pid and self._inc.species_vector[pid] < len(colours):
                    return  colours[self._inc.species_vector[pid]]
            except Exception as ex:
                print(repr(ex))
                return "black"
        return "black"
    
    def __describeArc__(self, startAngle, endAngle, radius,pid=None):
        colour = self.edge_colour(pid)
        return arc((startAngle,endAngle),radius, self._center,self._size, pid=pid, 
                   is_cut_line=(pid in self._dashed_edges),colour=colour,is_directed=self._is_directed).body    
  
        
    def get_line(self, pts, dash_array="",pid=None, colour=None):
        
        colour = self.edge_colour(pid) if colour == None else colour
        
        #todo hash array passed in is redundant remove
        
        midpoint = [(pts[0] + pts[2])/2.,(pts[1] + pts[3])/2.]            
        #l = """<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" stroke="{4}" stroke-width:2px"  {5}   />""".format(*pts, colour,dash_array )
        pidattr = "" if pid == None else """id='edge{}'""".format(pid)
        dasharray =  dash_array if pid not in self._dashed_edges else """stroke-dasharray='1,7'"""
        hl =""
        l = """<path d="M{0} {1} L {2} {3}" stroke="{4}" stroke-width="2px" {5} {6}  /> """.format(*pts, colour,pidattr,dasharray )
        #print(l)
        if pid not in self._dashed_edges and self._is_directed:#only display if not cut edge and it is a directed edge
            arrow = "arrow" if colour not in ["grey", "orange"] else "arrow"+colour
         
            hl = """<line x1="{0}" y1="{1}" x2="{2}" y2="{3}" stroke="{4}" stroke-width="2px"  {5} {6}  marker-end='url(#{7})' />""".format(pts[0],
                                                                                                                                          pts[1], 
                                                                                                                                          midpoint[0], 
                                                                                                                                          midpoint[1], 
                                                                                                                                          colour,
                                                                                                                                          pidattr,
                                                                                                                                          dasharray,
                                                                                                                                          arrow)
        return l+hl
    

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
                #e.g. a=0 and b=1 are connected like (0,-1) if b->a
                pq_cons = list(inc.get_all_connections(p,q))
                for i, edge_meta in enumerate(pq_cons):
                    #-1 should be EXITING a
                    a,b = (p,q) if edge_meta[-1] == -1 else (q,p) #check convention - this is only from the get_all_connections function
                    edge = [a,b]
                    edge_id = edge_meta[0]
                    
                    #print("i a,b, is some ring", i ,a, b, same_ring(a,b))
                    
                    if i > 1 or not same_ring(a,b): #max ring links or on a different ring
                        self.body +=self.get_line([self._metadata[a]["point"][0],self._metadata[a]["point"][1],
                                       self._metadata[b]["point"][0],self._metadata[b]["point"][1]],
                                       pid=edge_id)
                        
                        #print("direct line for ", a, b)
                    else: 
                        angles = [MD[a]["angle"],MD[b]["angle"]]
                        if i == 1:
                            angles = [angles[0]%360,angles[1]%360] 
                        #temp hack - not sure what the logic chould be but i know i should never have arcs bigger than 180
                        #for rings with 4 vertices, there might be a similar rule and i can generalise
                        if angles[0] - angles[1] >180 and angles[0] == 360:angles[0]=0
                        if angles[1] - angles[0] >180 and angles[1] == 360:angles[1]=0
                        
                        #print("edge and anggle for source", edge_id,angles,i)
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
            for vid, coord in enumerate(self.__arrange__(angles[i],radius , vset)):
                vtext_coords = self._metadata[vid]["text_point"] if "text_point" in self._metadata[vid] else None
                self.body+= """<circle cx="{0}" cy="{1}" r="3" stroke="black" stroke-width="1" fill="black"  /> """.format(*coord)
                ##add vertext text - optionally
                #self.body+= """ <text x="{0}" y="{1}" style="fill: #000000;stroke: none; font-size: 14px;" transform="rotate(-90,{0},{1}) translate(-7 5)">v{2}</text>""".format(*vtext_coords,vid)
        
                #get the star for this vertex and add all of the edges in the same arc without arrows for now
                #S = self._inc.get_star(vid)
                #ex_edges = S[S[:,-1]==-1][:,0]
                #num_ex = len(ex_edges)
                #sangle = angles[i] - (180-angles[i])
                
                #for j, e in enumerate(ex_edges):
                #    offset_angle = j*(90./(num_ex+1))
                #    pt1 = (coord[0], coord[1])
                #    pt = self.__polarToCartesian__(offset_angle+sangle, 30, coord)[:2]
                #    self.body +=self.get_line(pt1+pt) #pid=e
                #go to the point at an angle
                
        
        self.arrange_edges()    
        
        if "show_labels" in self._options:
                for e in range(len(self._inc.edges)):
                    self.body+= """<text x="10%" style="fill: #000000;stroke: none; font-size: 14px;" ><textPath xlink:href="#edge{0}">e{0}</textPath></text>""".format(e)
          
        if "show_ex_edges" in self._options: 
            #self.arrange_external_edges()
            #take the id rows from the residual - inc matrix knows what species these are
            entering = self._inc.residual(exiting=False,sort=True)[:,0]
            exiting =  self._inc.residual(entering=False,sort=True)[:,0]
            for i, e in enumerate(entering):
                angle = i*(180./(len(entering)+1))
                pt = self.__polarToCartesian__(angle, 50,self._center)[:2]
                self.body +=self.get_line(pt+self._center,pid=e)        
                #print("entering", pt, "ang", angle)
            for i, e in enumerate(exiting):
                angle = 180+i*(180./(len(exiting)+1))
                pt = self.__polarToCartesian__(angle, 50,self._center)[:2]
                self.body +=self.get_line(self._center+pt,pid=e)
                #print("exiting", pt, "ang", angle)
         
     
        rotation = 90
        if "add_rotation" in self._options:rotation+=self._options["add_rotation"]
        
        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate({5},100,100)"> {0}  </g> </svg>""".format(self.body,  self.x, self.y, *self._size,rotation)

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
    def __init__(self, inc_matrix, show_labels=False, max_width=None, min_split=5, alternate_background=True, is_directed=True):
        
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
                o = options={"cut_edges": edges_too_cut[counter], "is_directed" : is_directed}
                if show_labels: o["show_labels"] = True
                rd = ring_diagram(inc_matrix, options=o)
                row+="""<td>{}</td>""".format(rd.__display__() )
                counter+=1
            row+= "</tr>"
            rows+= row
                
        self.body = """<table >{}</table>""".format(rows)
        
    def _repr_html_(self): return self.body
    
class star(ring_diagram):
    def __init__(self, inc, options={}):
        self.root = (0,0)
        self._inc = inc
        self.body = ""
        self.radius = 35
        self._dashed_edges = []
        self._all_points = [[0,0]]
        self._is_directed = True
        self._vertex_depths = [0]
        self.x = 0
        self.y = 0
        self._size=(10,10) if "size" not in options else options["size"]
        
        vroot = inc.root_vertex()
        
        self.star_svg(vroot, center=(0,0))
        self._options= options
        
        
    def star_svg(self,vid,center,pid=None):
        marker_style = arc.get_header()
        self.body += marker_style
        
        v_pseudo_residual = self._inc.get_star(vid)
        entering = v_pseudo_residual[v_pseudo_residual[:,1]==1]
        exiting =  v_pseudo_residual[v_pseudo_residual[:,1]==-1]
        #sort so that we do internal links first
        exiting = exiting[exiting[:,-1].argsort()][::-1]
        
        for i, e in enumerate(entering[:,0]):
            angle = i*(180./(len(entering)+1))
            pt = self.__polarToCartesian__(angle, self.radius,center)[:2]
            self.body +=self.get_line(pt+center,pid=e)        
            #print("entering", pt, "ang", angle)
            self.update_bounds(pt)

        for i, e in enumerate(exiting[:,0]):
            pie_size = 180./(len(exiting)+1)
            angle = 180+((i)*pie_size)
            pt = self.__polarToCartesian__(angle, self.radius,center)[:2]    
            self.update_bounds(pt)
            self.body +=self.get_line(center+pt,pid=e)
            #if there is an internal v at the of e, then render that do
            ch = exiting[i][-1]
            if ch != -1:
                cpt = self.__polarToCartesian__(angle, 2*self.radius,center)[:2]
                self.update_bounds(cpt)
                self.star_svg(ch, center=cpt)
                
    #todo for dynamic diagrams, we should determine the size and center: center used for rotos
    def update_bounds(self,pt): 
        self._all_points.append(pt)
    
    #[0,15,20, 45, (60), 90] -> [0,15,20, 45, (58), 90]
    #[0,2,4,pid,3,1] - > [1,0,2,4,pid,3]
    def __shift_to__(l,p,i):
        _i = l.index(p)
        shift =  (i - _i) % len(l)
        return l[-shift:]+l[:-shift] 

    #def __get_star_children__(S):  return S[S[:,1]==-1]
        
    def __display__(self):
        padding = 50
        pts = np.array(self._all_points)
        x1,x2 = pts[:,0].min(),  pts[:,0].max()
        y1,y2 = pts[:,1].min(),  pts[:,1].max()

        self._size = (abs(x1-x2)+padding, abs(y1-y2)+padding)        
        self._center = (self._size[0]/2. ,self._size[1]/2. -padding/2.)
        trans = "translate({} {})".format(self._size[0]/2.,padding)
        rotation = 90 if "extra_rotation" not in self._options else 90 + self._options["extra_rotation"]
        rotation = 0 # test abs
        
        fstring = """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate({5},{6},{7}) {8} "> {0}  </g> </svg>"""
        return fstring.format(self.body,  self.x, self.y,*self._size, rotation, *self._center, trans)
    
        return self.body


# i think the debug grahp is fairly flexi. For one loops we can sort in a ring. for trees we can sort in a spine maybe in terms of out end points or bifuricate
class debug_graph(ring_diagram):
    def __init__(self, inc, options={}):
        self.root = (0,0)
        self._inc = inc
        self.body = ""
        self.radius = 35
        self._dashed_edges = [-1]
        self._all_points = [[0,0]]
        self._is_directed = True
        self._vertex_depths = [0]
        self.x = 0
        self.y = 0
        self._size=(500,250) if "size" not in options else options["size"]
        self._options= options
        self.rotation = 0
        self._center= (self._size[0]/2.,self._size[1]/2.)
        
        self.update_body()
        
    def update_body(self):
        marker_style = arc.get_header()
        self.body += marker_style
        vw = 100
        vinfty = self._inc.shape[0] - 1
        internals = vinfty - 1
        vs = {}
        vs_lu = {}
        w = vinfty * 100
        
        if "hide_vinfty" not in self._options:
            self.body+= """<circle cx="{0}" cy="{1}" r="2.5" stroke="grey" stroke-width="2" fill="black" /> """.format(w/2., 30)
                        
        offset = vw/2.
        for i, v in enumerate(self._inc[:-1]):
            x = offset + vw*i
            y = 100
            vs[i] = [x,y]
            self.body+= """<circle cx="{0}" cy="{1}" r="2.5" stroke="black" stroke-width="2" fill="black" /> """.format(x,y)
            
            vs_lu[i] = {}
            
            center = tuple(vs[i])
            v_pseudo_residual = self._inc.get_star(i)
            entering = v_pseudo_residual[v_pseudo_residual[:,1]==1]
            exiting =  v_pseudo_residual[v_pseudo_residual[:,1]==-1]

            for j, e in enumerate(entering[:,0]):
                pie_size = 180./(len(entering)+1)
                angle = 90+((j)*pie_size)
                pt = self.__polarToCartesian__(angle, self.radius,center)[:2]
                vs_lu[i][e] = pt
                self.body +=self.get_line(pt+center,pid=e)        

            for j, e in enumerate(exiting[:,0]):
                pie_size = 180./(len(exiting)+1)
                angle = 270+((j)*pie_size)
                pt = self.__polarToCartesian__(angle, self.radius,center)[:2]  
                vs_lu[i][e] = pt
                self.body +=self.get_line(center+pt,pid=e)

        #for all the edges create a midway point between the two vs and a directed circuit accorss them using dotted lines?
        for i, e in enumerate(self._inc.T):
            #add edge coords and for each vertex, for each edge, say where it ends
            try:
                ed = np.nonzero(e)[0]
                edge = e[ed]
                a,b = tuple(ed) if edge[0] == -1 else tuple(reversed(ed))        
                pa,pb = vs_lu[a][i],vs_lu[b][i] 

                if a < b: #do arc instead
                    ar = arc.__describeArc_PTS__(pa,pb,pid=-1,colour="grey")
                    self.body += ar
                else:
                    self.body +=self.get_line(( *pa,*pb ) ,colour="grey", dash_array="""stroke-dasharray='1,7'""" )
            except Exception as ex:
                #print(repr(ex))      
                pass
    def __display__(self):
        fstring = """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate({5},{6},{7}) "> {0}  </g> </svg>"""
        return fstring.format(self.body,  self.x, self.y, *self._size, self.rotation, *self._center)
        return self.body

class simple_ring_graph(ring_diagram):
    """
    example
    B = incidence_matrix(edges= [ [-1,3], [-1,3],[3,0],[3,1],[-1,1],[1,2],[2,0],[0,-1]],  species_vector=[0,1,0,1,0,0,0,0])
    """
    def __init__(self,inc, options={}):
        self.body=""
        self._dashed_edges = []
        self._is_directed=True
        self.sz = inc.shape[0] - 1
        self._center = (100,100)
        self._size = (200,200)
        self._radius=50
        self.locations = []
        self._inc = inc

        locs = [0,180]
        if self.sz == 3:locs = [0,90,180]
        if self.sz == 4:locs = [0,45,135,180]
        
        for v in range(self.sz):
            s= inc.get_star(v)
            ex,ex_dirs = [],[]
            e = -1
            try:
                e = s[np.where(s[:,1]==-1)[0]][0][0]
            except:
                #if this happens it is because it is not a proper ring but I need to fix how the drawing is done
                e = s[np.where(s[:,1]==1)[0]][0][0]
                
            if len(np.where(s[:,-1]==-1)[0]) > 0:#if has externals added the data
                ex = list(s[np.where(s[:,-1]==-1)[0]][:,0])
                ex_dirs = list(s[np.where(s[:,-1]==-1)[0]][:,1])   
            self.locations.append([locs[v], e, ex,ex_dirs])
                                
            #print(self.locations)
          
    def star_at(self, angle, out_edges=[], out_dirs = []):
        pt= self.__polarToCartesian__(angle, self._radius, self._center) [:2] 
        self.body+= """<circle cx="{0}" cy="{1}" r="2.5" stroke="black" stroke-width="2" fill="black" /> """.format(*pt)     
        if len(out_edges) == 0:return
        angles = [angle]
        if len(out_edges) == 2: angles = [angle-20, angle+20]
        if len(out_edges) == 3: angles = [angle-35, angle, angle+35]     
        #print(angles,out_dirs)
        for i, a in enumerate(angles):     
            pta = self.__polarToCartesian__(a, 30, pt)[:2]
            if out_dirs[i] == -1: self.body +=self.get_line((*pt, *pta ), pid=out_edges[i] )
            else: self.body +=self.get_line((*pta, *pt ), pid=out_edges[i] )
                
    def __display__(self, ):            
        marker_style = arc.get_header()
        self.body += marker_style
        loc = self.locations
        for i, l in enumerate(loc):
            nexti = (i + 1) % len(loc)
            end = loc[nexti][0] if loc[nexti][0] != 0 else 360
            if end == 360: 
                self.body += self.__describeArc__( l[0],end, radius=self._radius, pid=l[1])
            else:  self.body += self.__describeArc__(end,l[0],  radius=self._radius, pid=l[1])
            self.star_at(l[0], l[-2], l[-1])
        
        rotation = -90
        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
          <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate({5},100,100)"> {0}  </g> </svg>""".format(self.body,  10, 10, 200,200 ,rotation)
    
    def _repr_html_(self): return self.__display__()
    
    
class tabulate_graphs(object):

    def __init__(self, collection,  max_width=None, min_split=3, max_split=7, alternate_background=False, diagram_class=debug_graph, options={}):
        l = int(np.ceil(np.sqrt(len(collection))))
        width,height = l,l
        if width < min_split: width,height=len(collection),1
        if width > max_split: width,height = max_split,  int(np.ceil(len(collection)/max_split))
        style = "" if alternate_background == True and height > 1 else "style='background:white'"
        
        rows = ""
        counter = 0
        for i in range(height):
            row = """<tr {}>""".format(style)
            for j in range(width):
                if counter < len(collection):
                    cell_body = "ERROR AT " +str( counter )
                    try:
                        rd = diagram_class(collection[counter], options=options)
                        cell_body = rd.__display__()
                    except: pass
                    row+="""<td>{}</td>""".format(cell_body )
                counter+=1
            row+= "</tr>"
            rows+= row
                
        self.body = """<table >{}</table>""".format(rows)
        
    def _repr_html_(self): return self.body

    
    
#class melon_graph(ring_diagram):
#    def __init__(self,left_angles = [], right_angles =[], mids = []):
#        self.x = 10
#        self.y =10
#        self._size = (200,200)
#        self._dashed_edges = [0]
#        self.body = ""
#        self._is_directed = False
#        self.left_angles = left_angles
#       self.right_angles= right_angles
#       self.mids = mids
##        
#    def __display__(self, ):      
#        self.body += arc((135,45),50, (10,100) ,self._size,  is_cut_line=True,colour="black",is_directed=False).body    
#        self.body += arc((225,320),50, (80,100) ,self._size,  is_cut_line=True,colour="black",is_directed=False).body          
#        self.body +=self.get_line((30,64, 60,64), pid=1 )
#        self.body +=self.get_line((44,70, 45,140), pid=0 )
#        self.body+= """<circle cx="{0}" cy="{1}" r="3" stroke="black" stroke-width="1" fill="black"  /> """.format(44,64)        
#        proposed = self.__polarToCartesian__( 135, 50,(10,100)) 
#        self.body+= """<circle cx="{0}" cy="{1}" r="3" stroke="black" stroke-width="1" fill="black"  /> """.format(*proposed)
#                
#        for a in self.right_angles:
#            if a < 0:  self.body+= """<circle cx="{0}" cy="{1}" r="6" stroke="red" stroke-width="2" fill="none" stroke-dasharray='1,3' /> """.format(*self.__polarToCartesian__( abs(a), 50,(10,100)) #
#            self.body+= """<circle cx="{0}" cy="{1}" r="2" stroke="black" stroke-width="1" fill="black"  /> """.format(*self.__polarToCartesian__( abs(a), 50,(10,100)) )
#        
#        for a in self.left_angles:
#            if a < 0: self.body+= """<circle cx="{0}" cy="{1}" r="6" stroke="red" stroke-width="2" fill="none" stroke-dasharray='1,3' /> """.format(*self.__polarToCartesian__( abs(a), 50,(80,100)) )
#            self.body+= """<circle cx="{0}" cy="{1}" r="2" stroke="black" stroke-width="2" fill="black"  /> """.format(*self.__polarToCartesian__( abs(a), 51,(80,100)) ) 
#        for a in self.mids: self.body+= """<circle cx="{0}" cy="{1}" r="2" stroke="black" stroke-width="2" fill="black"  /> """.format(44.5,a)
#        rotation = 0
#        return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
#      <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate({5},100,100)"> {0}  </g> </svg>""".format(self.body,  self.x, self.y, *self._size,rotation)
    
#    def _repr_html_(self): return self.__display__()
    
#melon_graph(right_angles = [105,-80], left_angles=[255], mids=[100])
#melon_graph( left_angles=[270], mids=[100])
#melon_graph( left_angles=[285,-250], mids=[])


#class simple_ring_graph(ring_diagram):
#     def __init__(self,circ_dots = [], line_dots = [], is_contracted=True):
#         self.body=""
#         self._dashed_edges = [0]
#         self.circ_dots = circ_dots
#         self.line_dots = line_dots
#         self.is_contracted = is_contracted
#         self._is_directed=False
#         line_dots += [0, 60]
#         if not is_contracted: circ_dots += [90,270]
        
#     def __display__(self, ):      
#         rad = 30
#         center = (100,100)
#         dashed = "stroke-dasharray='1,5'" if self.is_contracted else ""
#         self.body+= """<circle cx="{0}" cy="{1}" r="{2}" stroke="black" stroke-width="2" fill="none" {3} /> """.format(100,100,rad,dashed)
#         self.body +=self.get_line((100,100-rad, 100,100+rad ), pid=0 if dashed else 1 )
        
#         for a in self.circ_dots:
            
#             pt= self.__polarToCartesian__(abs(a), rad, center)     
#             if a < 0: self.body+= """<circle cx="{0}" cy="{1}" r="4.5" stroke="red" stroke-width="1" fill="none" /> """.format(*pt)
#             self.body+= """<circle cx="{0}" cy="{1}" r="2.5" stroke="black" stroke-width="2" fill="black" /> """.format(*pt)
        
#         for a in self.line_dots:
#             if a < 0: self.body+= """<circle cx="{0}" cy="{1}" r="4.5" stroke="red" stroke-width="1" fill="none" /> """.format(100,100-rad+abs(a))
#             self.body+= """<circle cx="{0}" cy="{1}" r="2.5" stroke="black" stroke-width="2" fill="black" /> """.format(100,100-rad+abs(a))
          
#         if not self.is_contracted:
#             self.body += """<text x="65" y="75" style="fill: #000000;stroke: none; font-size: 12px;">k1</text>"""
#             self.body += """<text x="125" y="75" style="fill: #000000;stroke: none; font-size: 12px;">k2</text>"""  
#             self.body += """<text x="60" y="130" style="fill: #000000;stroke: none; font-size: 12px;">e3</text>"""
#             self.body += """<text x="125" y="130" style="fill: #000000;stroke: none; font-size: 12px;">e4</text>"""   
#             self.body += """<text x="102" y="110" style="fill: #000000;stroke: none; font-size: 12px;">e5</text>"""
        
        
#         if self.is_contracted:
#             self.body +=self.get_line((80,100-rad, 120,100-rad ) )
#         else:
#             self.body +=self.get_line((100-rad-20,100, 100-rad,100 ) )
#             self.body +=self.get_line((100+rad+20,100, 100+rad,100 ) )
            
            
#         rotation = 0
#         return  """<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" x="{1}" y="{2}" width="{3}" height="{4}">
#           <g fill="none" stroke="black" stroke-width="1.6" stroke-linecap="round" transform="rotate({5},100,100)"> {0}  </g> </svg>""".format(self.body,  10, 10, 200,200 ,rotation)
    
#     def _repr_html_(self): return self.__display__()
    
# simple_ring_graph(is_contracted=False, circ_dots=[320, 225])
# simple_ring_graph(is_contracted=False, circ_dots=[320, 135])
# simple_ring_graph(is_contracted=False, circ_dots=[320], line_dots=[30])
# simple_ring_graph(is_contracted=False, circ_dots=[45,135], line_dots=[])
# simple_ring_graph(is_contracted=False, circ_dots=[45,225], line_dots=[])
# simple_ring_graph(is_contracted=True, circ_dots=[270, -225])
# simple_ring_graph(is_contracted=True, circ_dots=[270, ], line_dots=[-30])
# simple_ring_graph(is_contracted=True, circ_dots=[90, -125])
# simple_ring_graph(is_contracted=True, circ_dots=[90, ], line_dots=[-30])
# simple_ring_graph(is_contracted=True, circ_dots=[270, ], line_dots=[-30])
# simple_ring_graph(is_contracted=True, circ_dots=[270, -90], line_dots=[])
# simple_ring_graph(is_contracted=True, circ_dots=[-90, ], line_dots=[30])
# simple_ring_graph(is_contracted=True, circ_dots=[ ], line_dots=[30, -45])
# simple_ring_graph(is_contracted=False, circ_dots=[135], line_dots=[30])
# simple_ring_graph(is_contracted=False, circ_dots=[225], line_dots=[30])
# simple_ring_graph(is_contracted=True, circ_dots=[270,90])
# simple_ring_graph(is_contracted=True, circ_dots=[120,60])
# simple_ring_graph(is_contracted=True, circ_dots=[360-60,360-120])
# simple_ring_graph(is_contracted=True, circ_dots=[270],line_dots=[30])