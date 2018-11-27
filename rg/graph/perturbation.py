

#maybe find somew
#combinations over the primitives choose 4           
def _shuffle_(A,B):
    for s in incidence_matrix.shuffle(A,B):
        if s != (0,0): #the null case     
             yield A.dprod(B,species_chord_set=s)
                
def brute_mul(A,B,C,D,l=1):
    genset = list(_shuffle_(A,B))
    news = []
    for g in genset:news += list(_shuffle_(g,C))
    genset+= news
    for g in genset:news += list(_shuffle_(g,D))
    genset+= news
    return  [g for g in genset if g.first_betti_number == l and g.bridge_count == 0]# and


