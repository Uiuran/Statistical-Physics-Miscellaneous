'''
@Fork by Uiuran 2013
@created on Jan 17, 2010
@author: mjbommar
'''
 
from GraphMovie import GraphMovie
import igraph
import numpy as np
 
if __name__ == "__main__":
    '''
    Create the GraphMovie object.
    '''
    m = GraphMovie(dimension = 3)
    '''
    Create the first graph.
    '''
    edges = [(0,1),(1,2),(1,3)]
    g = igraph.Graph(edges)
    for i,v in enumerate(g.vs):
        v['label'] = str(i)
        v['color'] = 1.0
    m.addGraph(g)
 
    '''
    Now add an edge.
    '''
    edges.append((4,1))
    g = igraph.Graph(edges)
#    g.add_vertex(4)
#    g.add_edge(4,1);
    for i,v in enumerate(g.vs):
        v['label'] = str(i)
    m.addGraph(g)
    for j in range(0,2):
        g = igraph.Graph(edges)
        for i,v in enumerate(g.vs):
            v['label'] = str(i)
            if np.mod(i,2) == 0:
                v['color'] = 0.0;
            else:
                v['color'] = 1.0; 
        m.addGraph(g)

    '''
    Now add a few edges!
    '''
    edges.extend([(5,2),(6,5),(7,5),(8,5),(8,7),(9,8)])
    g = igraph.Graph(edges)
    for i,v in enumerate(g.vs):
        v['label'] = str(i) 
    m.addGraph(g)
 
    '''
    Now process the layouts, render the frames, and generate the movie.
    '''
    m.doMovieLayout()
    m.renderMovie('output1')
