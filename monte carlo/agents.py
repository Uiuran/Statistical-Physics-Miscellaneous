#agent002, agent simulation in a General Topology with 1 float field for each agent (node) and random sampling of interactions (egdes) between time step, the process evaluation is a adition of the mean interaction to the fields { fieldA[time] =( fieldA[ time-1] + fieldB[time-1])/2 , the same is true for fieldB}, ie homophilic interaction, the agents start with random float number sampled from uniform [0,1). This process prints the visualization of the network.
#import paintGraph as pg
import numpy as np
import networkx as nx
import GraphMovie as gm
import igraph

class AgentModel(object):

    '''
    Agent models simulation
    '''

    def __init__(self, gml_file = "", model = { "model":"random", "node_number": 50, 'digraph': False, 'weighted': False, }, fields_dim = 1, time = 200):
        '''
        Constructor
        '''
        self.models = {
    
        "random":self.random,
        "scalefree":self.scalefree,
        "smallworld":self.smallworld,
    
        }
        self.time = time;
        a = self.topology( arquivo_gml = gml_file, mtype = model);
        self.fields = np.random.rand(fields_dim,len(a));      
        self.edges = [(v,j) for v in range(len(a)) for j in a[v]];
        for v in self.edges:                
            self.edges.pop(self.edges.index( ( v[1], v[0]) ) );                
        
    def random(self, param, digraph = False):
        ''' Returns a adjacency matrix for random graph in numpy array '''
        a = np.random.randint(0,2,(param["node_number"],param["node_number"]));
        if digraph == False:
            return (a + a.T)/2 + (a + a.T)%2 - np.diag(np.diag((a + a.T)/2 + (a + a.T)%2));
        else:
            return a;
   
    def scalefree(self, param):
        ''' Returns a adjacency list, not matrix, ie a list of lists with the number labels of the neighbors for each node. One must supply a dictionarie with the number of nodes of the final model and the number of edges to attach in each iteration of the preferential attachment rule: param["nodes_number"], param["edges_to_attach"] '''
        a = nx.generators.barabasi_albert_graph(param["nodes_number"], param["edges_to_attach"]);
        return np.array([v.keys() for v in a.adj.itervalues()]);
   
    def smallworld(self, param):
        ''' Returns a adjacency list, not matrix, ie a list of lists with the number labels of the neighbors for each node. One must supply a dictionarie with the number of nodes of the final model, the number of nearest neighbors to start attached in the ring topology and the probability of re-attachment: param["nodes_number"], param["k"], param["p"] '''
        a = nx.generators.watts_strogatz_graph(param["nodes_number"], param["k"], param["p"]);
        return np.array([v.keys() for v in a.adj.itervalues()]);
 
 
    def runSimulation(self, plot = True):
        full = self.time;
        if plot == True:
            self.movie = gm.GraphMovie();
#        attrs = {'fillcolor':fields ,'width':0.5,'height':0.5,'fontsize':0.0,'nodedegreesize':False, 'plotlabel':'agent0_'+str(full-ttime) }
#     A = nx.to_agraph(net);
#   1 A.layout(prog = 'twopi');
        while full > 0:
            if plot == True and self.time == full:  
                 self.net = igraph.read(self.string+".gml");          
                 for i,v in enumerate(self.net.vs):
                     v['label'] = str(i);
                     v['color'] = self.fields[0,i];
                 self.movie.addGraph(self.net);
            elif plot == True and self.time != full:
                 self.net = igraph.read(self.string+".gml")
                 for i,v in enumerate(self.net.vs):
                     v['label'] = str(i);
                     if i == self.edges[sampvar[0]][0] or i == self.edges[sampvar[0]][1]:
                         v['color'] = 0.0;
                         
                     else:
                         v['color'] = self.fields[0,i];                  
#                 e = self.net.es[int(sampvar[0])];
#                 print e;
#                 print "edges"
#                 print sampvar[0]
#                 print self.edges[sampvar[0]][0]; print self.edges[sampvar[0]][1];
#                 e['color'] = 0.0;                                   # set edge color 
                 self.movie.addGraph(self.net); 
            sampvar = np.random.randint(0,len(self.edges),(1)); # the edge(interaction) sampling process
            interact = float(self.fields[0,self.edges[sampvar[0]][0]]+self.fields[0,self.edges[sampvar[0]][1]])/2;
            self.fields[0,self.edges[sampvar[0]][0]] = float(self.fields[0,self.edges[sampvar[0]][0]] + interact)/2;
            self.fields[0,self.edges[sampvar[0]][1]] = float(self.fields[0,self.edges[sampvar[0]][1]] + interact)/2; # attrs['fillcolor'] = fields; 
            full = full - 1;
#            attrs['plotlabel'] = 'agent0_'+str(full-ttime);      
        self.movie.doMovieLayout();
        self.movie.renderMovie(name = self.string);
   
    def topology(self, arquivo_gml = "", mtype = { "model":"random", "node_number":50 , 'digraph' : False, 'weighted' : False, } ):
        if arquivo_gml == "":
            param = mtype.copy();
            param.pop("model");
            return self.models[mtype["model"]](param);
        else:
            self.net = igraph.read(arquivo_gml+".gml");
            self.string = arquivo_gml;                          
            return np.array([self.net.neighbors(v) for v in range(len(self.net.vs))]);
    
    
    
#    ttime = 5; # total time steps
#    a = topology("/home/humannoise/Desktop/UOne/academic/complex_netw/data/facebook/humannoise03_03_2013.gml");  #adjacency matrix for general topology .gml file
    # a = topology(mtype = {"model":"scalefree","nodes_number":600,"edges_to_attach":3,});
    # a = topology(mtype = {"model":"smallworld","nodes_number":600,"k":4,"p":0.6});
#    fields = np.random.rand(len(a)); # node fields, or agents atributs

    ### The part below account the edges if you use a random graph with adjacency matrix, if you use networkx and adjacency list dont need to do this
    # nonz = randomnet.nonzero();  2 dim-tuple with array counting indexes of connected nodes
    # nonznum = len(nonz[0][:]);  number of edges
    # edges = [(nonz[0][v],nonz[1][v]) for v in range(nonznum)];  the edges of the random graph
#    edges = [(v,j) for v in range(len(a)) for j in a[v]];
    
#    runSimulation(fields, edges, ttime, plot = True);
