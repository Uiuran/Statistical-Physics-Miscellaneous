'''
@Forked by Uiuran (c) 2013
@date Feb 21, 2011
@original author: mjbommar (c) 2011
@license Simplified BSD, (C) 2011.
 
Free to use for non-commerical purposes. Attribution appreciated :)
'''
 
import igraph
import os, os.path
import numpy 
import pylab as p
import time
import cairo
import operator
 
class GraphMovie(object):

    '''
    Graph Movie Object.
    '''
 
    def __init__(self, graphSequence = [], labelSequence = None, dimension = 3, **kwargs):
        '''
        Constructor

        '''
 
        '''
        This is the number of interpolating frames per sequence element.
        Larger values mean smoother transitions but larger videos.
        '''
        self.interpFrames = 5
        '''
        Number of delay frames for start and end of movie.
        '''
        self.delayFrames = 5
        '''
        This is the number of frames per second in the resulting video.
        '''
        self.fps = 10 
 
        '''
        Boolean for node labelling.
        '''
        self.labelNodes = False
        '''
        Invert colors?
        '''
        self.invertColors = False
        '''
        Number of iterations for KK.
        '''
        self.kkIterations = 1000
        '''
        Frame counter.
        '''
        self.frame = 0
        '''
        Movie dimensions and margins.
        '''
        self.dim = dimension
        self.pixelsX = 800
        self.pixelsY = 600
        self.margin = self.pixelsY / 10.0
           
        '''
        Output path.
        Trailing slash!
        '''
        self.frameDirectory = 'frames/'
        if not os.path.exists(self.frameDirectory):
            os.mkdir(self.frameDirectory)
        '''
        This is the list of graph objects.
        '''
        self.graphSequence = graphSequence
        if labelSequence:
            self.labelSequence = labelSequence
        else:
            self.labelSequence = map(str, range(len(graphSequence)))
        if len(self.graphSequence) > 0:
            self.checkGraphSequence()



    def project2D(self,layout, alpha, beta):
        '''
        This method will project a set of points in 3D to 2D based on the given
        angles alpha and beta.
        '''
# Calculate the rotation matrices based on the given angles.
        c = numpy.matrix([[1, 0, 0], [0, numpy.cos(alpha), numpy.sin(alpha)], [0, -numpy.sin(alpha), numpy.cos(alpha)]])
        c = c * numpy.matrix([[numpy.cos(beta), 0, -numpy.sin(beta)], [0, 1, 0], [numpy.sin(beta), 0, numpy.cos(beta)]])
        b = numpy.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
# Hit the layout, rotate, and kill a dimension
        layout = numpy.matrix(layout)
        X = (b * (c * layout.transpose())).transpose()
        return [[X[i,0],X[i,1],X[i,2]] for i in range(X.shape[0])]
 
    def drawGraph3D(self,graph, layout, angle, fileName):
        '''
        Draw a graph in 3D with the given layout, angle, and filename.
        '''
# Setup some vertex attributes and calculate the projection
        graph.vs['degree'] = graph.degree()
        vertexRadius = 0.2 * (0.9 * 0.9) / numpy.sqrt(graph.vcount())
        graph.vs['x3'], graph.vs['y3'], graph.vs['z3'] = zip(*layout)
        layout2D = self.project2D(layout, angle[0], angle[1])
        graph.vs['x2'], graph.vs['y2'], graph.vs['z2'] = zip(*layout2D)
        minX, maxX = min(graph.vs['x2']), max(graph.vs['x2'])
        minY, maxY = min(graph.vs['y2']), max(graph.vs['y2'])
        minZ, maxZ = min(graph.vs['z2']), max(graph.vs['z2'])
# Calculate the draw order. This is important if we want this to look
# realistically 3D.
        zVal, zOrder = zip(*sorted(zip(graph.vs['z3'], range(graph.vcount()))))
# Setup the cairo surface
        surf = cairo.ImageSurface(cairo.FORMAT_ARGB32, 1280, 800)
        con = cairo.Context(surf)
        con.scale(1280.0, 800.0)
# Draw the background
        con.set_source_rgba(255.0, 255.0, 255.0, 1.0)
        con.rectangle(0.0, 0.0, 1.0, 1.0)
        con.fill()
# Draw the edges without respect to z-order but set their alpha along
# a linear gradient to represent depth.
        for e in graph.es:
# Get the first vertex info
            v0 = graph.vs[e.source]
            x0 = (v0['x2'] - minX) / (maxX - minX)
            y0 = (v0['y2'] - minY) / (maxY - minY)
            alpha0 = (v0['z2'] - minZ) / (maxZ - minZ)
            alpha0 = max(0.1, alpha0)
# Get the second vertex info
            v1 = graph.vs[e.target]
            x1 = (v1['x2'] - minX) / (maxX - minX)
            y1 = (v1['y2'] - minY) / (maxY - minY)
            alpha1 = (v1['z2'] - minZ) / (maxZ - minZ)	
            alpha1 = max(0.1, alpha1)
# Setup the pattern info
            pat = cairo.LinearGradient(x0, y0, x1, y1)
            pat.add_color_stop_rgba(1.0, 0.0, 0.0, 0.0, alpha0/ 1.0)  # 1.0,0,0,0 for background with white rgba, 1.0,1.0,1.0,1.0 for black background
#            try:
#                if e.attributes().has_key('color'):
#                    try:
#                        pat.add_color_stop_rgba(1, e['color'], e['color'], e['color'], alpha1/ 1.0)    # one must use smaller alpha with dark background a bigger with white   
#                    except TypeError:
#                        print "Type Error"
#            except NameError:                
            pat.add_color_stop_rgba(1,1.0,1.0,1.0,alpha1/1.0);
            con.set_source(pat)
# Draw the line
            con.set_line_width(vertexRadius / 4.0)
            con.move_to(x0, y0)	
            con.line_to(x1, y1)
            con.stroke()
# Draw vertices in z-order
        for i in zOrder:
            v = graph.vs[i]
#            print v['color']
#            print i
            alpha = (v['z2'] - minZ) / (maxZ - minZ)
            alpha = max(0.1, alpha)
            if v['color'] != 0.0:
                radius = 0.5 * ((-numpy.log(numpy.sqrt(2*numpy.pi)*v['color']))* 0.5) / numpy.sqrt(graph.vcount()) # use degree(use log degree) or fields(use a larg multiplier) to vertex radius
            else:
                radius = 0.02 
            x = (v['x2'] - minX) / (maxX - minX)
            y = (v['y2'] - minY) / (maxY - minY)
# Setup the radial pattern for 3D lighting effect
            pat = cairo.RadialGradient(x, y, radius / 4.0, x, y, radius)
            try:
                arrombado = [v['color'],0,0]
#                arrombado = [alpha*v['color'],alpha*v['color'],alpha*v['color']]
#                arrombado = map(operator.div,map(operator.add,p.cm.YlGn(int(v['color']*255))[:-1],p.cm.autumn(int(v['color']*255))[:-1]),(2,2,2))
            except KeyError:
                arrombado = map(operator.div,map(operator.add,p.cm.YlGn(int(alpha*255))[:-1],p.cm.autumn(int(alpha*255))[:-1]),(2,2,2))
            if 'arrombado' not in locals():     
                arrombado = map(operator.div,map(operator.add,p.cm.YlGn(int(alpha*255))[:-1],p.cm.autumn(int(alpha*255))[:-1]),(2,2,2))
            pat.add_color_stop_rgba(0,arrombado[0]*alpha,arrombado[1]*alpha,arrombado[2]*alpha,1)
            del arrombado
            pat.add_color_stop_rgba(1,0,0,0,1)
            con.set_source(pat)
# Draw the vertex sphere
            con.move_to(x, y)
            con.arc(x, y, radius, 0, 2 * numpy.pi)	
            con.fill()
# Output the surface
        surf.write_to_png(fileName)

    def checkGraphSequence(self):
        '''
        Check to make sure that our current graph sequence is OK.
        * Every element must be an igraph object with at least 1 vertex.

        '''
        for i, graph in enumerate(self.graphSequence):
            if not isinstance(graph, igraph.Graph):
                raise Exception("GraphMovie::checkGraphSequence: element {0} is not an igraph.Graph object.".format(i))
            if graph.vcount() == 0:
                raise Exception("GraphMovie::checkGraphSequence: element {0} has 0 vertices.".format(i))
 
    def addGraph(self, graph, frameLabel = None):
        '''
        Add a single graph element to the sequence.
        '''
        if not isinstance(graph, igraph.Graph):
            raise Exception("GraphMovie::addGraph: graph is not an igraph.Graph object.")
        if graph.vcount() == 0:
            raise Exception("GraphMovie::addGraph: graph has 0 vertices.")
        self.graphSequence.append(graph)
        if frameLabel:
            self.labelSequence.append(frameLabel)
        else:
            self.labelSequence.append(str(len(self.graphSequence)))
 
    def interpolateLayout(self, lastLabels, lastLayout, g):
        '''
        Calculate the interpolated layout.
        '''
        labels = [v['label'] for v in g.vs]
        seedLayout = []
        meanX = sum([p[0] for p in lastLayout]) / float(len(lastLayout))
        meanY = sum([p[1] for p in lastLayout]) / float(len(lastLayout))
        if self.dim == 3:
            meanZ = sum([p[2] for p in lastLayout]) / float(len(lastLayout))  
        '''
        This is old code and should be optimized, but oh well.
        Maybe next version...
        '''
        for label in labels:
            if label in lastLabels:
                index = lastLabels.index(label)
                seedLayout.append(lastLayout[index])
            else:
                if self.dim == 3:
                    seedLayout.append([meanX,meanY,meanZ])
                else:
                    seedLayout.append([meanX,meanY])    
        if self.dim == 3:
            layout = g.layout_kamada_kawai_3d(seed = seedLayout, maxiter = self.kkIterations) 
            layoutDiff =[[layout[i][0] - seedLayout[i][0], layout[i][1] - seedLayout[i][1],layout[i][2]-seedLayout[i][2]] for i in range(len(seedLayout))] 
        else:
            layout = g.layout_kamada_kawai(seed = seedLayout, maxiter = self.kkIterations)
            layoutDiff =[[layout[i][0] - seedLayout[i][0], layout[i][1] - seedLayout[i][1]] for i in range(len(seedLayout))]
        interpLayout = []
 
        '''
        This is just simple linear interpolation.
        I would like to try something else, but for now this works!
        '''
        for i in range(self.interpFrames):
            c = float(i) / (self.interpFrames - 1)
            if self.dim == 3:
	        interpLayout.append([[seedLayout[i][0] + c*layoutDiff[i][0],seedLayout[i][1] + c*layoutDiff[i][1], seedLayout[i][2] + c*layoutDiff[i][2]] for i in range(len(seedLayout))]) 
            else:
                interpLayout.append([[seedLayout[i][0] + c*layoutDiff[i][0],seedLayout[i][1] + c*layoutDiff[i][1]] for i in range(len(seedLayout))])
        return interpLayout

    def plotFrame(self, g, label, layout):
        '''
        Plot a frame of the movie.
        '''
        if self.dim == 3:
            frameFile = self.frameDirectory + '%08d.png' % (self.frame)
            if not self.labelNodes:
                for v in g.vs:
                    v['label_size'] = 2. 
            self.frame += 1
            self.alpha = float(self.frame*numpy.pi)/(2*self.delayFrames+self.interpFrames*float(len(self.labelSequence)))
            self.beta = float(self.frame*numpy.pi)/(2*self.delayFrames+self.interpFrames*float(len(self.labelSequence)))     
            self.drawGraph3D(g, layout, (self.alpha,self.beta), frameFile)
           
        else:
            frameFile = self.frameDirectory + '%08d.png' % (self.frame)
            if not self.labelNodes:
                for v in g.vs:
                    v['label_size'] = 0.001 
            p = igraph.drawing.plot(g, layout = layout, bbox = (self.pixelsX, self.pixelsY), margin = self.margin, target = frameFile)
            '''
            Now label and possibly invert the frames.
            '''
            self.frame += 1

    def doMovieLayout(self):
        '''
        Calculate the graph movie layout.
        '''
        lastLayout = None
        lastLabels = None
        for i, graph in enumerate(self.graphSequence):
            label = self.labelSequence[i]
            if i == 0:
                '''
                This is the first element in the graph sequence, so we have no meaningful initial conditions.
                '''
                if self.dim == 3:
                    layout = graph.layout_grid_3d()
                    layout = graph.layout_kamada_kawai_3d(maxiter = self.kkIterations * 10, seed = layout) 
                else:
                    layout = graph.layout_grid()
		    layout = graph.layout_kamada_kawai(maxiter = self.kkIterations * 10, seed = layout)
                for j in range(self.delayFrames):
                    self.plotFrame(graph, label, layout)
                lastLayout = layout
                lastLabels = [v['label'] for v in graph.vs]
            else:
                '''
                We need to use the previous graph element initial conditions so the movie makes sense.
                '''
                interpLayout = self.interpolateLayout(lastLabels, lastLayout, graph)
                for layout in interpLayout:
                    self.plotFrame(graph, label, layout)
                lastLayout = interpLayout[-1]
                lastLabels = [v['label'] for v in graph.vs]
        for j in range(self.delayFrames):
            self.plotFrame(graph, label, layout)

    def renderMovie(self, name='output'):
        '''  
        Render the movie with mencoder.
        '''
        cmdString = 'cd frames/ && mencoder -noskip mf://*.png -mf fps='+str(self.fps)+':type=png -ovc lavc -lavcopts vcodec=mpeg4 -o '+name+'.avi'.format(self.frameDirectory, self.fps)
        os.system(cmdString);
