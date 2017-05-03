# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 08:45:07 2017

@author: Eyal
"""

# for geometry and plotting
from matplotlib import pyplot as plt
import numpy as np
import statistics
from shapely.geometry import MultiPoint, GeometryCollection, LineString
from shapely.geometry import Point as shapelyPoint
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from ClusterProcessor import ClusterProcessor
from Projector import EquirectangularProjector

# for clustering
from IMCObjects import Segment, Trajectory, LSegment #GLine
from networkx import DiGraph, shortest_path, descendants
from networkx import exception as nxe
import itertools

class SegmentsClusterer:
    
    def __init__(self, segments_list_of_trajectory_collection, eps, MinLns):
        self.segmentsList = segments_list_of_trajectory_collection
        self.lsegmentsList = None
        self.eps = eps
        self.MinLns = MinLns
        self.segmentsMatrix = None
        self.directReachablityGraph = None
        self.projector = None
#        print("SegmentsClusterer Init: building graph...")
#        self.directReachablityGraph = self.computeDirectReachabilityGraph()
    
    def initActions(self):
        print("InitActions of SegmentsClusterer:")
        
        print("Starting Matrix' projector...")
        self.startMyProjector()
        
        print("Convert all segments to local xy")
        self.convertSegmentsToLocalXY()
        
        #### TODO: debug   #########
        if False:
            fig, axs = plt.subplots()
            self.plotAllLocalSegments(axs)
        
        print("New instance of SegmentsMatrix...")
        self.segmentsMatrix = SegmentsMatrix(self.lsegmentsList, self.eps)
        
        print("Constructing Matrix...")
        self.segmentsMatrix.constructMatrix()
        
        print("Matrix constructed!")
#        for inx, val in self.segmentsMatrix.spInxMap.items():
#            print("{}: {} items".format(inx, len(val)))
        
        print("Building the graph...")
        self.directReachablityGraph = self.computeDirectReachabilityGraph()
        
        print("SegmentsClusterer initialization is complete!")
    
    def plotAllLocalSegments(self, axs):
        for gseg in self.lsegmentsList:
            axs.plot([gseg.end1[0], gseg.end2[0]], [gseg.end1[1], gseg.end2[1]])
            
    def startMyProjector(self):
        self.projector = EquirectangularProjector(self.segmentsList)
        
    def convertSegmentsToLocalXY(self):
        self.lsegmentsList = [LSegment(self.projector.latLongToXY(gseg.end1[0], gseg.end1[1]),
                                       self.projector.latLongToXY(gseg.end2[0], gseg.end2[1]),
                                       gseg.status, gseg.trajIndex) for gseg in self.segmentsList]
        
        
    def computeDirectReachabilityGraph(self):
        '''
        an edge in the graph v1->v2 tells us that v2 is dirReachable from v1
        an edge in this graph connects from a core lineSeg to a dirReach node
        '''
        dig = DiGraph() #Directed Graph
        dig.add_nodes_from(self.lsegmentsList)
        i = 0
        for anode in dig.nodes(): #TODO: might be sped up with thread workers?
#            print("computeDirectReachabilityGraph: calculating edges of node #{}".format(i))
            i += 1
            a_epsN = self.eps_neighborhood_of_seg(anode)
            if len(a_epsN) < self.MinLns: # anode is not a core segment
                anode.status = 2 # noise
                continue
            
            anode.status = 1 # signal(core)
            for bnode in a_epsN:
                if anode is bnode:
                    continue
                dig.add_edge(anode , bnode)
        return dig
    
    
    def LineSegmentClustering(self):
        '''
        Algorithm Line Segment Clustering
        works after the computeDirectReachabilityGraph was called
        returns a dictionary of <cid, cluster>
        '''
        clusters = {}
        cid = 0
        possibleNodesList = list(self.directReachablityGraph.nodes())
        
        while len(possibleNodesList) > 0:
            node = possibleNodesList[0]
            possibleNodesList = possibleNodesList[1:]
            if node.status == 1:
                clusters[cid] = list(descendants(self.directReachablityGraph, node)) + [node]
#                print("LineSegmentClustering: added cluster id {}".format(cid))
                # remove nodes from possibility:
                newPossibleNodesList = [nd for nd in possibleNodesList if not nd in clusters[cid]]
                possibleNodesList = newPossibleNodesList
                cid += 1
#       remove clusters with no variaty of trajectories
#        for cIndex, cluster in list(clusters.items()):
#            trajIndexesSet = set([seg.trajIndex for seg in cluster]) # removes duplicates
#            if len(trajIndexesSet) < self.MinLns:
#                del clusters[cIndex]
        return clusters

        
    def eps_neighborhood_of_seg(self, Li):
        '''
        Traclus article Definition 4
        '''
        matrixIndex = self.segmentsMatrix.segToMatrixInx(Li)
        segmentsFrom3x3List = self.segmentsMatrix.getSegmentsFrom3x3(matrixIndex)
        return [Lj for Lj in segmentsFrom3x3List if Li.myDistance(Lj) < self.eps] # used Gline's myDistance


    def plotClusters(self, axs, clusterList):   
        '''
        shows in local xy coords
        '''
        axs.set_title("Clustering eps = {} MinLns = {}".format(self.eps, self.MinLns))
        colors = iter(cm.rainbow(np.linspace(0, 1, len(clusterList))))
        clusterProcessor = ClusterProcessor(self.MinLns) 
        print("plotClusters: {} clusters to plot.".format(len(clusterList)))
        for i, cluster in enumerate(clusterList):
#            if i != 2:
#                continue
#            print("plotting cluster i = {}".format(i))
            ccolor = next(colors)
            clusterProcessor.loadLocalCluster_lsegList(cluster)
            
            #################################
            #       Representative          #
            #################################
            
            repr_And_walls = clusterProcessor.calcRepresentativeLocal(30)
            rxs = [raw[0][0] for raw in repr_And_walls]
            rys = [raw[0][1] for raw in repr_And_walls]
            if len(rxs) > 1:
                axs.plot(rxs, rys, color = ccolor, linewidth = 5)
#                cols = np.arange(len(rxs))
#                axs.scatter(rxs, rys, c = cols, s =150)
            else:
                print("cluster {} had no good representative".format(i))
                continue
            
            rxs = [raw[1][0] for raw in repr_And_walls] + [raw[2][0] for raw in repr_And_walls][::-1]
            rxs += [rxs[0]]
            rys = [raw[1][1] for raw in repr_And_walls] + [raw[2][1] for raw in repr_And_walls][::-1]
            rys += [rys[0]]
            axs.plot(rxs, rys, color = ccolor, linewidth = 5)
#            axs.scatter(rxs, rys, color = 'b', s =150)
#            
#            rxs = [raw[2][0] for raw in repr_And_walls]
#            rys = [raw[2][1] for raw in repr_And_walls]
#            axs.plot(rxs, rys, color = ccolor, linewidth = 5)
#            axs.scatter(rxs, rys, color = 'g', s =150)
            
            
            for seg in cluster:
                e1x, e1y = seg.end1
                e2x, e2y = seg.end2
                xs = [e1x, e2x]
                ys = [e1y, e2y]
                axs.plot(xs, ys, color = ccolor, alpha=0.5, linestyle = '--')
#            
#            # make convex hull:
#            hull = MultiPoint([shapelyPoint(seg.end1[0], seg.end1[1]) for seg in cluster] +\
#                                [shapelyPoint(seg.end2[0], seg.end2[1]) for seg in cluster]).convex_hull
#
#            if type(hull) == GeometryCollection or type(hull) == LineString:
#                continue
#            else:
##                print("made hull")
#                xs, ys = hull.exterior.xy
#                axs.plot(xs, ys, color = ccolor, linewidth = 2, alpha=0.5) #, linestyle = '--')
            

            #################################
            #   AvgDir Plotting             #
            #################################
            
#            avgDirLLine = clusterProcessor.calcAverageDirection()
#            axs.plot([avgDirLLine.end1[0], avgDirLLine.end2[0]], [avgDirLLine.end1[1], avgDirLLine.end2[1]], color = ccolor, linewidth = 5)
#            axs.scatter([avgDirLLine.end1[0]],[avgDirLLine.end1[1]], color = 'r', s =150)
            
            
            
#            print("cluster i={}, reprAndWals = {}".format(i, repr_And_walls))
            
            
            
                
                
    
    # UNUSED METHODS #
    def isCoreLineSegment(self, Li):
        '''
        Traclus article Definition 5
        '''
        raise Exception("Don't use this highly inefficient function")
        return len(self.eps_neighborhood_of_seg(Li, self.eps)) >= self.MinLns
    
    def isDirectlyDensityReachable(self, line, fromLine):
        '''
        Traclus article Definition 6
        '''
        raise Exception("Don't use this highly inefficient function")
        if line in self.eps_neighborhood_of_seg(fromLine, self.eps):
            if len(self.eps_neighborhood_of_seg(fromLine, self.eps)) >= self.MinLns:
                return True
        return False
     
    def isDensityReachable(self, line, fromLine):
        '''
        Traclus article Definition 7
        TODO: currently inefficient (read about R-Tree for optimization)
        '''
        try:
            sp = shortest_path(self.directReachablityGraph, source = fromLine, target = line)
            if len(sp)>0:
                return True
            else:
                return False
        except nxe.NetworkXNoPath as e:
            return False
        
    def isDensityConnected(self, line, toLine):
        '''
        Traclus article Definition 8
        '''
        if any(self.isDirectlyDensityReachable(line, Lk) 
                and self.isDirectlyDensityReachable(toLine, Lk) 
                    for Lk in self.directReachablityGraph.nodes()):
            return True
        return False
    
    def isDensityConnectedSet(self, clusterListOfSegments):
        '''
        Traclus article Definition 9
        '''
        if len(clusterListOfSegments) == 0:
            return False
        
        connectivity = all(self.isDensityConneced(prod[0], prod[1]) for prod in itertools.product(clusterListOfSegments, clusterListOfSegments))
        if not connectivity:
            return False
        
        Maximality =  all(prod[1] in clusterListOfSegments for
            prod in itertools.product(self.directReachablityGraph.nodes(), self.directReachablityGraph.nodes())
            if prod[0] in clusterListOfSegments and self.isDensityReachable(prod[1], prod[0]))
        return Maximality
       
    
class SegmentsMatrix:
    def __init__(self, lsegments_list, eps):
        self.lsegmentList = lsegments_list
        '''
        <spatialIndex, Segment Class Object List> dictionary
        where spatialIndex is a tuple of (i, j) - 
        the index in the matrix of the slot which starts at
        (originX + i * eps, originY + j * eps)
        
            j
            
            |
            |
            |
            |___________ i
        originXY    
        '''
        self.spInxMap = {}
        self.eps = eps
    
    
    def segToMatrixInx(self, lseg):
        '''
        returns the index of the matrix the segment belongs to
        '''
        seg_end1_xy = lseg.end1 #self.projector.latLongToXY(seg[0][0], seg[0][1])
        seg_end2_xy = lseg.end2 #self.projector.latLongToXY(seg[1][0], seg[1][1])
        xm = (seg_end1_xy[0] + seg_end2_xy[0]) / 2.0
        ym = (seg_end1_xy[1] + seg_end2_xy[1]) / 2.0
        segMXY = (xm, ym)
        segMatrixIndex = self.segMXYToMatrixInx(segMXY)
        return segMatrixIndex
    
    def segMXYToMatrixInx(self, segMXY):
        '''
        segMXY inx must be greater than originXY in both dimensions
        '''
        mx, my = segMXY
        i = int(mx / self.eps)
        j = int(my / self.eps)
#        i = int((mx - self.projector.originXY[0]) / self.eps)
#        j = int((my - self.projector.originXY[1]) / self.eps)
        return (i,j)
        
    def constructMatrix(self):
        myMap = {}
        for seg in self.lsegmentList:
            segMatrixIndex = self.segToMatrixInx(seg)
            if segMatrixIndex in myMap:
                myMap[segMatrixIndex].append(seg)
            else:
                myMap[segMatrixIndex] = [seg]
                
        self.spInxMap = myMap
        
    def getSegmentsFrom3x3(self, spatialIndex):
        i,j = spatialIndex
        if i < 0 or j < 0:
            raise ValueError("Indices are not negative!")
        
        result = []
        if spatialIndex in self.spInxMap:
            for ii in range(i-1, i+2):
                for jj in range(j-1, j+2):
                    if ii < 0 or jj < 0:
                        continue
                    if (ii,jj) in self.spInxMap:
                        result += self.spInxMap[(ii,jj)]
        return result
    

