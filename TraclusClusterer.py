# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 08:45:07 2017

@author: Eyal
"""

# for geometry and plotting
from shapely.geometry import MultiPoint, GeometryCollection, LineString
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from ClusterProcessor import ClusterProcessor
from Projector import EquirectangularProjector

# for clustering
from IMCObjects import Segment, Trajectory #GLine
from networkx import DiGraph, shortest_path, descendants
from networkx import exception as nxe
import itertools

class SegmentsClusterer:
    
    def __init__(self, segments_list_of_trajectory_collection, eps, MinLns):
        self.segmentsList = segments_list_of_trajectory_collection
        self.eps = eps
        self.MinLns = MinLns
        self.segmentsMatrix = None
        self.directReachablityGraph = None
#        print("SegmentsClusterer Init: building graph...")
#        self.directReachablityGraph = self.computeDirectReachabilityGraph()
    
    def initActions(self):
        print("InitActions of SegmentsClusterer:")
        print("New instance of SegmentsMatrix...")
        self.segmentsMatrix = SegmentsMatrix(self.segmentsList, self.eps)
        print("Starting Matrix' projector...")
        self.segmentsMatrix.startMyProjector()
        print("Constructing Matrix...")
        self.segmentsMatrix.constructMatrix()
        print("Building the graph...")
        self.directReachablityGraph = self.computeDirectReachabilityGraph()
        print("SegmentsClusterer initialization is complete!")
        
    def computeDirectReachabilityGraph(self):
        '''
        an edge in the graph v1->v2 tells us that v2 is dirReachable from v1
        an edge in this graph connects from a core lineSeg to a dirReach node
        '''
        dig = DiGraph() #Directed Graph
        dig.add_nodes_from(self.segmentsList)
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
                print("LineSegmentClustering: added cluster id {}".format(cid))
                # remove nodes from possibility:
                newPossibleNodesList = [nd for nd in possibleNodesList if not nd in clusters[cid]]
                possibleNodesList = newPossibleNodesList
                cid += 1
        # todo: remove clusters with no variaty of trajectories
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
        proj = self.segmentsMatrix.projector
        colors = iter(cm.rainbow(np.linspace(0, 1, len(clusterList))))
        clusterProcessor = ClusterProcessor(proj) 
#        i = 0
        for cluster in clusterList:
            ccolor = next(colors)
            for seg in cluster:
                xs = [seg[0].shapelyPoint.x, seg[1].shapelyPoint.x]
                ys = [seg[0].shapelyPoint.y, seg[1].shapelyPoint.y]
                axs.plot(xs, ys, color = ccolor, alpha=0.5, linestyle = '--')
            
            # make convex hull:
            hull = MultiPoint([seg[0].shapelyPoint for seg in cluster] + [seg[1].shapelyPoint for seg in cluster]).convex_hull
            if type(hull) == GeometryCollection or type(hull) == LineString:
                continue
            else:
                xs, ys = hull.exterior.xy
                axs.plot(xs, ys, color = ccolor, linewidth = 2, alpha=0.5, linestyle = '--')
            
            clusterProcessor.loadGeoCluster_segmentsList(cluster)
            start, end = clusterProcessor.calcAverageDirection()
            lat1, long1 = proj.XYToLatLong(start[0], start[1])
            lat2, long2 = proj.XYToLatLong(end[0], end[1])
            axs.plot([long1, long2], [lat1, lat2], color = ccolor, linewidth = 5)
    
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
    def __init__(self, segments_list_of_trajectory_collection, eps):
        self.segmentList = segments_list_of_trajectory_collection
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
        self.projector = None
        self.eps = eps
        
    def startMyProjector(self):
        self.projector = EquirectangularProjector(self.segmentList)
    
    
    def segToMatrixInx(self, seg):
        '''
        returns the index of the matrix the segment belongs to
        '''
        seg_end1_xy = self.projector.latLongToXY(seg[0][0], seg[0][1])
        seg_end2_xy = self.projector.latLongToXY(seg[1][0], seg[1][1])
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
        i = int((mx - self.projector.originXY[0]) / self.eps)
        j = int((my - self.projector.originXY[1]) / self.eps)
        return (i,j)
        
    def constructMatrix(self):
        myMap = {}
        for seg in self.segmentList:
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
    

