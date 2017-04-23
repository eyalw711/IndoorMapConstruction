# -*- coding: utf-8 -*-
"""
Created on Sun Apr 23 08:45:07 2017

@author: Eyal
"""

from objects import GLine
from networkx import DiGraph, shortest_path
from networkx import exception as nxe
import itertools

class SegmentsClusterer:
    
    def __init__(self, GLineList, eps, MinLns):
        self.segmentsList = [segment(gline[0], gline[1]) for gline in GLineList]
        self.eps = eps
        self.MinLns = MinLns
        self.directReachablityGraph = self.computeDirectReachabilityGraph()
    
    def computeDirectReachabilityGraph(self):
        '''
        an edge in the graph v1->v2 tells us that v2 is dirReachable from v1
        '''
        dig = DiGraph() #Directed Graph
        dig.add_nodes_from(self.segmentsList)
        for prod in itertools.product(dig.nodes(), dig.nodes()):
            if prod[0] is prod[1]:
                continue
            if self.isDirectlyDensityReachable(prod[0], prod[1]):
                dig.add_edge(prod[1], prod[2])
        return dig
        
    def eps_neighborhood_of_seg(self, Li):
        '''
        Traclus article Definition 4
        '''
        return [Lj for Lj in self.segmentsList if segment.distance(Li,Lj) < self.eps]
        
    def isCoreLineSegments(self, Li):
        '''
        Traclus article Definition 5
        '''
        return len(self.eps_neighborhood_of_seg(Li, self.eps)) > self.MinLns
    
    def isDirectlyDensityReachable(self, line, fromLine):
        '''
        Traclus article Definition 6
        '''
        if line in self.eps_neighborhood_of_seg(fromLine, self.eps):
            if len(self.eps_neighborhood_of_seg(fromLine, self.eps)) > self.MinLns:
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
        
        cond1 = all(self.isDensityConneced(prod[0], prod[1]) for prod in itertools.product(clusterListOfSegments, clusterListOfSegments))
        if not cond1:
            return False
        
        cond2 =  all(prod[1] in clusterListOfSegments for
            prod in itertools.product(self.directReachablityGraph.nodes(), self.directReachablityGraph.nodes())
            if prod[0] in clusterListOfSegments and self.isDensityReachable(prod[1], prod[0]))
        return cond2
       
    
    
class segment(GLine):
    def __init__(self, gpt1, gpt2, stat = 0):
        GLine.__init__(self, gpt1, gpt2)
        
        ''' status values:
            0: unclassified
            1: signal
            2: noise
        '''
        self.status = stat    
    
