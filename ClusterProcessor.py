# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 00:38:50 2017

@author: Eyal
"""

from matplotlib import pyplot as plt
from shapely.geometry import LineString, GeometryCollection
from shapely.geometry import  Point as shapelyPoint
#from shapely.ops import cascaded_union
#from geopy.distance import vincenty, VincentyDistance
#from geopy.distance import Point as gpt
from IMCObjects import LLine
from Projector import RotationalProjector
import copy
from itertools import chain
import math
import statistics



class ClusterProcessor:
    '''
    works each time on a different single cluster
    works in xy-local
    '''
    def __init__(self, MinLns):
        self.localLLinesList = None
        self.MinLns = MinLns
        self.rotationalProjector = RotationalProjector()
        
    def loadLocalCluster_lsegList(self, localClusterLlineList):
        self.localLLinesList = localClusterLlineList
            
    def calcAverageDirection(self):
        '''
        assume all LLines are with bearing 0 <= bearing < 180
        returns xytup_start xytup_end
        '''
        if len(self.localLLinesList) < 1:
            raise ValueError("calcMainDirection: cluster has to have segments!")
        origin = self.localLLinesList[0].end1 #(x,y) tuple
        
        currPosition = origin
        for lline in self.localLLinesList:
            nextPosition = lline.llineTravel(currPosition)
            currPosition = nextPosition
        
        totalWay = LLine(origin,currPosition)
        
        originCopy = copy.deepcopy(origin)
        avgDirEnd2 = totalWay.llineTravel(originCopy, 2/len(self.localLLinesList))
        
        avgDir = LLine(originCopy, avgDirEnd2)
        return avgDir
    
    def avgDirPhi(avgDirEnd1, avgDirEnd2):
        dx = avgDirEnd2[0] - avgDirEnd1[0]
        dy = avgDirEnd2[1] - avgDirEnd1[1]
        if dx == 0 and dy == 0:
            raise Exception("calcRepresentativeLocal: main dir of cluster is a dot!")
        elif dx == 0 and dy > 0:
            phi = 90
        elif dx == 0 and dy < 0:
            phi = -90
        elif dy == 0:
            phi = 0
        elif dy > 0:
            phi = math.degrees(math.atan(dy/dx))
        else:
            phi = - math.degrees(math.atan(-dy/dx))
        return phi
    
    
    def rotateLLineCoords(self, lline):
        e1x, e1y = lline.end1
        e1xytag = self.rotationalProjector.xyToRotatedXY(e1x, e1y)
        
        e2x, e2y = lline.end2
        e2xytag = self.rotationalProjector.xyToRotatedXY(e2x, e2y)
        return LLine(e1xytag, e2xytag)
    
    
    def calcRepresentativeLocal(self, gamma):
        '''
        calculates the representative trajectory for the cluster
        currently loaded to the processor
        '''
        
        def pltLLines(axs, llinesList, color = 'r'):
            for lline in llinesList:
                xs = [lline.end1[0], lline.end2[0]]
                ys = [lline.end1[1], lline.end2[1]]
                axs.plot(xs, ys, color = color)
        
        
        avgDir = self.calcAverageDirection()
        avgDirEnd1, avgDirEnd2 = avgDir[0], avgDir[1]

        phi = ClusterProcessor.avgDirPhi(avgDirEnd1, avgDirEnd2)
        self.rotationalProjector.rotate(phi)
               
        rotLLinesList = [self.rotateLLineCoords(lline) for lline in self.localLLinesList]
        rotLLinesDxYy = [(lline.end2[0] - lline.end1[0], lline.end2[1] - lline.end1[1]) for lline in rotLLinesList]
        
        # sorting llines once by starting points and once by ending points
        startingPoint = lambda l: l.end1[0]
        endingPoint = lambda l: l.end2[0]
        
        sortedByStartings = sorted(rotLLinesList, key = startingPoint)
        sortedByEndings = sorted(rotLLinesList, key = endingPoint)
        
        calculationsPool = [] #list of elements of format [x-val, [y-vals]]
        xval = float("-inf")
        while len(sortedByStartings) > 0 or len(sortedByEndings) > 0:
            if len(sortedByStartings) == 0 and len(sortedByEndings) > 0:
                a_starting_point = False
            elif len(sortedByEndings) == 0 and len(sortedByStartings ) > 0:
                a_starting_point = True
            else:
                if sortedByStartings[0].end1[0] == sortedByEndings[0].end2[0]:
                    sortedByEndings = sortedByEndings[1:]
                    a_starting_point = True
                else:
                    a_starting_point = sortedByStartings[0].end1[0] <= sortedByEndings[0].end2[0]
            
            if a_starting_point:
                currLLine = sortedByStartings[0]
                sortedByStartings = sortedByStartings[1:]
                currPt = currLLine.end1
            else:
                currLLine = sortedByEndings[0]
                sortedByEndings = sortedByEndings[1:]
                currPt = currLLine.end2
            
            if currPt[0] < xval + gamma:
                continue
            
            xval = currPt[0]
            yvals = []
            # add to not look at too many segments
            for i, lline in enumerate(rotLLinesList):
                if not (lline.end1[0] <= xval and xval <= lline.end2[0]):
                    continue
                
                dx, dy = rotLLinesDxYy[i]
                if dx == 0 and lline.end1[0] == xval:
                    # edge case exactly perpendicular to axis
                    yvals.append((lline.end1[1] + lline.end2[1]) / 2) 
                elif dx != 0:
                    # need to find the intersection with the line.
                    m = dy/dx
                    u,v = lline.end2
                    y = m*(xval - u) + v
                    yvals.append(y)
                    
            if len(yvals) < 2: #self.MinLns:
                continue
            
            calculationsPool.append([xval, statistics.mean(yvals), (statistics.stdev(yvals))])
            
        representative_and_walls = [(self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]),\
                                     self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]+ math.sqrt(3)*pe[2]),\
                                     self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]- math.sqrt(3)*pe[2])) for pe in calculationsPool]
             
        if len(representative_and_walls) < 2:
            return []
        return representative_and_walls