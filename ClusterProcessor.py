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
        
#        fig, axs = plt.subplots()
        
        def pltLLines(axs, llinesList, color = 'r'):
            for lline in llinesList:
                xs = [lline.end1[0], lline.end2[0]]
                ys = [lline.end1[1], lline.end2[1]]
                axs.plot(xs, ys, color = color)
        
        
        avgDir = self.calcAverageDirection()
        avgDirEnd1, avgDirEnd2 = avgDir[0], avgDir[1]
#        print("calcRepresentativeLocal: calculated AvgDir")
        # need to find phi and rotate the projector:
        phi = ClusterProcessor.avgDirPhi(avgDirEnd1, avgDirEnd2)
#        print("calcRepresentativeLocal: phi = {}".format(phi))
        self.rotationalProjector.rotate(phi)
        
#        print("calcRepresentativeLocal: input - approx mean x val = {}, approx mean y val = {}".format(
#                    statistics.mean([lline.end1[0] for lline in self.localLLinesList]),
#                    statistics.mean([lline.end1[1] for lline in self.localLLinesList])))
        
#        pltLLines(axs, self.localLLinesList, color = 'b')
        
        rotLLinesList = [self.rotateLLineCoords(lline) for lline in self.localLLinesList]
        
#        pltLLines(axs, rotLLinesList, color = 'r')
        
#        ############### TESTTTTTTT ################333
#        print("############### TESTTTTTTT 1111111111111111", self.projector.rotatedXYToXY(160898298.15807286, 205006006.25297567))
        
        
#        print("calcRepresentativeLocal: after rotation - approx mean x val = {}, approx mean y val = {}".format(
#                    statistics.mean([lline.end1[0] for lline in rotLLinesList]),
#                    statistics.mean([lline.end1[1] for lline in rotLLinesList])))
#        
        rotLLinesDxYy = [(lline.end2[0] - lline.end1[0], lline.end2[1] - lline.end1[1]) for lline in rotLLinesList]
        
        # sorting llines once by starting points and once by ending points
        startingPoint = lambda l: l.end1[0]
        endingPoint = lambda l: l.end2[0]
        
        sortedByStartings = sorted(rotLLinesList, key = startingPoint)
        sortedByEndings = sorted(rotLLinesList, key = endingPoint)
        
#        print("sorted all segments")
        
        calculationsPool = [] #list of elements of format [x-val, [y-vals]]
        xval = float("-inf")
        while len(sortedByStartings) > 0 or len(sortedByEndings) > 0:
#            print("sortedByStartings: {}, sortedByEndings: {}".format(len(sortedByStartings), len(sortedByEndings)))
            if len(sortedByStartings) == 0 and len(sortedByEndings) > 0:
                a_starting_point = False
            elif len(sortedByEndings) == 0 and len(sortedByStartings ) > 0:
                a_starting_point = True
            else:
                if sortedByStartings[0].end1[0] == sortedByEndings[0].end1[0]:
                    sortedByEndings = sortedByEndings[1:]
                    a_starting_point = True
                else:
                    a_starting_point = sortedByStartings[0].end1[0] <= sortedByEndings[0].end1[0]
            
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
                if not (lline.end1[0] >= xval and xval <= lline.end2[0]):
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
                    
            if len(yvals) < self.MinLns:
                yvals = []
                continue
            
#            print(xval ,statistics.mean(yvals), "back to original", self.projector.rotatedXYToXY(xval, statistics.mean(yvals)))
            
            calculationsPool.append([xval, statistics.mean(yvals), (statistics.stdev(yvals))])
            yvals = []
        

#        self.projector.rotate(-phi)
#        print("############### TESTTTTTTT 222222222222222222", self.projector.xyToRotatedXY(160898298.15807286, 205006006.25297567))
        representative_and_walls = [(self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]),\
                                     self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]+ math.sqrt(3)*pe[2]),\
                                     self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]- math.sqrt(3)*pe[2])) for pe in calculationsPool]
        
        if len(representative_and_walls) < 2:
            return []
#        
#        representative_and_walls = [((pe[0], pe[1]),\
#                                     (pe[0], pe[1]+ math.sqrt(3)*pe[2]),\
#                                     (pe[0], pe[1]- math.sqrt(3)*pe[2])) for pe in calculationsPool]
#        
#        print("calcRepresentativeLocal output - approx mean x val = {}, approx mean y val = {}".format(
#                    statistics.mean([raw[0][0] for raw in representative_and_walls]),
#                    statistics.mean([raw[0][1] for raw in representative_and_walls])))
#        axs.scatter([e[0][0] for e in representative_and_walls], [e[0][1] for e in representative_and_walls], color = 'b')
#        axs.scatter([e[1][0] for e in representative_and_walls], [e[1][1] for e in representative_and_walls], color = 'r')
        return representative_and_walls
        
                    
        
    
#    def calcRepresentative(self, clusterListOfSegments, gamma):
#        raise Exception("dont use this function")
#        maxDist = 1000 #suppose no lines in cluster will have 1km perpendicular distance
#        
#        averageDirection = ClusterProcessor.calcAverageDirection(clusterListOfSegments)
#        projector = EquirectangularProjector(clusterListOfSegments)
#        phi = averageDirection.azimuth()
#        projector.rotate(phi)
#        
#        firstEnd = lambda x: x.end1
#        secndEnd = lambda x: x.end2
#        
#        getXtag = lambda p: projector.latLongToRotatedXY(p.getLat(), p.getLong())[0]
#        getYtag = lambda p: projector.latLongToRotatedXY(p.getLat(), p.getLong())[1]
#        
#        groupOfAllGPoints = list(chain.from_iterable((firstEnd(x), secndEnd(x)) for x in clusterListOfSegments))
#        k = lambda p: projector.latLongToRotatedXY(p.getLat(), p.getLong())[0]
#        sortedGroupOfAllGPoints = sorted(groupOfAllGPoints, key = k)
#        
#        shapelySortedPointsList_local = [shapelyPoint(getXtag(p), getYtag(p)) for p in sortedGroupOfAllGPoints]
#        shapelyLineStringSegments_local = GeometryCollection([LineString([(getXtag(firstEnd(seg)), getYtag(firstEnd(seg))), (getXtag(secndEnd(seg)), getYtag(secndEnd(seg)))]) \
#                                     for seg in clusterListOfSegments])
#        represntative_local = []
#        for p in shapelySortedPointsList_local:
#            sweeper_local = LineString([(p.x, p.y - maxDist), (p.x, p.y - maxDist)])
#            intersection_local = sweeper_local.intersection(shapelyLineStringSegments_local)
#            ys = []
#            for ob in intersection_local:
#                x, y = ob.xy
#                if len(y) == 1:
#                    ys.append(y)
#            represntative_local.append((p.x, sum(ys)/len(ys)))
#        
#        representative_latlong_tup = [projector.rotatedXYToLatLong(p[0],p[1]) for p in represntative_local]
#        representativeFixList = [Fix(tup[0],tup[1]) for tup in representative_latlong_tup]
#        return representativeFixList
    
#        representativeGLineList = []
#        for i in range(len(representative_latlong_tup) - 1):
#            llTup_s = representative_latlong_tup[i]
#            llTup_e = representative_latlong_tup[i+1]
#            gline = GLine(GPoint(llTup_s[0], llTup_s[1]),GPoint(llTup_e[0], llTup_e[1]))
#            representativeGLineList.append(gline)
#        return representativeGLineList
        
        
        
#        firstEnd = lambda x: x.end1
#        secndEnd = lambda x: x.end2
#        
#        groupOfAllGPoints = list(chain.from_iterable((firstEnd(x), secndEnd(x)) for x in clusterListOfSegments))
#        sortedGroupOfAllGPoints = sorted(groupOfAllGPoints)
#        
#        pToLatLong = lambda p: (p.getLat(), p.getLong())
#        
#        for p in sortedGroupOfAllGPoints:
#            cntr = 0
#            for seg in clusterListOfSegments:
#                e1_lat, e1_long = pToLatLong(seg.end1)
#                e1_xtag, e1_ytag = projector.latLongToRotatedXY(e1_lat, e1_long)
#                
#                e2_lat, e2_long = pToLatLong(seg.end2)
#                e2_xtag, e2_ytag = projector.latLongToRotatedXY(e2_lat, e2_long)
#                
#                p_lat, p_long = pToLatLong(p)
#                p_xtag, p_ytag = projector.latLongToRotatedXY(p_lat, p_long)
#                
#                if e1_xtag <= p_xtag and p_xtag <= e2_xtag:
#                    cntr+=1
#                
#                
        
        
    