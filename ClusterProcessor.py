# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 00:38:50 2017

@author: Eyal
"""
from shapely.geometry import LineString, GeometryCollection
from shapely.geometry import  Point as shapelyPoint
from shapely.ops import cascaded_union
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from IMCObjects import GPoint, Fix, Trajectory, GLine
from Projector import EquirectangularProjector
import copy
from itertools import chain
import random
import math

class LLine:
    def LLineTest():
        gl = GLine(GPoint(0,0.0001), GPoint(-0.00011, 0.0002))
        prj = EquirectangularProjector([gl])
        ll = LLine(gl, prj)
        print(ll, ll.length(), ll.azimuth())
        
    def __init__(self, gline, projector, coords = None):
        if coords == None:
            xe1, ye1 = projector.latLongToXY(gline.end1.getLat(), gline.end1.getLong())
            self.end1 = (xe1, ye1)
            
            xe2, ye2 = projector.latLongToXY(gline.end2.getLat(), gline.end2.getLong())
            self.end2 = (xe2, ye2)
            if xe2 - xe1 < 0:
                self.end1, self.end2 = self.end2, self.end1
        else:
            self.end1 = coords[0]
            self.end2 = coords[1]
            
    def __repr__(self):
        return "(LLine ({},{}))".format(self[0], self[1])
        
    def __getitem__(self, key):
        if key == 0:
            return self.end1
        elif key == 1:
            return self.end2
        else:
            raise IndexError("GLine has only two indices!")
    
    def length(self):
        if self[0][0] == self[1][0] and self[0][1] == self[1][1]:
            return 0
        return math.sqrt((self.end2[0] - self.end1[0])**2 + (self.end2[1] - self.end1[1])**2)
    
    def azimuth(self):
        deltaX = self.end2[0] - self.end1[0]
        if deltaX < 0:
            self.end1, self.end2 = self.end2, self.end1
            return self.azimuth()
        deltaY = self.end2[1] - self.end1[1]
        if deltaY == 0: #az 90
            return 90
        elif deltaY > 0: #az [0,90):
            return math.degrees(math.atan(deltaX/deltaY))
        else:   #az (90, 180)
            complemet = math.degrees(math.atan(deltaX/(-deltaY)))
            return 180 - complemet
    
    def llineTravel(self, xytup, frac = 1):
        '''
        returns the xytup coords of the travel destination
        '''
        
        if len(xytup) != 2:
            raise ValueError("xytup should have 2 elemets!")
        deltaX = (self.end2[0] - self.end1[0])*frac
        deltaY = (self.end2[1] - self.end1[1])*frac
        return (xytup[0] + deltaX, xytup[1] + deltaY)

class ClusterProcessor:
    ''' works each time on a different single cluster'''
    def __init__(self, projector):
        self.projector = projector
        self.localLLinesList = None
        
    def loadGeoCluster_segmentsList(self, geoClusterGlineList):
        self.localLLinesList = []
        for glineseg in geoClusterGlineList:
            self.localLLinesList.append(LLine(glineseg, self.projector))
            
        
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
            glineAz, glineLength = lline.azimuth(), lline.length()
            print("glineAz = {}, glineLength = {}".format(glineAz, glineLength))
            nextPosition = lline.llineTravel(currPosition)
            currPosition = nextPosition
            print("calcAverageDirection: currPosition = {}".format(currPosition))
        
        totalWay = LLine(None, None, coords = [origin,currPosition])
#        avgDirAzimuth = totalWay.azimuth()
#        avgDirLength = totalWay.length() / len(self.localLLinesList)
        
        originCopy = copy.deepcopy(origin)
        avgDirEnd2 = totalWay.llineTravel(originCopy, 2/len(self.localLLinesList))
        
        avgDir = LLine(None, None, coords = [originCopy, avgDirEnd2])
        print("plotClusters: avgDir: {}".format(avgDir))
        return avgDir.end1, avgDir.end2
    
    def calcRepresentative(self, clusterListOfSegments, gamma):
        
        maxDist = 1000 #suppose no lines in cluster will have 1km perpendicular distance
        
        averageDirection = ClusterProcessor.calcAverageDirection(clusterListOfSegments)
        projector = EquirectangularProjector(clusterListOfSegments)
        phi = averageDirection.azimuth()
        projector.rotate(phi)
        
        firstEnd = lambda x: x.end1
        secndEnd = lambda x: x.end2
        
        getXtag = lambda p: projector.latLongToRotatedXY(p.getLat(), p.getLong())[0]
        getYtag = lambda p: projector.latLongToRotatedXY(p.getLat(), p.getLong())[1]
        
        groupOfAllGPoints = list(chain.from_iterable((firstEnd(x), secndEnd(x)) for x in clusterListOfSegments))
        k = lambda p: projector.latLongToRotatedXY(p.getLat(), p.getLong())[0]
        sortedGroupOfAllGPoints = sorted(groupOfAllGPoints, key = k)
        
        shapelySortedPointsList_local = [shapelyPoint(getXtag(p), getYtag(p)) for p in sortedGroupOfAllGPoints]
        shapelyLineStringSegments_local = GeometryCollection([LineString([(getXtag(firstEnd(seg)), getYtag(firstEnd(seg))), (getXtag(secndEnd(seg)), getYtag(secndEnd(seg)))]) \
                                     for seg in clusterListOfSegments])
        represntative_local = []
        for p in shapelySortedPointsList_local:
            sweeper_local = LineString([(p.x, p.y - maxDist), (p.x, p.y - maxDist)])
            intersection_local = sweeper_local.intersection(shapelyLineStringSegments_local)
            ys = []
            for ob in intersection_local:
                x, y = ob.xy
                if len(y) == 1:
                    ys.append(y)
            represntative_local.append((p.x, sum(ys)/len(ys)))
        
        representative_latlong_tup = [projector.rotatedXYToLatLong(p[0],p[1]) for p in represntative_local]
        representativeFixList = [Fix(tup[0],tup[1]) for tup in representative_latlong_tup]
        return representativeFixList
    
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
        
        
    