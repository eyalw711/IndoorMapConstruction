# -*- coding: utf-8 -*-
'''This file contains all the basic classes for geometry, geography and navigation'''
import traceback
import nvector as nv
from shapely.geometry import Polygon, MultiPoint, MultiLineString, LinearRing, LineString, MultiPolygon #later will use for polygons unification
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from shapely.geometry import Point as spt
import copy
import math
import random

class GPoint:
    """
    A GPoint is a Geographical point and a Geometric one
    param lat: latitude in degrees [90, -90]
    param long: longitude in degrees [0, 360]
    param alt: altitude in meters
    """
    def __init__(self, lat, long, alt = 0):
        self.geopyPoint = gpt(lat, long, alt)
        self.shapelyPoint = spt(long, lat, alt)
        
    def __repr__(self):
        return "(GPoint(lat: {}, long: {}, alt: {}))".format(self.geopyPoint[0],
                       self.geopyPoint[1], self.geopyPoint[2])
    def __getitem__(self, key):
        if key in [0,1,2]:
            return self.geopyPoint[key]
        else:
            raise IndexError("Invalid GPoint coord")
            
    def getLat(self):
        return self.geopyPoint[0]
    
    def getLong(self):
        return self.geopyPoint[1]
    
    def distance(self, other):
        '''returns the vincenty distance between the two points in meters
        Doesn't consider elevation'''
        
        if other == None:
            return -1
        else:
            frame_E = nv.FrameE(a=6371e3, f=0)
            positionA = frame_E.GeoPoint(latitude=self.getLat(), longitude=self.getLong(), degrees=False)
            positionB = frame_E.GeoPoint(latitude=other.getLat(), longitude=other.getLong(), degrees=False)
            path = nv.GeoPath(positionA, positionB)
            s_AB2 = path.track_distance(method='greatcircle').ravel()
            return s_AB2
#            return vincenty(self.geopyPoint,other.geopyPoint).meters
        
    def travel(self, bearing, distance):
        '''returns the destination of a travel from self for <distance> meters in <bearing> direction
        param bearing in deg
        param distance in meters'''
        nextPt = VincentyDistance(distance/1000).destination(self.geopyPoint, bearing)
        return GPoint(nextPt[0], nextPt[1], nextPt[2])
#        '''returns the destination of a travel from self for <distance> meters in <bearing> direction
#        param bearing in deg
#        param distance in meters'''
#        frame = nv.FrameE(a=6371e3, f=0)
#        pointA = frame.GeoPoint(latitude=self.getLat(), longitude=self.getLong(), degrees=True)
#        pointB, _azimuthb = pointA.geo_point(distance=dist, azimuth=bearing, degrees=True)
#        lat, lon = pointB.latitude_deg, pointB.longitude_deg
##        nextPt = VincentyDistance(distance/1000).destination(self.geopyPoint, bearing)
#        return GPoint(lat, lon, self[2])
        
    def noisy(self, sigma):
        bearingNoise = random.uniform(0, 360)
        distanceNoise = random.gauss(0, sigma)
        return self.travel(bearingNoise, distanceNoise)
    
    def within(self, obj):
        if type(obj) == Polygon:
            return self.shapelyPoint.within(obj)
        if type(obj) == Building:
            return self.within(obj.getPolygon())
        else:
            raise TypeError("expecting Polygon or Building but got {}".format(type(obj)))
    
    def toFix(self, time):
        return Fix(self.geopyPoint[0], self.geopyPoint[1], self.geopyPoint[2], time)
    
    def getCoords(self):
        return (self.shapelyPoint.x, self.shapelyPoint.y)

class GLine:
    ''' A line with two GPoints at its ends'''
    def __init__(self, gpt1, gpt2):
        # too early to use self functions:
        y1,x1 = gpt1.getCoords()
        y2,x2 = gpt2.getCoords()
        p1 = nv.GeoPoint(x1, y1)
        p2 = nv.GeoPoint(x2, y2)
        az = math.degrees(p1.distance_and_azimuth(p2, degrees = False)[1])
#        print("GLine ctor: type of az is {}".format(type(az)))
        if az < 0 or az > 180:
            gpt1, gpt2 = gpt2, gpt1
        self.end1 = gpt1
        self.end2 = gpt2
        
    def __repr__(self):
        return "(Gline ({},{}))".format(self[0], self[1])
        
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
        return self[0].distance(self[1])
    
    def minLatOfLineSeg(self):
        return min(self[0][0], self[1][0])
    
    def maxLatOfLineSeg(self):
        return max(self[0][0], self[1][0])
    
    def minLongOfLineSeg(self):
        return min(self[0][1], self[1][1])
    
    def maxLongOfLineSeg(self):
        return max(self[0][1], self[1][1])
    
    def isEqual(self, other):
        if self is other:
            return True
        
        me1x, me1y = self[0][0], self[0][1]
        oe1x, oe1y = other[0][0], other[0][1]
        
        me2x, me2y = self[1][0], self[1][1]
        oe2x, oe2y = other[1][0], other[1][1]
        
        #compare both ends both coords
        if me1x == oe1x and me1y == oe1y and me2x == oe2x and me2y == oe2y: 
            return True
        
        #compare backwards
        if me1x == oe2x and me1y == oe2y and me2x == oe1x and me2y == oe1y:
            return True
    
    def azimuth(self):  # azimuth between nvector GeoPoints
        '''returns azimuth in degrees.  
        Arbitrary in +-pi'''
        p1 = nv.GeoPoint(self.end1.getLat(), self.end1.getLong())
        p2 = nv.GeoPoint(self.end2.getLat(), self.end2.getLong())
        az = math.degrees(p1.distance_and_azimuth(p2, degrees = False)[1])
#        if math.isnan(az):
#            raise ValueError("azimuth is nan for linesegment {}".format(self))
        return az
            
    def myDistance(self, other):
        '''
        my distance function
        computes the sum of three lengths:
        start to start, middle to middle, end to end
        can argue this is a metric:
        * positive-definite
        * symmetric
        * complies with triangle inequality
        '''
        if self.isEqual(other):
            return 0
        
        # See example in https://pypi.python.org/pypi/nvector
        frame = nv.FrameE(a=6371e3, f= 0)
        
        myStart = frame.GeoPoint(self[0][0], self[0][1], degrees = True)
        myEnd = frame.GeoPoint(self[1][0], self[1][1], degrees = True)
        
        otherStart = frame.GeoPoint(other[0][0], other[0][1], degrees = True)
        otherEnd = frame.GeoPoint(other[1][0], other[1][1], degrees = True)
        
        #myMiddle
        myStartVector = myStart.to_nvector()
        myEndVector = myEnd.to_nvector()
        myPath = nv.GeoPath(myStartVector, myEndVector)
        myMiddle = myPath.interpolate(0.5).to_geo_point()
        #otherMiddle
        otherStartVector = otherStart.to_nvector()
        otherEndVector = otherEnd.to_nvector()
        otherPath = nv.GeoPath(otherStartVector, otherEndVector)
        otherMiddle = otherPath.interpolate(0.5).to_geo_point()
        
        startToStartPath = nv.GeoPath(myStart, otherStart)
        endToEndPath = nv.GeoPath(myEnd, otherEnd)
        
        try:
            intersectingPointObject = startToStartPath.intersect(endToEndPath) #.to_geo_point()
            
            if startToStartPath.on_path(intersectingPointObject) or endToEndPath.on_path(intersectingPointObject):
                # chosen points incorrectly for start, end
                otherStart, otherEnd = otherEnd, otherStart
                startToStartPath = nv.GeoPath(myStart, otherStart)
                endToEndPath = nv.GeoPath(myEnd, otherEnd)
        except RuntimeWarning:
            print("startToStart (x1,y1) = ({},{}) (x2,y2) = ({},{})".format(startToStartPath.positionA.latitude,
                  startToStartPath.positionA.longitude, startToStartPath.positionB.latitude,
                  startToStartPath.positionB.longitude))
            print("endToEndPath (x1,y1) = ({},{}) (x2,y2) = ({},{})".format(endToEndPath.positionA.latitude,
                  endToEndPath.positionA.longitude, endToEndPath.positionB.latitude,
                  endToEndPath.positionB.longitude))
            traceback.print_exc()
        
        def _dist(p1, p2): # distance between nvector GeoPoints
            return p1.distance_and_azimuth(p2)[0]
        
        return _dist(myStart, otherStart) +\
                _dist(myMiddle, otherMiddle) +\
                _dist(myEnd, otherEnd)
 

class Segment(GLine):
    def __init__(self, gpt1, gpt2, stat = 0, trajIndex = 0, clusterId = -1):
        GLine.__init__(self, gpt1, gpt2)
        
        ''' status values:
            0: unclassified
            1: signal
            2: noise
        '''
        self.status = stat    
        self.trajIndex = trajIndex


class LLine:
    '''
    in opposition to GLine, a LLine is local xy coordinates
    '''
    
#    def LLineTest():
#        gl = GLine(GPoint(0,0.0001), GPoint(-0.00011, 0.0002))
#        prj = EquirectangularProjector([gl])
#        ll = LLine(gl, prj)
#        print(ll, ll.length(), ll.azimuth())
        
    def __init__(self, xy1, xy2):
        x1, y1 = xy1
        x2, y2 = xy2
        self.end1 = (x1, y1)
        self.end2 = (x2, y2)
        if self.end1[0] > self.end2[0]:
            self.end1, self.end2 = self.end2, self.end1
#            print("Tried creating a LLine in reverse!")
            
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
    
    def minXOfLineSeg(self):
        return min(self[0][0], self[1][0])
    
    def maxXOfLineSeg(self):
        return max(self[0][0], self[1][0])
    
    def minYOfLineSeg(self):
        return min(self[0][1], self[1][1])
    
    def maxYOfLineSeg(self):
        return max(self[0][1], self[1][1])
    
    def isEqual(self, other):
        if self is other:
            return True
        
        me1x, me1y = self[0][0], self[0][1]
        oe1x, oe1y = other[0][0], other[0][1]
        
        me2x, me2y = self[1][0], self[1][1]
        oe2x, oe2y = other[1][0], other[1][1]
        
        #compare both ends both coords
        if me1x == oe1x and me1y == oe1y and me2x == oe2x and me2y == oe2y: 
            return True
        
        #compare backwards
        if me1x == oe2x and me1y == oe2y and me2x == oe1x and me2y == oe1y:
            return True
    
    def llineTravel(self, xytup, frac = 1, reverse = False):
        '''
        returns the xytup coords of the travel destination
        '''
        
        if len(xytup) != 2:
            raise ValueError("xytup should have 2 elements!")
        deltaX = (self.end2[0] - self.end1[0])*frac
        deltaY = (self.end2[1] - self.end1[1])*frac

        if reverse:
            deltaX, deltaY = -deltaX, -deltaY
        return (xytup[0] + deltaX, xytup[1] + deltaY)

    def myDistance(self, other):
        '''
        my distance function
        computes the sum of three lengths:
        start to start, middle to middle, end to end
        can argue this is a metric:
        * positive-definite
        * symmetric
        * complies with triangle inequality
        '''
        if self.isEqual(other):
            return 0
        
        # See example in https://pypi.python.org/pypi/nvector
#        frame = nv.FrameE(a=6371e3, f= 0)
        
        myStart = self.end1
        myEnd = self.end2
        
        otherStart = other.end1 
        otherEnd = other.end2 
        
        #myMiddle
        myMiddle = ((myStart[0] + myEnd[0])/2 , (myStart[1] + myEnd[1])/2)
        #otherMiddle
        otherMiddle = ((otherStart[0] + otherEnd[0])/2 , (otherStart[1] + otherEnd[1])/2)
        
#        startToStartPath = LLine(myStart, otherStart)
#        endToEndPath = LLine(myEnd, otherEnd)
        # check intersection between self, other
        # cehck intersection between startToStartPath, endToEndPath
        # deal with both cases
        
        def _dist(p1, p2):
            return math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2)

        factor = 5
        d1 = factor*_dist(myStart, otherStart) + _dist(myMiddle, otherMiddle) + factor*_dist(myEnd, otherEnd)
                
        d2 = factor*_dist(myStart, otherEnd) + _dist(myMiddle, otherMiddle) + factor*_dist(myEnd, otherStart)
        return min(d1, d2)

class LSegment(LLine):
    def __init__(self, xy1, xy2, stat, trajIndex):
        LLine.__init__(self,  xy1, xy2)
        
        ''' status values:
            0: unclassified
            1: signal
            2: noise
            3: medium
        '''
        self.status = stat    
        self.trajIndex = trajIndex
        
    
    
class Fix(GPoint):
    def __init__(self, lat, long, alt = 0, time = 0):
        GPoint.__init__(self, lat, long, alt)
        self.time = time
        
    def __repr__(self):
        return "Fix (lat: {}, long: {}, alt: {}, time: {})".format(self.geopyPoint[0],
                       self.geopyPoint[1], self.geopyPoint[2], self.time)
    
    def append(self, pts):
        '''returns a new trajectory'''
        if type(pts) == Fix:
            if self.time <= pts.time:
                return Trajectory([self,pts])
            else:
                raise Exception('Cannot append ' + pts + ' since times are inconsistent!')
        elif type(pts) == Trajectory:
            if len(pts.FixList) == 0:
                return Trajectory([self])
            else:
                if self.time <= pts.FixList[0].time:
                    return Trajectory([self]+pts.FixList)
                else:
                    raise Exception('Cannot append ' + pts + ' since times are inconsistent!')
        else:
            raise Exception('Expecting Fix or Trajectory and got {}'.format(type(pts)))

            
            
class Trajectory:
    def __init__(self, FixList):
        self.FixList = FixList
        self.traj = None
        self.trajValid = False
        
    def __repr__(self):
        return "aTrajectory with {} points".format(len(self.FixList))
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            return Trajectory(copy.deepcopy(self.FixList[key])) #returns new Trajectory with copy of Fixs
        return copy.deepcopy(self.FixList[key]) #returns a Fix
    
    def toLineSegmentList(trajectory):
        return [GLine(trajectory[i], trajectory[i+1]) for i in range(len(trajectory.FixList)-1)]
    
    def append(self, traj):
        '''appends in place'''
        if type(traj) == Fix:
            if all(FixListPoint.time <= traj.time for FixListPoint in self.FixList):
                self.FixList += [traj]
                self.trajValid = False
            else:
                raise Exception('Cannot append pt=' + traj + 'to Trajectory since times are inconsistent!')
                
        elif type(traj) == Trajectory:
            if len(self.FixList) == 0:
                self.FixList = traj.FixList
            else:
                if all(FixListPoint.time >= self.FixList[-1].time for FixListPoint in traj.FixList):
                    self.FixList += traj.FixList
                    self.trajValid = False
                else:
                    raise Exception('Cannot append trajectory=' + traj + 'to Trajectory since times are inconsistent!')
        else:
             raise Exception('Expecting Fix or Trajectory and got {}'.format(type(traj)))
             
    def scatterXY(self):
        x = [fx.shapelyPoint.x for fx in self.FixList]
        y = [fx.shapelyPoint.y for fx in self.FixList]
        return (x,y)

    def addNoise(self, sigma, dist = 3):
        newtraj = [p.noisy(dist) for p in self.FixList]
        self.traj = newtraj
        self.trajValid = True
    
    def plot(self, axs, color = '#6699cc', linestyle = '', linewidth = 3):
        xs, ys = self.scatterXY()
        if linestyle != '':
            axs.plot(xs, ys, color, alpha=0.7, linestyle = linestyle, linewidth=linewidth, solid_capstyle='round', zorder=2)
        else:
            axs.plot(xs, ys, color, alpha=0.7, linewidth=linewidth, solid_capstyle='round', zorder=2)
