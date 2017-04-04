# -*- coding: utf-8 -*-
'''This file contains all the basic classes for geometry, geography and navigation'''

import nvector as nv
from shapely.geometry import Polygon, MultiPoint, MultiLineString, LinearRing, LineString, MultiPolygon #later will use for polygons unification
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from shapely.geometry import Point as spt
import copy
from math import sin, cos, acos


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
        return "GPoint (lat: {}, long: {}, alt: {})".format(self.geopyPoint[0],
                       self.geopyPoint[1], self.geopyPoint[2])
    def __getitem__(self, key):
        if key in [0,1,2]:
            return self.geopyPoint[key]
        else:
            raise IndexError("Invalid GPoint coord")
            
    def distance(self, other):
        '''returns the vincenty distance between the two points in meters
        Doesn't consider elevation'''
        if other == None:
            return -1
        else:
            return vincenty(self.geopyPoint,other.geopyPoint).meters
        
    def travel(self, bearing, distance):
        '''returns the destination of a travel from self for <distance> meters in <bearing> direction
        param bearing in deg
        param distance in meters'''
        nextPt = VincentyDistance(distance/1000).destination(self.geopyPoint, bearing)
        return GPoint(nextPt[0], nextPt[1], nextPt[2])
        
    def noisy(self, sigma):
        bearingNoise = random.uniform(0, 360)
        disanceNoise = random.gauss(0, sigma)
        return self.travel(bearingNoise, disanceNoise)
    
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
        self.end1 = gpt1
        self.end2 = gpt2
        
    def __getitem__(self, key):
        if key == 0:
            return self.end1
        elif key == 1:
            return self.end2
        else:
            raise IndexError("GLine has only two indices!")
    
    def length(self):
        return self[0].distance(self[1])
        
    def perpendicularDist(line1, line2):
        '''
        Perpendicular distance function. See Article.
        In this function L1 is being Li and L2 is Lj, meaning
        the projection points are on L1.
        
        returns value in meters
        '''
        eps = 1e-2
        
        # See example in https://pypi.python.org/pypi/nvector
        frame = nv.FrameE(a=6371e3, f= 0)
        pointA1 = frame.GeoPoint(line1[0][0], line1[0][1], degrees = True)#TODO: hoping GeoPoint is in lat long alt
        pointA2 = frame.GeoPoint(line1[1][0], line1[1][1], degrees = True)
        
        pointB1 = frame.GeoPoint(line2[0][0], line2[0][1], degrees=True)
        pointB2 = frame.GeoPoint(line2[1][0], line2[1][1], degrees=True)
        
        pathA = nv.GeoPath(pointA1, pointA2)
        s_xt1 = abs(pathA.cross_track_distance(pointB1, method = 'greatcircle'))#.ravel()
        s_xt2 = abs(pathA.cross_track_distance(pointB2, method = 'greatcircle'))#.ravel()
        
        if s_xt1 <= eps and s_xt2 <= eps:
            return 0.0
        else:
            d_perpen = (s_xt1**2 + s_xt2**2) / (s_xt1 + s_xt2)
            return d_perpen
    
    def angularDist(line1, line2):
        frame = nv.FrameE(name='WGS84')#a=6371e3, f= 0)
        pointA1 = frame.GeoPoint(line1[0][0], line1[0][1], degrees = True)#TODO: hoping GeoPoint is in lat long alt
        pointA2 = frame.GeoPoint(line1[1][0], line1[1][1], degrees = True)
        
        pointB1 = frame.GeoPoint(line2[0][0], line2[0][1], degrees=True)
        pointB2 = frame.GeoPoint(line2[1][0], line2[1][1], degrees=True)
        
        def _azimuth(p1, p2):
            '''returns azimuth in radians'''
            return p1.distance_and_azimuth(p2, degrees = False)[1]
        
        def _almostEquals(az1, az2):
            '''gets angles in radians'''
            eps = 1e-5
#            az1 = radians(az1)
#            az2 = radians(az2)
            deltaDist = (cos(az1) - cos(az2))**2 + (sin(az1) - sin(az2))**2
            angleDiff = acos((2.0 - deltaDist) / 2.0)
            if (abs(angleDiff) < eps):
                return True
            else:
                return False
                
        #if lines are almost the same bearing - return automatic 0
        azA = _azimuth(pointA1, pointA2)
        azB = _azimuth(pointB1, pointB2)
        
        if (_almostEquals(azA, azB)):
            return 0.0
        
        pathA = nv.GeoPath(pointA1, pointA2)
        pathB = nv.GeoPath(pointB1, pointB2)
        
        intersectingPointObject = pathA.intersection(pathB).to_geo_point()
        pointC = frame.GeoPoint(intersectingPointObject.latitude_deg[0], intersectingPointObject.longitude_deg[0], degrees = True)
#        print("pointC was made with", pointC.latitude, pointC.longitude, "but in deg is", pointC.latitude_deg, pointC.longitude_deg)
        #in Radians
        minAngle = min(abs(_azimuth(pointC, pointA1) - _azimuth(pointC, pointB1)),
                       abs(_azimuth(pointC, pointB1) - _azimuth(pointC, pointA2)),
                       abs(_azimuth(pointC, pointA2) - _azimuth(pointC, pointB2)),
                       abs(_azimuth(pointC, pointB2) - _azimuth(pointC, pointA1)))
        
        return line1.length() * sin(minAngle)
    
    def parallelDist(line1, line2):
        '''
        Parallel distance function. See Article.
        In this function L1 is being Li and L2 is Lj, meaning
        the projection points are on L1.
        
        returns value in meters
        '''
        # See example in https://pypi.python.org/pypi/nvector
        frame = nv.FrameE(name = 'WGS84', a=6371e3, f= 0)
        
        pointA1 = frame.GeoPoint(latitude = line1[0][0], longitude = line1[0][1], degrees = True)#TODO: hoping GeoPoint is in lat long alt
        pointA2 = frame.GeoPoint(latitude = line1[1][0], longitude = line1[1][1], degrees = True)
        pathA = nv.GeoPath(pointA1, pointA2)
        
        pointB1 = frame.GeoPoint(line2[0][0], line2[0][1], degrees=True)
        pointB2 = frame.GeoPoint(line2[1][0], line2[1][1], degrees=True)
        
        #projection points:
        pointC1 = pathA.closest_point_on_great_circle(pointB1)
        pointC2 = pathA.closest_point_on_great_circle(pointB2)
        
        def _dist(p1, p2):
            return p1.distance_and_azimuth(p2)[0]
        
        #not like in the reference since here we don't care who's shorter.
        return min([_dist(pointA1, pointC1), _dist(pointA1, pointC2),
                     _dist(pointA2, pointC1), _dist(pointA2, pointC2)])
 
    
    
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
    
    def plot(self, axs, color = '#6699cc', linestyle = ''):
        xs, ys = self.scatterXY()
        if linestyle != '':
            axs.plot(xs, ys, color, alpha=0.7, linestyle = linestyle, linewidth=3, solid_capstyle='round', zorder=2)
        else:
            axs.plot(xs, ys, color, alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
