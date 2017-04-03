# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:10:18 2017

@author: Eyal
""" 
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
import nvector as nv
from shapely.geometry import Polygon, MultiPoint, MultiLineString, LinearRing, LineString, MultiPolygon
from shapely.geometry import Point as spt
from matplotlib import pyplot as plt
import random
import pandas as pd
import copy
from Convertions import LLHtoECEF

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
            
    def perpendicularDist(line1, line2):
        '''
        Perpendicular distance function. See Article.
        In this function L1 is being Li and L2 is Lj, meaning
        the projection points are on L1.
        
        returns value in meters
        '''
        
        # See example in https://pypi.python.org/pypi/nvector
        frame = nv.FrameE(a=6371e3, f= 0)
        pointA1 = frame.GeoPoint(line1[0][0], line1[0][1], degrees = True)#TODO: hoping GeoPoint is in lat long alt
        pointA2 = frame.GeoPoint(line1[1][0], line1[1][1], degrees = True)
        
        pointB1 = frame.GeoPoint(line2[0][0], line2[0][1], degrees=True)
        pointB2 = frame.GeoPoint(line2[1][0], line2[1][1], degrees=True)
        
        pathA = nv.GeoPath(pointA1, pointA2)
        s_xt1 = pathA.cross_track_distance(pointB1, method = 'greatcircle')#.ravel()
        s_xt2 = pathA.cross_track_distance(pointB2, method = 'greatcircle')#.ravel()
                
        d_perpen = (s_xt1**2 + s_xt2**2) / (s_xt1 + s_xt2)
        return d_perpen
    
    def angularDist(line1, line2):
        return 0
    
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
        
        #scenario like article:
        if pathA.on_path(pointC1) and pathA.on_path(pointC2):
            #find out which point is closer to which end:
            if pointC1.distance_and_azimuth(pointA1)[0] >= pointC2.distance_and_azimuth(pointA1)[0]:
                pointC1tag = pointC2
                pointC2tag = pointC1
            else:
                pointC1tag = pointC1
                pointC2tag = pointC2
            return min(pointC1tag.distance_and_azimuth(pointA1)[0],
                       pointC2tag.distance_and_azimuth(pointA2)[0])
        #scenario like notebook
        else:
            #case a is closer to c, and b is closer to d:
            if (pointC1.distance_and_azimuth(pointA1)[0] < pointC2.distance_and_azimuth(pointA1)[0]
                and pointC2.distance_and_azimuth(pointA2)[0] < pointC1.distance_and_azimuth(pointA2)):
                #staright on
                pass
            elif (pointC1.distance_and_azimuth(pointA1)[0] < pointC2.distance_and_azimuth(pointA1)[0]
                and pointC2.distance_and_azimuth(pointA2)[0] < pointC1.distance_and_azimuth(pointA2)):
                #swap
                pass
            else: #case a is closer to both c,d than b is:
                #give b its closer point and a gets what's left
                pass
        return 0
    
    
class Building:
    def __init__(self, outerWallsLinearRing, innerWallsMultiLineString = None, holesMultiPolygon = None):
        """
        param outerWallsRing: a LinearRing class object for outer walls. include points for walls connection.
        param innerWallsList: a list of LineString class objects for inner walls of the building.
        """
        self.outerWallsLinearRing = outerWallsLinearRing
        self.innerWallsMultiLineString = innerWallsMultiLineString
        self.holesMultiPolygon = holesMultiPolygon
    
    def getPolygon(self):
        outerWallsCoords = self.outerWallsLinearRing.coords
        polygon = Polygon(outerWallsCoords)
        if not self.holesMultiPolygon == None:
            polygon = polygon.difference(self.holesMultiPolygon)
        return polygon
    
    def containsPoint(self, gpoint):
        polygon = self.getPolygon()
        return gpoint.within(polygon)
        
    def legalTravel(self, pt1, pt2):
        line = LineString([pt1.getCoords(), pt2.getCoords()])
        illegal = line.intersects(self.outerWallsLinearRing) or any(line.intersects(geo) for geo in [self.innerWallsMultiLineString, self.holesMultiPolygon] if geo != None)
        return not illegal

    def buildingBuilder(typenum):
        if typenum == 1:
            outerWallsCrds = list()
            meetingInnerOuter = list()
            startPoint = GPoint(32.178074, 34.911721)
            currPoint = GPoint(32.178074, 34.911721) #Start as north westren point
            outerWallsCrds.append(currPoint.getCoords())
            for i in range(5):
                nextPoint = currPoint.travel(90,5)
                outerWallsCrds.append(nextPoint.getCoords())
                if not (i == 4):
                    meetingInnerOuter.append(nextPoint)
                currPoint = nextPoint
                
            nextPoint = currPoint.travel(180, 10)
            outerWallsCrds.append(nextPoint.getCoords())
            currPoint = nextPoint
            
            nextPoint = currPoint.travel(270, 25)
            outerWallsCrds.append(nextPoint.getCoords())
            outerWallsCrds.append(startPoint.getCoords())
            
            outerWallsLinRing = LinearRing(outerWallsCrds)
            
            innerWallsCrds = list()
            for i in [0,2]:
                innerWallsCrds.append(
                        [meetingInnerOuter[i].getCoords(),
                         meetingInnerOuter[i].travel(180,8).getCoords(),
                         meetingInnerOuter[i+1].travel(180,8).getCoords(),
                         meetingInnerOuter[i+1].getCoords()])
            innerWallsMultiPolygon = MultiPolygon([Polygon(crds) for crds in innerWallsCrds])
            
            return Building(outerWallsLinRing, holesMultiPolygon = innerWallsMultiPolygon)

        else:
            print("No implementation for typenum = {}".format(typenum))
            return None
        
    def plot(self, axs):
        polygon = self.getPolygon()#Polygon(self.outerWallsLinearRing)
        xs, ys = polygon.exterior.xy
        axs.fill(xs, ys, alpha=0.5, fc='r', ec='none')
        if self.innerWallsMultiLineString != None:
            for line in self.innerWallsMultiLineString:
                x, y = line.xy
                axs.plot(x, y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)
    
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
            return Trajectory(copy.deepcopy(self.FixList[key])) #returns copy
        return copy.deepcopy(self.FixList[key])
    
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
    
    def plot(self, axs):
        xs, ys = self.scatterXY()
        axs.plot(xs, ys, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)


class TrajectoryMaker:
    def __init__(self, building):
        self.building = building
        self.trajectoryCollection = []
        
    def selectStartingFix(self):
        bnds = self.building.getPolygon().bounds
        while True: 
            long = random.uniform(bnds[0], bnds[2]) #easting
            lat = random.uniform(bnds[1],bnds[3]) #northing
            aPnt = Fix(lat,long)
            if aPnt.within(self.building.getPolygon()):
               break
        return aPnt
    
    def makeFixList(self, dt):
        '''returns a trajectory of Fixes inside the TrajectoryMaker's Polygon'''
        #walk 5km/h = 1.389 m/s
        velocity = 1.389
        #walk for 2 minutes
        T = 120
        
        t = 0
        
        currFix = self.selectStartingFix()
        aTraj = Trajectory([currFix])
        currBrng = random.uniform(0,360)
        while t < T:
            while True:
                bearingNoise = random.uniform(-25,25)
                nextFix = currFix.travel(currBrng + bearingNoise, velocity * dt).toFix(t+dt)
                if self.building.legalTravel(currFix, nextFix):
                #if nextFix.within(self.building):
                    t += dt
                    aTraj.append(nextFix)
                    currFix = nextFix
                    break
                else:
                    currBrng += random.uniform(-45,45)
        self.trajectoryCollection.append(aTraj)
        return aTraj
    
    def makeDataSet(self, filename, numberOfTrajectories):
        data = []
        for i in range(numberOfTrajectories):
            self.makeFixList(2)
        
        for i in range(len(self.trajectoryCollection)):
            traj = self.trajectoryCollection[i]
            for fix in traj.FixList:
                data.append({"trajIndex": i,
                             "time": fix.time,
                             "lat": fix.geopyPoint[0],
                             "long": fix.geopyPoint[1],
                             "alt": fix.geopyPoint[2]})
        df = pd.DataFrame(data)
        df.to_csv('Data//' + filename +'.csv')

                
class TrajectoryCollectionCSVLoader:
    def __init__(self):
        '''each element in the collection is of a Trajectory class object'''
        self.trajectoryCollection = {}
    
    def loadTrajectoryCollectionFromCSV(self, filename):
        data = pd.read_csv('data//' + filename +'.csv' )
        for index, row in data.iterrows():
            trajectoryIndex = row['trajIndex']
            if not trajectoryIndex in self.trajectoryCollection:
                self.trajectoryCollection[trajectoryIndex] = Trajectory(list())
            self.trajectoryCollection[trajectoryIndex].append(Fix(row['lat'], row['long'], row['alt'], row['time']))

#        for key, traj in self.trajectoryCollection.items():
#            print("trajIndex:{}, {}".format(key, traj))
                    
'''unused'''
def plotPoly(p):
    x,y = p.exterior.xy
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    ax.plot(x, y, color='#6699cc', alpha=0.7,
    linewidth=3, solid_capstyle='round', zorder=2)
    ax.set_title('Polygon')


def test1():
    bldng1 = [(32.178028, 34.912227),
            (32.178022, 34.912530),
            (32.177977, 34.912527),
            (32.177972, 34.912886),
            (32.177927, 34.912889),
            (32.177929, 34.912444),
            (32.177976, 34.912443),
            (32.177985, 34.912226)]
    
    
    ks_p1 = (32.177916, 34.911915)
    ks_p2 = (32.177925, 34.911857)
    d12 = vincenty(ks_p1,ks_p2).meters
    ptk_p1 = (32.101759, 34.850336)
    ks_p3 = VincentyDistance(200/1000).destination(ks_p1, 90)
    
    #pol1 = Polygon(bldng1)
    pol1 = Polygon([(lon,lat) for (lat,lon) in bldng1])
    pol2 = MultiPoint(bldng1).convex_hull    
    
    #plotPoly(pol1)
    #plot(pol2)
    
    tm = TrajectoryMaker(pol1)
    
    fig, axs = plt.subplots()
    xs, ys = pol1.exterior.xy
    axs.fill(xs, ys, alpha=0.5, fc='r', ec='none')
    
    traj = tm.makeFixList(2)
    x,y = traj.scatterXY()
    axs.plot(x,y)
    
    tm.makeDataSet('BarIlan', 20)

fig, axs = plt.subplots()  
bd = Building.buildingBuilder(1)
tm = TrajectoryMaker(bd)
t1 = tm.makeFixList(2)
bd.plot(axs)
t1.plot(axs)
tm.makeDataSet('KfarSaba', 20)
#frm = nv.FrameE(name = 'WGS84', a=6371e3, f= 0)
#pta1 = frm.GeoPoint(32,34, degrees = True)
#pta2 = frm.GeoPoint(31.999, 34.001, degrees = True)
#d, b1, b2 = pta1.distance_and_azimuth(pta2, degrees = True)
#print("dist {}m, bearing1 {}, bearing2 {}".format(d,b1,b2))

