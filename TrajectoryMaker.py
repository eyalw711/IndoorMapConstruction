# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:10:18 2017

@author: Eyal
""" 
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from shapely.geometry import Polygon, MultiPoint
from shapely.geometry import Point as spt
from matplotlib import pyplot as plt
import random
import pandas as pd

class GPoint:
    def __init__(self, lat, long, alt = 0):
        self.geopyPoint = gpt(lat, long, alt)
        self.shapelyPoint = spt(long, lat, alt)
        
    def __repr__(self):
        return "GPoint (lat: {}, long: {}, alt: {}, time: {})".format(self.geopyPoint[0],
                       self.geopyPoint[1], self.geopyPoint[2])
        
    def distance(self, other):
        '''returns the vincenty distance between the two points in meters
        Doesn't consider elevation'''
        if other == None:
            return -1
        else:
            return vincenty(self.geopyPoint,other.geopyPoint).meters
        
    def next_GPoint(self, bearing, distance):
        '''returns the destination of a travel from self for <distance> meters in <bearing> direction
        param bearing in deg
        param distance in meters'''
        nextPt = VincentyDistance(distance/1000).destination(self.geopyPoint, bearing)
        return GPoint(nextPt[0], nextPt[1], nextPt[2])
        
    def noisy(self, sigma):
        bearingNoise = random.uniform(0, 360)
        disanceNoise = random.gauss(0, sigma)
        return self.next_GPoint(bearingNoise, disanceNoise)
    
    def within(self, poly):
        return self.shapelyPoint.within(poly)
    
    def toFix(self, time):
        return Fix(self.geopyPoint[0], self.geopyPoint[1], self.geopyPoint[2], time)
   

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
            if len(pts.groundTruth) == 0:
                return Trajectory([self])
            else:
                if self.time <= pts.groundTruth[0].time:
                    return Trajectory([self]+pts.groundTruth)
                else:
                    raise Exception('Cannot append ' + pts + ' since times are inconsistent!')
        else:
            raise Exception('Expecting Fix or Trajectory and got {}'.format(type(pts)))

            
            
class Trajectory:
    def __init__(self, FixList):
        self.groundTruth = FixList
        self.traj = None
        self.trajValid = False
        
    def __repr__(self):
        return "aTrajectory with {} points".format(len(self.groundTruth))
    
    def append(self, traj):
        '''appends in place'''
        if type(traj) == Fix:
            if all(groundTruthPoint.time <= traj.time for groundTruthPoint in self.groundTruth):
                self.groundTruth += [traj]
                self.trajValid = False
            else:
                raise Exception('Cannot append pt=' + traj + 'to Trajectory since times are inconsistent!')
                
        elif type(traj) == Trajectory:
            if len(self.groundTruth) == 0:
                self.groundTruth = traj.groundTruth
            else:
                if all(groundTruthPoint.time >= self.groundTruth[-1].time for groundTruthPoint in traj.groundTruth):
                    self.groundTruth += traj.groundTruth
                    self.trajValid = False
                else:
                    raise Exception('Cannot append trajectory=' + traj + 'to Trajectory since times are inconsistent!')
        else:
             raise Exception('Expecting Fix or Trajectory and got {}'.format(type(traj)))
             
    def scatterXY(self):
        x = [fx.shapelyPoint.x for fx in self.groundTruth]
        y = [fx.shapelyPoint.y for fx in self.groundTruth]
        return (x,y)

    def addNoise(self, sigma, dist = 3):
        newtraj = [p.noisy(dist) for p in self.groundTruth]
        self.traj = newtraj
        self.trajValid = True
        
        
class TrajectoryMaker:
    def __init__(self, polygon):
        self.polygon = polygon
        self.trajectoryCollection = []
        
    def selectStartingFix(self):
        bnds = self.polygon.bounds
        while True: 
            long = random.uniform(bnds[0], bnds[2]) #easting
            lat = random.uniform(bnds[1],bnds[3]) #northing
            aPnt = Fix(lat,long)
            if aPnt.within(self.polygon):
               break
        return aPnt
    
    def makeGroundTruth(self, dt):
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
                nextFix = currFix.next_GPoint(currBrng, velocity * dt).toFix(t+dt)
                if nextFix.within(self.polygon):
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
            self.makeGroundTruth(2)
        
        for i in range(len(self.trajectoryCollection)):
            traj = self.trajectoryCollection[i]
            for fix in traj.groundTruth:
                data.append({"trajIndex": i,
                             "time": fix.time,
                             "lat": fix.geopyPoint[0],
                             "long": fix.geopyPoint[1],
                             "alt": fix.geopyPoint[2]})
        df = pd.DataFrame(data)
        df.to_csv('Data//' + filename +'.csv')

                
                    
   
def plotPoly(p):
    x,y = p.exterior.xy
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    ax.plot(x, y, color='#6699cc', alpha=0.7,
    linewidth=3, solid_capstyle='round', zorder=2)
    ax.set_title('Polygon')



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

xs, ys = pol1.exterior.xy
fig, axs = plt.subplots()
axs.fill(xs, ys, alpha=0.5, fc='r', ec='none')

traj = tm.makeGroundTruth(2)
x,y = traj.scatterXY()
axs.plot(x,y)

tm.makeDataSet('BarIlan', 20)
    



