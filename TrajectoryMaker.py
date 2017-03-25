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

def shapePt(geoPt):
    return spt(geoPt[0], geoPt[1])
    
class TrajectoryMaker:
    def __init__(self, polygon):
        self.poly = polygon
    
    def selectStartingPoint(self):
        bnds = self.poly.bounds
        while True: 
            x = random.uniform(bnds[0], bnds[2])
            y = random.uniform(bnds[1],bnds[3])
            pnt = spt(x,y)
            if pnt.within(self.poly):
               break
        return gpt(pnt.x, pnt.y)
    
    def makeTrajWithoutNoise(self, dt):
        '''returns a trajectory of geoPoints'''
        #walk 5km/h = 1.389 m/s
        velocity = 1.389
        #walk for 2 minutes
        T = 30
        
        t = 0
        
        geo_currPt = self.selectStartingPoint()
        shape_traj = [shapePt(geo_currPt)]
        currBrng = random.uniform(0,360)
        while t < T:
            while True:
                geo_nextPt = VincentyDistance((velocity * dt)/1000).destination(geo_currPt, currBrng)
                if shapePt(geo_nextPt).within(self.poly):
                    t += dt
                    shape_traj.append(shapePt(geo_nextPt))
                    geo_currPt = geo_nextPt
                    break
                else:
                    currBrng += random.uniform(-45,45)
        return shape_traj
                    
                

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

pol1 = Polygon(bldng1)
pol2 = MultiPoint(bldng1).convex_hull    

#plotPoly(pol1)
#plot(pol2)

tm = TrajectoryMaker(pol1)
pt1 = tm.selectStartingPoint()
print(pt1)

xs, ys = pol1.exterior.xy
fig, axs = plt.subplots()
axs.fill(xs, ys, alpha=0.5, fc='r', ec='none')
#axs.scatter(pt1.x, pt1.y)

shape_traj = tm.makeTrajWithoutNoise(2)
axs.plot([p.x for p in shape_traj], [p.y for p in shape_traj])
plt.show()

    



