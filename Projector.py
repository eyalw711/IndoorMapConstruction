# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 01:01:02 2017

@author: Eyal
"""
from matplotlib import pyplot as plt
import math
from IMCObjects import GPoint, GLine, LLine

class EquirectangularProjector:
    '''
    simple class for projections of small areas on earth to x,y
    reference: http://stackoverflow.com/a/16271669
    '''
    radius = 6371e3
    def __init__(self, segments_list_of_trajectory_collection):
        minLat = min(seg.minLatOfLineSeg() for seg in segments_list_of_trajectory_collection)
        maxLat = max(seg.maxLatOfLineSeg() for seg in segments_list_of_trajectory_collection)
        
        minLong = min(seg.minLongOfLineSeg() for seg in segments_list_of_trajectory_collection)
        maxLong = max(seg.maxLongOfLineSeg() for seg in segments_list_of_trajectory_collection)

        self.meanLat = (minLat + maxLat) / 2.0
        self.meanLong = (minLong + maxLong) / 2.0
        self.cosMeanLat = math.cos(self.meanLat)
        
        ox = EquirectangularProjector.radius * minLong * self.cosMeanLat
        oy = EquirectangularProjector.radius * minLat
        self.originXY = (ox, oy)
        
    def latLongToXY(self, lat, long):
        
        '''
        x = R * long * cos(meanLat)
        y = R * lat
        '''
        x = EquirectangularProjector.radius * long * self.cosMeanLat - self.originXY[0]
        y = EquirectangularProjector.radius * lat - self.originXY[1]
        return (x,y)
        
    def XYToLatLong(self, x, y):
        return ((y+self.latLongToXY[1])/ EquirectangularProjector.radius ,(x + self.latLongToXY[0])/(EquirectangularProjector.radius * self.cosMeanLat))        

    


class RotationalProjector:
    
    def __init__(self, phi = 0):
        self.phi = math.radians(phi)
    
    def rotate(self, phi):
        self.phi = math.radians(phi)
    
    def xyToRotatedXY(self, x, y):
        xtag = math.cos(self.phi) * x + math.sin(self.phi) * y
        ytag = -math.sin(self.phi) * x + math.cos(self.phi) * y
        return (xtag, ytag)
    
    def rotatedXYToXY(self, xtag, ytag):
        x = math.cos(-self.phi) * xtag + math.sin(-self.phi) * ytag
        y = -math.sin(-self.phi) * xtag + math.cos(-self.phi) * ytag
        return (x,y)


def testRotationalProjector():
    rproj = RotationalProjector()
    
    fig, axs = plt.subplots()
    ll = LLine((1,1), (2,2))
    axs.plot([ll.end1[0], ll.end2[0]],
             [ll.end1[1], ll.end2[1]], color = 'r')
    
    rproj.rotate(45)
    
    ll2 = LLine(rproj.xyToRotatedXY(ll[0][0], ll[0][1]),
                rproj.xyToRotatedXY(ll[1][0], ll[1][1]))
    
    axs.plot([ll2.end1[0], ll2.end2[0]],
             [ll2.end1[1], ll2.end2[1]], color = 'b')
    
        

    
#    def latLongToRotatedXY(self, lat, long):
#        raise Exception("Dont rotate in Equirectangular projector!")
#        x,y = self.latLongToXY(lat, long)
#        return self.xyToRotatedXY(x,y)
#        xtag = math.cos(self.phi) * x + math.sin(self.phi) * y
#        ytag = -math.sin(self.phi) * x + math.cos(self.phi) * y
#        return (xtag, ytag)
    
    
        

    
#    def rotatedXYToLatLong(self, xtag, ytag):
#        raise Exception("Dont rotate in Equirectangular projector!")
#        x,y = self.rotatedXYToXY(xtag, ytag)
#        return self.XYToLatLong(x,y)
    
#def testProjector():
#    proj = EquirectangularProjector([GLine(GPoint(34,32.005), GPoint(34.005,32.01))])
#    print("lat long" , 34, 32.005)
#    x, y = proj.latLongToXY(34,32.005)
#    print("x,y", x,y)
#    proj.rotate(45)
#    xr, yr = proj.xyToRotatedXY(x,y)
#    print("xr,yr", xr,yr)
#    xo, yo = proj.rotatedXYToXY(xr, yr)
#    print("xo,yo", xo,yo)
#    lato, longo = proj.XYToLatLong(xo, yo)
#    print("lato,longo", lato,longo)
        