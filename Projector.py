# -*- coding: utf-8 -*-
"""
Created on Sat Apr 29 01:01:02 2017

@author: Eyal
"""
import math

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
        self.originXY = self.latLongToXY(minLat, minLong)
        self.phi = 0
        
    def rotate(self, phi):
        self.phi = phi
    
    def latLongToRotatedXY(self, lat, long):
        x,y = self.latLongToXY(lat, long)
        xtag = math.cos(self.phi) * x + math.sin(self.phi) * y
        ytag = -math.sin(self.phi) * x + math.cos(self.phi) * y
        return (xtag, ytag)
    
    def rotatedXYToXY(self, xtag, ytag):
        x = math.cos(-self.phi) * xtag + math.sin(-self.phi) * ytag
        y = -math.sin(-self.phi) * xtag + math.cos(-self.phi) * ytag
        return (x,y)
        
    def latLongToXY(self, lat, long):
        '''
        x = R * long * cos(meanLat)
        y = R * lat
        '''
        return (EquirectangularProjector.radius * long * self.cosMeanLat, EquirectangularProjector.radius * lat)
        
    def XYToLatLong(self, x, y):
        return (y/ EquirectangularProjector.radius ,x/(EquirectangularProjector.radius * self.cosMeanLat))
    
    def rotatedXYToLatLong(self, xtag, ytag):
        x,y = self.rotatedXYToXY(xtag, ytag)
        return self.XYToLatLong(x,y)