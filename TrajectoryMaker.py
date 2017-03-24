# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:10:18 2017

@author: Eyal
""" 
from geopy.distance import vincenty, Point, VincentyDistance

ks_p1 = (32.177916, 34.911915)
ks_p2 = (32.177925, 34.911857)
ptk_p1 = (32.101759, 34.850336)
print(vincenty(ks_p1, ks_p2).meters)
print(vincenty(ks_p1, ptk_p1).meters)

ks_p3 = VincentyDistance(2).destination(ks_p1, 90)
print(ks_p3[0], ks_p3[1])


    