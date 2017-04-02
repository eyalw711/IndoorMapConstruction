# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 15:33:07 2017

@author: Eyal
"""
from TrajectoryMaker import Trajectory, TrajectoryCollectionCSVLoader, GLine, GPoint
from math import log2

class TrajectorySegmenter:
    def __init__(self, trajectoryCollection):
        self.trajectoryCollection = trajectoryCollection
        
    def approximateTrajectoryPartitioning(trajectory):
        CPs = [] #Characteristic Points
        CPs += trajectory[0] #starting point
        trajLen = len(trajectory.FixList)
        startIndex = 0
        length = 1
        while startIndex + length <= trajLen - 1:
            currIndex = startIndex + length
            costPar = TrajectorySegmenter.MDL_par(trajectory, startIndex, currIndex)
            costNoPar = TrajectorySegmenter.MDL_noPar(trajectory, startIndex, currIndex)
            #check if partitioning at the current point makes
            #the MDL cost larger than not partitioning
            if costPar > costNoPar:
                #partition at previous point
                CPs += trajectory[currIndex - 1]
                startIndex = currIndex - 1
                length = 1
            else:
                length += 1
        CPs += trajectory[-1]
        return CPs
    
    def MDL_par(trajectory, startIndex, currIndex):
        '''
        denotes the MDL cost L(H)+L(D|H) of trajectory between PstartIndex and PcurrIndex
        (PstartIndex < PcurrIndex) when assuming that PstartIndex and PcurrIndex are only
        characteristic points. (one line end to end)
        '''
        
        LH = log2(trajectory[startIndex].distance(trajectory[currIndex]))
        
        endToEndLine = GLine(trajectory[startIndex], trajectory[currIndex])
        lineSegments = Trajectory.toLineSegmentList(trajectory[startIndex : currIndex + 1])
        
        LDH = sum(log2(GLine.perpendicularDist(endToEndLine, lineSeg)) for lineSeg in lineSegments) \
            + sum(log2(GLine.angularDist(endToEndLine, lineSeg)) for lineSeg in lineSegments)
        
        return LH + LDH
    
    def MDL_noPar(trajectory, startIndex, currIndex):
        '''
        denotes the MDL cost L(H)+L(D|H) of trajectory between PstartIndex and PcurrIndex
        (PstartIndex < PcurrIndex) when assuming that there is no characteristic
        point between pi and pj , i.e., when preserving the original trajectory.
        (trajectory as is)
        '''
        lineSegments = Trajectory.toLineSegmentList(trajectory[startIndex : currIndex + 1])
        LH = sum(log2(lineSeg.length()) for lineSeg in lineSegments)
        # LDH = 0 as seen in document
        return LH

class Segment:
    def __init__(self, trajectory):
        self.trajectory = trajectory

class Cluster:
    def __init__(self, segmentsList):
        self.segments = segmentsList


loader = TrajectoryCollectionCSVLoader()
loader.loadTrajectoryCollectionFromCSV('KfarSaba')

