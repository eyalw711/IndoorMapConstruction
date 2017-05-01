# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 15:33:07 2017

@author: Eyal
"""
from IMCObjects import Trajectory, GLine, Segment
from math import log2

class TrajectorySegmenter:
    def __init__(self, trajectoryCollection):
        ''' a dictionary with <trajectoryIndex, TrajectoryClassObject> '''
        self.trajectoryCollection = trajectoryCollection
    
    def segmentsOfTrajectoryCollection(self, max_num_of_trajectories):
        ''' returns a list of Segment class objects '''
        segments = []
        i = 0
        for trajIndex, trajObject in self.trajectoryCollection.items():
#            print("segmentsOfTrajectoryCollection: {}/{}".format(i, len(self.trajectoryCollection.items())))
            if i == max_num_of_trajectories:  # only take part of DB
                return segments
            i+=1
            trajCPsFixList = TrajectorySegmenter.approximateTrajectoryPartitioning(trajObject)
#            print("partitioned {}/{} and now building segments".format(i, len(self.trajectoryCollection.items())))
            segmentsOfTraj = TrajectorySegmenter.getSegmentsFromCPs(trajIndex, trajCPsFixList)
            segments += segmentsOfTraj
        return segments
    
    def getSegmentsFromCPs(trajIndex, trajCPsFixList):
        ''' 
        param trajCPs needs to be a list of GPoints (could be fixes since inheritance)
        returns a list of Segment objects 
        '''
        segments = []
        for i in range(len(trajCPsFixList) - 1):
            cpstart = trajCPsFixList[i]
            cpend = trajCPsFixList[i+1]
            segment = Segment(cpstart, cpend, trajIndex)
            segments.append(segment)
        return segments
            
    
    def approximateTrajectoryPartitioning(trajectory):
        ''' 
        Returns Fix list of Characteristic Points 
        of the segmented trajectory
        '''
        CPs = [] #Characteristic Points
        CPs += [trajectory[0]] #starting point
        trajLen = len(trajectory.FixList)
        startIndex = 0
        length = 1
#        print("approximateTrajectoryPartitioning: length {}".format(trajLen - 1))
        while startIndex + length <= trajLen - 1:
            currIndex = startIndex + length
#            print("approximateTrajectoryPartitioning: going into MDL_par - {}".format(currIndex))
            costPar = TrajectorySegmenter.MDL_par(trajectory, startIndex, currIndex)
#            print("approximateTrajectoryPartitioning: going into MDL_noPar - {}".format(currIndex))
            costNoPar = TrajectorySegmenter.MDL_noPar(trajectory, startIndex, currIndex)
#            print("currIndex: {}, costPar: {}, costNoPar: {}".format(currIndex, costPar, costNoPar))
            #check if partitioning at the current point makes
            #the MDL cost larger than not partitioning
            if costPar > costNoPar:
                #partition at previous point
#                print("added to CPs - index:{}, fix:{}".format(currIndex - 1, trajectory[currIndex - 1]))
                CPs += [trajectory[currIndex - 1]]
                startIndex = currIndex - 1
                length = 1
            else:
                length += 1
        CPs += [trajectory[trajLen - 1]]
        return CPs
    
    def MDL_par(trajectory, startIndex, currIndex):
        '''
        denotes the MDL cost L(H)+L(D|H) of trajectory between PstartIndex and PcurrIndex
        (PstartIndex < PcurrIndex) when assuming that PstartIndex and PcurrIndex are only
        characteristic points. (one line end to end)
        '''
        
        endToEndLine = GLine(trajectory[startIndex], trajectory[currIndex])
        
        LH = log2(1 + endToEndLine.length())
        
        lineSegments = Trajectory.toLineSegmentList(trajectory[startIndex : currIndex + 1])
        
        #TODO: changed here to new distance metric
        try:
            LDH = log2(1 + sum(endToEndLine.myDistance(lineSeg) for lineSeg in lineSegments))
#        try:
#            LDH = log2(1 + sum(GLine.perpendicularDist(endToEndLine, lineSeg) for lineSeg in lineSegments)) +\
#                log2(1 + sum(GLine.angularDist(endToEndLine, lineSeg) for lineSeg in lineSegments))
        except ValueError:
            print("startInx= {}, currInx= {}".format(startIndex, currIndex))
            print("myDistances", [endToEndLine.myDistance(lineSeg) for lineSeg in lineSegments])
#            print("Perpendicular distances", [GLine.perpendicularDist(endToEndLine, lineSeg) for lineSeg in lineSegments])
#            print("Angular distances", [GLine.angularDist(endToEndLine, lineSeg) for lineSeg in lineSegments])
            raise SystemExit(0)
        return LH + LDH
    
    def MDL_noPar(trajectory, startIndex, currIndex):
        '''
        denotes the MDL cost L(H)+L(D|H) of trajectory between PstartIndex and PcurrIndex
        (PstartIndex < PcurrIndex) when assuming that there is no characteristic
        point between pi and pj , i.e., when preserving the original trajectory.
        (trajectory as is)
        '''
        
        #TODO: article says add small constant to cost_NoPar to get better clustering - make experiements
        lineSegments = Trajectory.toLineSegmentList(trajectory[startIndex : currIndex + 1])
        LH = log2(1 + sum(lineSeg.length() for lineSeg in lineSegments))
#        print("MDL_noPar - {}".format(LH))
        # LDH = 0 as seen in document
        return LH
        


