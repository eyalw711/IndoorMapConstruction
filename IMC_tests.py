# -*- coding: utf-8 -*-
import sys, traceback
from objects import Trajectory
from TrajectoryMaker import TrajectoryCollectionCSVLoader, Building, TrajectoryMaker
from TraclusSegmenter import TrajectorySegmenter
from matplotlib import pyplot as plt



def testTrajectoryMaking():
    fig, axs = plt.subplots()  
    bd = Building.buildingBuilder(1)
    tm = TrajectoryMaker(bd)
    t1 = tm.makeFixList(2)
    bd.plot(axs)
    t1.plot(axs)
    tm.makeDataSet('KfarSaba', 20)

def testTrajectorySegmentation():
    loader = TrajectoryCollectionCSVLoader()
    loader.loadTrajectoryCollectionFromCSV('KfarSaba')
    
    trajectory = loader.trajectoryCollection[3]

    fixlist = TrajectorySegmenter.approximateTrajectoryPartitioning(trajectory)
    segmentedTrajectory = Trajectory(fixlist)

    
    fig, axs = plt.subplots(2,1)
    segmentedTrajectory.plot(axs[1], color = '#ea2ce1') #, linestyle = '--')
    trajectory.plot(axs[0])

try:
#    testTrajectoryMaking()    
    testTrajectorySegmentation()    
except:
    print("Exception:")
    print('-'*60)
    traceback.print_exc(file=sys.stdout)
    print('-'*60)
    
#frm = nv.FrameE(name = 'WGS84', a=6371e3, f= 0)
#pta1 = frm.GeoPoint(32,34, degrees = True)
#pta2 = frm.GeoPoint(31.999, 34.001, degrees = True)
#d, b1, b2 = pta1.distance_and_azimuth(pta2, degrees = True)
#print("dist {}m, bearing1 {}, bearing2 {}".format(d,b1,b2))