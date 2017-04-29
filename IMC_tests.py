# -*- coding: utf-8 -*-
import sys, traceback
import pickle
from IMCObjects import Trajectory, Segment, GPoint
from TrajectoryMaker import TrajectoryCollectionCSVLoader, Building, TrajectoryMaker
from TraclusSegmenter import TrajectorySegmenter
from TraclusClusterer import SegmentsClusterer
from matplotlib import pyplot as plt
from shapely.geometry import MultiPoint
import time


def testTrajectoryMaking():
    print("testTrajectoryMaking: Run")
    fig, axs = plt.subplots()  
    bd = Building.buildingBuilder(1)
    tm = TrajectoryMaker(bd)
    t1 = tm.makeFixList(2)
    bd.plot(axs)
    t1.plot(axs)
    tm.makeDataSet('KfarSaba', 20)
    print("testTrajectoryMaking: End")

def testTrajectorySegmentation():
    loader = TrajectoryCollectionCSVLoader()
    loader.loadTrajectoryCollectionFromCSV('KfarSaba')
    
    trajectory = loader.trajectoryCollection[3]

    fixlist = TrajectorySegmenter.approximateTrajectoryPartitioning(trajectory)
    segmentedTrajectory = Trajectory(fixlist)

    
    fig, axs = plt.subplots(2,1)
    segmentedTrajectory.plot(axs[1], color = '#ea2ce1') #, linestyle = '--')
    trajectory.plot(axs[0])

def testPlotClusters():
    clusters = {}
    c0 = [Segment(GPoint(0,0), GPoint(0.0045, 0)),
          Segment(GPoint(0,0.001), GPoint(0.0045, 0.001)),
          Segment(GPoint(0,0), GPoint(0.001, 0.001))]
    
    clusters[0] = c0
    
    c1 = [Segment(GPoint(0.005,0), GPoint(0.005, 0.001)),
          Segment(GPoint(0.006,0), GPoint(0.006, 0.001)),
          Segment(GPoint(0.0055, 0.0005), GPoint(0.0065, 0.0005))]
    
    clusters[1] = c1
    fig, axs = plt.subplots()
    SegmentsClusterer.plotClusters(axs, list(clusters.values()))
    
    

def testTrajectoryClustering(withPickle = False):
    start = time.time()
    if not withPickle:
        print("{}: This is a long test, please be patient...".format(time.time() - start))
        loader = TrajectoryCollectionCSVLoader()
        print("{}: Loading Trajectory collection...".format(time.time() - start))
        trajectoryDict = loader.loadTrajectoryCollectionFromCSV('KfarSaba')
        print("{}: Trajectory collection loaded.".format(time.time() - start))
    
        print("{}: Starting a TrajectorySegmenter...".format(time.time() - start))    
        segmenter = TrajectorySegmenter(trajectoryDict)
        print("{}: TrajectorySegmenter started.".format(time.time() - start))
        print("{}: Starting to segment the trajectory collection...".format(time.time() - start))    
        segmentsList = segmenter.segmentsOfTrajectoryCollection(-1) #list of segment objects (GLine)
        print("{}: Segmentation of trajectory collection done. {} segments extracted.".format(time.time() - start, len(segmentsList)))
        
        print("{}: Starting a SegmentsClusterer. This will build a graph, might be very slow...".format(time.time() - start))
        clusterer = SegmentsClusterer(segmentsList, 12, 3)
        print("{}: SegmentsClusterer starting initActions...".format(time.time() - start))
        clusterer.initActions()
        with open("SegmentClusterer_15_1.p", "wb") as pickleFile:
            pickle.dump(clusterer, pickleFile)
        print("Saved clusterer as a pickle file, I did a lot of work to get it!")
    else:
        with open("SegmentClusterer_15_1.p", "rb") as pickleFile:
            clusterer = pickle.load(pickleFile)
    
    
    print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(clusterer.directReachablityGraph.nodes()),
          len(clusterer.directReachablityGraph.edges())))
    clusters = clusterer.LineSegmentClustering()
    print("{}: Clustering process ended.".format(time.time() - start))
    
    print("{}: Plotting the clusters...".format(time.time() - start))
#    fig, axs = plt.subplots(2,1)
    fig, axs = plt.subplots()
#    axs[0].set_title("Building")
#    axs[1].set_title("Reconstruction")
#    bd = Building.buildingBuilder(1)
#    bd.plot(axs[0])
#    clusterer.plotClusters(axs[1], list(clusters.values()))
    clusterer.plotClusters(axs, list(clusters.values()))
    print("{}: Plot ended.".format(time.time() - start))
    print("{}: Test ended.".format(time.time() - start))
    
def testTrajectoryClusteringWithPickle():
    start = time.time()
    with open("SegmentClusterer.p", "rb") as pickleFile:
        clusterer = pickle.load(pickleFile)
    
    print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(clusterer.directReachablityGraph.nodes()),
          len(clusterer.directReachablityGraph.edges())))
    clusters = clusterer.LineSegmentClustering()
    print("{}: Clustering process ended. Number of clusters is {}".format(time.time() - start, len(clusters)))
    
    print("{}: Plotting the clusters...".format(time.time() - start))
    fig, axs = plt.subplots()
    SegmentsClusterer.plotClusters(axs, clusters)
    print("{}: Plot ended.".format(time.time() - start))
    print("{}: Test ended.".format(time.time() - start))
    
try:
#    testTrajectoryMaking()    
##    testTrajectorySegmentation()  
    testTrajectoryClustering(withPickle=True)
##    testTrajectoryClusteringWithPickle()  
##    testPlotClusters()
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