# -*- coding: utf-8 -*-
import sys, traceback
import pickle
from IMCObjects import Trajectory, Segment, GPoint
from TrajectoryMaker import TrajectoryCollectionCSVLoader, Building, TrajectoryMaker
from TraclusSegmenter import TrajectorySegmenter
from TraclusClusterer import SegmentsClusterer
from matplotlib import pyplot as plt
from shapely.geometry import MultiPoint
from mapper import IndoorMapper
import time


def testMapper():
    eps = 1.2
    minLns = 5
    scenes = {
        1: (Building.buildBuilding1, 'KfarSaba'),
        2: (Building.buildBuilding2, 'PTK')
    }
    sceneNumber = 1  # 'KfarSaba'
    builderfunc, sceneName = scenes[sceneNumber]
    fig, axs = plt.subplots(1, 2)
    fig.suptitle("Mapping {} with eps = {} minLns = {}".format(sceneName, eps, minLns))
    axs[0].set_title('Simulated Building')
    bd = builderfunc()
    bd.plot(axs[0])
    axs[1].set_title('Mapping Result')
    IndoorMapper.testMapperOn(axs[1], sceneName, eps, minLns)


def testTrajectoryMaking(buildingNum):
    siteNames = {1: 'KfarSaba', 2: 'PTK'}
    print("testTrajectoryMaking: Run")
    fig, axs = plt.subplots()
    fig.suptitle("testTrajectoryMaking")
    bd = Building.buildingBuilder(buildingNum)
    bd.plotOnGoogleMap(buildingNum)
    tm = TrajectoryMaker(bd)
    t1 = tm.makeFixList(1)
    bd.plot(axs)
    t1.plot(axs)
    tm.makeDataSet(siteNames[buildingNum], 25, 1)
    print("testTrajectoryMaking: End")


def testTrajectorySegmentation():
    loader = TrajectoryCollectionCSVLoader()
    loader.loadTrajectoryCollectionFromCSV('KfarSaba')

    trajectory = loader.trajectoryCollection[3]

    fixlist = TrajectorySegmenter.approximateTrajectoryPartitioning(trajectory)
    segmentedTrajectory = Trajectory(fixlist)

    fig, axs = plt.subplots(2, 1)
    segmentedTrajectory.plot(axs[1], color='#ea2ce1')  # , linestyle = '--')
    trajectory.plot(axs[0])


def testPlotClusters():
    clusters = {}
    c0 = [Segment(GPoint(0, 0), GPoint(0.0045, 0)),
          Segment(GPoint(0, 0.001), GPoint(0.0045, 0.001)),
          Segment(GPoint(0, 0), GPoint(0.001, 0.001))]

    clusters[0] = c0

    c1 = [Segment(GPoint(0.005, 0), GPoint(0.005, 0.001)),
          Segment(GPoint(0.006, 0), GPoint(0.006, 0.001)),
          Segment(GPoint(0.0055, 0.0005), GPoint(0.0065, 0.0005))]

    clusters[1] = c1
    fig, axs = plt.subplots()
    SegmentsClusterer.plotClusters(axs, list(clusters.values()))


def testTrajectoryClustering(withPickle=False):
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
        segmentsList = segmenter.segmentsOfTrajectoryCollection(-1)  # list of segment objects (GLine)
        print("{}: Segmentation of trajectory collection done. {} segments extracted.".format(time.time() - start,
                                                                                              len(segmentsList)))

        print("{}: Starting a SegmentsClusterer. This will build a graph, might be very slow...".format(
            time.time() - start))
        clusterer = SegmentsClusterer(segmentsList, 12, 3)
        print("{}: SegmentsClusterer starting initActions...".format(time.time() - start))
        clusterer.initActions()
        with open("SegmentClusterer_15_1.p", "wb") as pickleFile:
            pickle.dump(clusterer, pickleFile)
        print("Saved clusterer as a pickle file, I did a lot of work to get it!")
    else:
        with open("SegmentClusterer_15_1.p", "rb") as pickleFile:
            clusterer = pickle.load(pickleFile)

    print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(
        clusterer.directReachabilityGraph.nodes()),
                                                                                      len(
                                                                                          clusterer.directReachabilityGraph.edges())))
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

    print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(
        clusterer.directReachabilityGraph.nodes()),
                                                                                      len(
                                                                                          clusterer.directReachabilityGraph.edges())))
    clusters = clusterer.LineSegmentClustering()
    print("{}: Clustering process ended. Number of clusters is {}".format(time.time() - start, len(clusters)))

    print("{}: Plotting the clusters...".format(time.time() - start))
    fig, axs = plt.subplots()
    SegmentsClusterer.plotClusters(axs, clusters)
    print("{}: Plot ended.".format(time.time() - start))
    print("{}: Test ended.".format(time.time() - start))


try:
    testMapper()
    plt.show()
# testTrajectoryMaking(2)
##    testTrajectorySegmentation()  
#    testTrajectoryClustering(withPickle=True)
##    testTrajectoryClusteringWithPickle()  
##    testPlotClusters()
except:
    print("Exception:")
    print('-' * 60)
    traceback.print_exc(file=sys.stdout)
    print('-' * 60)
