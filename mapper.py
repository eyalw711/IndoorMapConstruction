# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 18:19:19 2017

@author: Eyal
"""
import time
from TrajectoryMaker import TrajectoryCollectionCSVLoader
from TraclusSegmenter import TrajectorySegmenter
from TraclusClusterer import SegmentsClusterer
import pickle
from matplotlib import pyplot as plt

class IndoorMapper:
    
    def __init__(self, csvName, eps, MinLns):
        self.csvName = csvName
        self.eps = eps
        self.MinLns = MinLns
        
    def run(self, withPickle):
        start = time.time()
        pickle_success = True
        if withPickle:
            try:
                with open("SegmentClusterer_eps_{}_minlns_{}.p".format(self.eps, self.MinLns), "rb") as pickleFile:
                    clusterer = pickle.load(pickleFile)
            except Exception:
                pickle_success = False  
        if (not withPickle) or (withPickle and not pickle_success):
            print("{}: This is a long operation, please be patient...".format(time.time() - start))
            loader = TrajectoryCollectionCSVLoader()
            print("{}: Loading Trajectory collection...".format(time.time() - start))
            trajectoryDict = loader.loadTrajectoryCollectionFromCSV(self.csvName)
            print("{}: Trajectory collection loaded.".format(time.time() - start))
        
            print("{}: Starting a TrajectorySegmenter...".format(time.time() - start))    
            segmenter = TrajectorySegmenter(trajectoryDict)
            print("{}: TrajectorySegmenter started.".format(time.time() - start))
            print("{}: Starting to segment the trajectory collection...".format(time.time() - start))    
            segmentsList = segmenter.segmentsOfTrajectoryCollection(-1) #list of segment objects (GLine)
            print("{}: Segmentation of trajectory collection done. {} segments extracted.".format(time.time() - start, len(segmentsList)))
            
            print("{}: Starting a SegmentsClusterer. This will build a graph, might be very slow...".format(time.time() - start))
            
            #############################################
            #   From Here all operations are Local      #
            #############################################
            
            clusterer = SegmentsClusterer(segmentsList, self.eps, self.MinLns)
            
            print("{}: SegmentsClusterer starting initActions...".format(time.time() - start))
            clusterer.initActions()
            with open("SegmentClusterer_eps_{}_minlns_{}.p".format(self.eps, self.MinLns), "wb") as pickleFile:
                pickle.dump(clusterer, pickleFile)
            print("Saved clusterer as a pickle file, I did a lot of work to get it!")

        print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(clusterer.directReachablityGraph.nodes()),
              len(clusterer.directReachablityGraph.edges())))
        clusters = clusterer.LineSegmentClustering()
        print("{}: Clustering process ended. Got {} clusters.".format(time.time() - start, len(clusters.values())))
        
        print("{}: Plotting the clusters...".format(time.time() - start))
        fig, axs = plt.subplots()
        
        # temp operation...
        clusterer.plotClusters(axs, list(clusters.values()))
        print("{}: Plot ended.".format(time.time() - start))
        print("{}: Test ended.".format(time.time() - start))
    
    
def testKfarSaba(eps, MinLns):
    csvName = 'KfarSaba'
    mapper = IndoorMapper(csvName, eps, MinLns)
    mapper.run(False)
    return

testKfarSaba(40, 3)
