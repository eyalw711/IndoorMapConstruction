# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 18:19:19 2017

@author: Eyal
"""
import time
from TrajectoryMaker import TrajectoryCollectionCSVLoader
from TraclusSegmenter import TrajectorySegmenter
from TraclusClusterer import SegmentsClusterer
from ClusterProcessor import ClusterProcessor
import pickle
from matplotlib import pyplot as plt

class IndoorMapper:
    
    def __init__(self, csvName, eps, MinLns):
        self.csvName = csvName
        self.eps = eps
        self.MinLns = MinLns
        
    def run(self, axs, withPickle):
        print("IndoorMapper: run")
        pickleName = "pickles//SegmentClusterer_csvName_{}_eps_{}_minlns_{}".format(\
            self.csvName, self.eps, self.MinLns).replace(".", "_") + ".p"
        start = time.time()
        pickle_success = True
        if withPickle:
            try:
                with open(pickleName, "rb") as pickleFile:
                    print("{}: Loading a pickle".format(time.time() - start))
                    clusterer = pickle.load(pickleFile)
            except Exception:
                pickle_success = False  
        if (not withPickle) or (withPickle and not pickle_success):
            print("{}: Running without pickle".format(time.time() - start))
            print("{}: This is a long operation, please be patient...".format(time.time() - start))
            loader = TrajectoryCollectionCSVLoader()
            print("{}: Loading Trajectory collection...".format(time.time() - start))
            trajectoryDict = loader.loadTrajectoryCollectionFromCSV(self.csvName)
            print("{}: Trajectory collection loaded.".format(time.time() - start))
        
            print("{}: Starting a TrajectorySegmenter...".format(time.time() - start))    
            segmenter = TrajectorySegmenter(trajectoryDict)
            print("{}: TrajectorySegmenter started.".format(time.time() - start))
            print("{}: Starting to segment the trajectory collection...".format(time.time() - start))    
            segmentsList = segmenter.segmentsOfTrajectoryCollection(25) #list of segment objects (GLine)
            print("{}: Segmentation of trajectory collection done. {} segments extracted.".format(time.time() - start, len(segmentsList)))
            
            print("{}: Starting a SegmentsClusterer. This will build a graph, might be very slow...".format(time.time() - start))
            
            #############################################
            #   From Here all operations are Local      #
            #############################################
            
            clusterer = SegmentsClusterer(segmentsList, self.eps, self.MinLns)
            
            print("{}: SegmentsClusterer starting initActions...".format(time.time() - start))
            clusterer.initActions()
            with open(pickleName, "wb") as pickleFile:
                pickle.dump(clusterer, pickleFile)
            print("Saved clusterer as a pickle file, I did a lot of work to get it!")

        print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(clusterer.directReachablityGraph.nodes()),
              len(clusterer.directReachablityGraph.edges())))
        clustersDict = clusterer.LineSegmentClustering()
        print("{}: Clustering process ended. Got {} clusters.".format(time.time() - start, len(clustersDict.values())))
        
#        clusterer.plotClusters(axs, clustersDict.values())
        
        print("{}: Starting Cluster processing...".format(time.time() - start))
        processor = ClusterProcessor(self.MinLns)
        polygonsAndReprTrajs = processor.process(clustersDict.values())

        print("{}: Plotting the clusters...".format(time.time() - start))
        axs.set_title("Clustering eps = {} MinLns = {}".format(self.eps, self.MinLns))
        processor.plotMap(axs, polygonsAndReprTrajs)

       # temp operation...
#        clusterer.plotClusters(axs, list(clustersDict.values()))
                
        print("{}: Plot ended.".format(time.time() - start))
        print("{}: Test ended.".format(time.time() - start))
    
    def testMapperOn(axs, csvName, eps, MinLns):
        mapper = IndoorMapper(csvName, eps, MinLns)
        mapper.run(axs, True)
        return


#test('KfarSaba', 60, 3)

#testKfarSaba(60, 4)
#testKfarSaba(60, 3)