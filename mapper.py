# -*- coding: utf-8 -*-
"""
Created on Sun Apr 30 18:19:19 2017

@author: Eyal
"""
import time
import copy
from TrajectoryMaker import TrajectoryCollectionCSVLoader
from TraclusSegmenter import TrajectorySegmenter
from TraclusClusterer import SegmentsClusterer
from ClusterProcessor import ClusterProcessor
import pickle
from matplotlib import pyplot as plt
import json
import geojson
from shapely.geometry import mapping
import networkx as nx

class IndoorMapper:
    
    def __init__(self, csvName, eps, MinLns):
        self.csvName = csvName
        self.eps = eps
        self.MinLns = MinLns
        
    def run(self, axs, withPickle):
        print("IndoorMapper: run")
        pickle_cluseterer_name = "pickles//SegmentClusterer_csvName_{}_minlns_{}".format(\
            self.csvName, self.MinLns).replace(".", "_") + ".p"
        start = time.time()
        pickle_success = True
        if withPickle:
            try:
                with open(pickle_cluseterer_name, "rb") as pickleFile:
                    print("{}: Loading a pickle".format(time.time() - start))
                    clusterer = pickle.load(pickleFile)
            except Exception:
                pickle_success = False  
        if (not withPickle) or (withPickle and not pickle_success): # all pickle data saved is unrelated to epsilon
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

            with open(pickle_cluseterer_name, "wb") as pickleFile:
                pickle.dump(clusterer, pickleFile)
            print("Saved clusterer as a pickle file, I did a lot of work to get it!")

        #############################################
        #   From Now - Epsilon Related Operations   #
        #############################################

        opt_eps = clusterer.set_epsilon(self.eps)
        # self.eps = opt_eps

        # self.MinLns = int(self.eps) + 2

        pickle_graph_name = "pickles//graph_csvName_{}_eps_{}_minlns_{}".format(\
            self.csvName, self.eps, self.MinLns).replace(".", "_") + ".p"

        pickle_success = True
        if withPickle:
            try:
                with open(pickle_graph_name, "rb") as pickleFile:
                    print("{}: Loading a pickle".format(time.time() - start))
                    graph = pickle.load(pickleFile)
            except Exception:
                pickle_success = False

        if (not withPickle) or (withPickle and not pickle_success):
            print("{}: Computing Direct Reachability Graph".format(time.time() - start))
            clusterer.computeDirectReachabilityGraph()
            with open(pickle_graph_name, "wb") as pickleFile:
                pickle.dump(clusterer.directReachabilityGraph, pickleFile)

        else:
            clusterer.directReachabilityGraph = graph

        print("{}: Starting clustering process... Graph has {} nodes and {} edges".format(time.time() - start, len(clusterer.directReachabilityGraph.nodes()),
              len(clusterer.directReachabilityGraph.edges())))
        clustersDict = clusterer.LineSegmentClustering()
        print("{}: Clustering process ended. Got {} clusters.".format(time.time() - start, len(clustersDict.values())))
        
        print("{}: Starting Cluster processing...".format(time.time() - start))
        processor = ClusterProcessor(self.MinLns)
        polygonsAndReprTrajs = processor.process(clustersDict.values())

        print("{}: Plotting the clusters...".format(time.time() - start))
        axs.set_title("Clustering eps = {} MinLns = {}".format(self.eps, self.MinLns))
        zip_for_map = copy.deepcopy(polygonsAndReprTrajs)
        processor.plot_clusters_and_reprs(axs, polygonsAndReprTrajs)
        map_poly, routing_graph = processor.make_map_from_zip_polygons_trajs(zip_for_map)

        geojson_map_string, geo_routing_graph = IndoorMapper.convert_local_to_geo(map_poly, routing_graph, clusterer.projector)

        floor_plan_file_name = "geojson//floor_plan_{}_{}.geojson".format(self.csvName, time.strftime("%d-%m-%Y_%H-%M-%S", time.gmtime()))
        with open(floor_plan_file_name, "w") as geojson_file:
            geojson_file.write(geojson_map_string)

        print("{}: Plot ended.".format(time.time() - start))
        print("{}: Test ended.".format(time.time() - start))

    def convert_local_to_geo(map_poly, routing_graph, projector):
        f_local_to_geo = lambda tup: projector.XYToLatLong(*tup)[::-1]
        f_local_edge_to_geo = lambda two_tups: (f_local_to_geo(two_tups[0]), f_local_to_geo(two_tups[1]))[::-1]

        local_map_mapping = mapping(map_poly)
        geo_map_mapping = IndoorMapper.geojson_map_coords(f_local_to_geo, local_map_mapping)
        geojson_map_string = geojson.dumps(geo_map_mapping)

        geo_graph = nx.Graph()

        local_nodes = list(routing_graph.nodes())
        geo_nodes = map(f_local_to_geo, local_nodes)

        local_edges = list(routing_graph.edges())
        geo_edges = map(f_local_edge_to_geo, local_edges)

        geo_graph.add_nodes_from(geo_nodes)
        geo_graph.add_edges_from(geo_edges)

        return geojson_map_string, geo_graph

    def testMapperOn(axs, csvName, eps, MinLns):
        mapper = IndoorMapper(csvName, eps, MinLns)
        mapper.run(axs, True)
        return

    def geojson_map_coords(func, obj):
        """
        :param func: Function to apply to tuples
        :type func: function
        :param obj: A geometry or feature to extract the coordinates from.
        :type obj: Point, LineString, MultiPoint, MultiLineString, Polygon,
        MultiPolygon
        :return: The result of applying the function to each tuple in the
        array.
        :rtype: list
        :raises ValueError: if the provided object is not a Geometry.
        """

        if obj['type'] == 'Point':
            coordinates = tuple(func(obj['coordinates']))
        elif obj['type'] in ['LineString', 'MultiPoint']:
            coordinates = [tuple(func(c)) for c in obj['coordinates']]
        elif obj['type'] in ['MultiLineString', 'Polygon']:
            coordinates = [[
                tuple(func(c)) for c in curve]
                for curve in obj['coordinates']]
        elif obj['type'] == 'MultiPolygon':
            coordinates = [[[
                tuple(func(c)) for c in curve]
                for curve in part]
                for part in obj['coordinates']]
        else:
            raise ValueError("Invalid geometry object %s" % repr(obj))
        return {'type': obj['type'], 'coordinates': coordinates}