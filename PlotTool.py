from matplotlib import pyplot as plt
import time
import copy
from TraclusSegmenter import TrajectorySegmenter
from TraclusClusterer import SegmentsClusterer
from ClusterProcessor import ClusterProcessor
import pickle
from matplotlib import pyplot as plt
import json
import geojson
from shapely.geometry import mapping
import networkx as nx
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from shapely.geometry import Polygon, MultiPoint, LinearRing, LineString, MultiPolygon
from matplotlib import pyplot as plt
import random
import pandas as pd
from IMCObjects import GPoint, Fix, Trajectory
from Projector import EquirectangularProjector
import gmplot


class PlotUtility:
    @staticmethod
    def PlotAndBlok(X,Y):
        fig, axs = plt.subplots(1,1)
        axs.plot(X,Y, 'ro')
        plt.show()

    @staticmethod
    def PlotLatLongAndBlock(Lat, Long):
        minLat = min(lat for lat in Lat)
        maxLat = max(lat for lat in Lat)
        minLong = min(fix for fix in Long)
        maxLong = max(fix for fix in Long)

        projector = EquirectangularProjector(None, False, (minLat, minLong), (maxLat, maxLong))
        X = []
        Y = []
        for lat, long in zip(Lat,Long):
            x,y = projector.latLongToXY(lat, long)
            X.append(x)
            Y.append(y)
        PlotUtility.PlotAndBlok(X,Y)

