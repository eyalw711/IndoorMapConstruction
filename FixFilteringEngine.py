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
import pandas as pd
import pandas as pd
import datetime
from IMCObjects import GPoint
from Projector import EquirectangularProjector
from Structs import Responder
from matplotlib import pyplot as plt
from TrajectoryMaker import TrajectoryCollectionCSVLoader
import pandas as pd
from IMCObjects import GPoint, Fix, Trajectory
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from shapely.geometry import Polygon, MultiPoint, LinearRing, LineString, MultiPolygon
from matplotlib import pyplot as plt
import random
import pandas as pd
from IMCObjects import GPoint, Fix, Trajectory
from Projector import EquirectangularProjector
import gmplot
from BaruchCode import trilateration
from Utility import CircleUtility
from GoogleMapsTool import GoogleMapsConfig, GoogleMapFactory
import xml.etree.ElementTree as ET
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


class CsvTool:
    @staticmethod
    def DumpToCsv(data, fileName):
        df = pd.DataFrame(data)
        df.to_csv(fileName)


class GoogleEarthTool:
    @staticmethod
    def ParseToCsvPointsFromKML(kmlFileName='Data\\Border\\Polygon.kml', csvFileName='Data\\Border\\DumpPolygon.csv'):
        PolygonEndPoints = GoogleEarthTool.KMLPointsParser(kmlFileName)
        CsvTool.DumpToCsv(GoogleEarthTool.GenerateDataToDump(PolygonEndPoints), csvFileName)
        return PolygonEndPoints

    ''' Return List of type GPoint'''

    @staticmethod
    def KMLPointsParser(kmlFileName):
        PointList = list()
        tree = ET.parse(kmlFileName)
        root = tree.getroot()
        points = root.findall(".//Placemark/Point/coordinates")
        i = 0
        for point in points:
            i += 1
            temp = str(point.text).split(',')
            PointList.append(GPoint(float(temp[1]), float(temp[0])))
            print("{} {}".format(i, temp))

        return PointList

    @staticmethod
    def GenerateDataToDump(endPoints):
        dataToDump = []
        for p in endPoints:
            assert isinstance(p, GPoint)
            dataToDump.append({"lat": p.getLat(), "long": p.getLong()})
        return dataToDump


class OuterWallsFiltering:
    def __init__(self, gPointsList, projector=None):
        self.GEndPointsList = gPointsList
        self.Projector = None
        self.Polygon = None

        assert isinstance(gPointsList, list)
        assert len(gPointsList) > 0

        if projector is None:
            self.InitProjector()
        else:
            assert isinstance(projector, EquirectangularProjector)
            self.Projector = projector

        self.GeneratePolygon()

    def InitProjector(self):
        print("Generate projector for OuterWallsFiltering")
        minLat = min(fix.getLat() for fix in self.GEndPointsList)
        maxLat = max(fix.getLat() for fix in self.GEndPointsList)
        minLong = min(fix.getLong() for fix in self.GEndPointsList)
        maxLong = max(fix.getLong() for fix in self.GEndPointsList)

        self.Projector = EquirectangularProjector(None, False, (minLat, minLong), (maxLat, maxLong))

    def GeneratePolygon(self):
        """ Close the polygon"""
        self.GEndPointsList.append(self.GEndPointsList[0])
        Coords = []
        for p in self.GEndPointsList:
            assert isinstance(p, GPoint)
            x, y = self.Projector.latLongToXY(p.getLat(), p.getLong())
            Coords.append((x, y))
        self.Polygon = Polygon(Coords)

    def IsInPolygon(self, point):
        assert isinstance(point, Point)
        return self.Polygon.contains(point)

    def Plot(self, axs):
        axs.plot(*self.Polygon.exterior.xy, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round',
                 zorder=2)


class FilteringFactoryGenerator:
    @staticmethod
    def GenerateFilter(projector):
        assert projector is not None
        assert isinstance(projector, EquirectangularProjector)
        EndPoints = GoogleEarthTool.ParseToCsvPointsFromKML()
        return OuterWallsFiltering(EndPoints, projector)

    @staticmethod
    def PureGenerateFilter():
        EndPoints = GoogleEarthTool.ParseToCsvPointsFromKML()
        return OuterWallsFiltering(EndPoints)


def main():
    filterEngine = FilteringFactoryGenerator.PureGenerateFilter()

    fig, axs = plt.subplots(1, 1)

    point = Point(0.5, 0.5)
    print(filterEngine.IsInPolygon(point))

    point = Point(40, 25)
    print(filterEngine.IsInPolygon(point))

    filterEngine.Plot(axs)
    plt.show()


if __name__ == "__main__":
    main()
