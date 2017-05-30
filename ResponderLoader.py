import pandas as pd

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


class ResponderLoader:
    def __init__(self, fileName, floor):
        assert (fileName is not None)
        assert (floor is not None and isinstance(floor, int))
        self.FileName = fileName
        self.Floor = floor
        self.Projector = None
        self.RespondersCollection = []

    def LoadRespondersCoords(self):
        Coords = []
        data = pd.read_csv('data//' + self.FileName + '.csv')
        for index, row in data.iterrows():
            Coords.append(Responder(row['lat'], row['long'], 0, 0, row['mac'], 0))
        return Coords

    def InitResponders(self):
        respondersGCoords = self.LoadRespondersCoords()

        minLat = min(fix.getLat() for fix in respondersGCoords)
        maxLat = max(fix.getLat() for fix in respondersGCoords)
        minLong = min(fix.getLong() for fix in respondersGCoords)
        maxLong = max(fix.getLong() for fix in respondersGCoords)

        self.Projector = EquirectangularProjector(None, False, (minLat, minLong), (maxLat, maxLong))
        for fix in respondersGCoords:
            newX, newY = self.Projector.latLongToXY(fix.getLat(), fix.getLong())
            self.RespondersCollection.append(Responder(fix.getLat(), fix.getLong(), newX, newY, fix.Mac, self.Floor))

    def getProjector(self):
        return self.Projector


class RangeLoader:
    def __init__(self, fileName, respondersMacList, measurmentInCycle=43):
        self.FileName = fileName
        self.RespondersMacList = respondersMacList
        # Collection of type LightRangeInfo
        self.RangeCollection = []
        self.NumOfMeasuInCycle = measurmentInCycle

    def LoadRanges(self):
        data = pd.read_csv('data//' + self.FileName + '.csv')
        internal = []
        for index, row in data.iterrows():
            if index != 0 and index % self.NumOfMeasuInCycle == 0:
                print("{} valid ranges was found".format(len(internal)))
                if len(internal) >= 3:
                    self.RangeCollection.append(internal)
                else:
                    print("Dropping group due to not enough valid measurments")

                internal = []
            status = row['status']
            r = row['range']
            if (status == "SUCCESS") and (r >= 0) and row['mac'] in self.RespondersMacList:
                internal.append(LightRangeInfo(row['mac'], row['time'], row['range']))

        if len(internal) > 0:
            print("{} valid ranges was found".format(len(internal)))
            self.RangeCollection.append(internal)

        print("{} fixes was found".format(len(self.RangeCollection)))

    def DumpRanges(self):
        for l in self.RangeCollection:
            data = [{"mac": r.Mac, "time": r.Time, "range": r.Range} for r in l]
        df = pd.DataFrame(data)
        df.to_csv('data//' + self.FileName + 'Dump.csv')


class LightRangeInfo:
    def __init__(self, mac, time, distance):
        self.Mac = mac.lower()
        self.Time = time

        # convert from cm to m
        self.Range = distance / 100


class LocationFinder:
    def __init__(self, responderSet):
        self.ResponderDic = self.GetResponderDic(responderSet)

    def GetResponderDic(self, responderSet):
        dic = {}
        for resp in responderSet:
            assert (isinstance(resp, Responder))
            dic[resp.Mac] = [resp.X, resp.Y]
        return dic

    def GetLocation(self, locationSet):
        locationDic = {}
        for loc in locationSet:
            assert (isinstance(loc, LightRangeInfo))
            locationDic[loc.Mac] = loc.Range
        x, y = trilateration(self.ResponderDic, locationDic)

        return x, y


if __name__ == "__main__":

    RangeFileName = "RangesClasses"
    loader = ResponderLoader("Respodenrs", 2)
    loader.InitResponders()
    #   for rep in loader.RespondersCollection:
    # print(rep)

    fig, axs = plt.subplots(1, 1)
    X = [res.X for res in loader.RespondersCollection]
    Y = [res.Y for res in loader.RespondersCollection]

    axs.plot(X, Y, "ro", color="b")
    # plt.show()
    print("Starting Loading ranges")

    mac_respondersList = [r.Mac for r in loader.RespondersCollection]
    rangeLoader = RangeLoader(RangeFileName, mac_respondersList)
    rangeLoader.LoadRanges()

    rangeCollection = rangeLoader.RangeCollection

    rangeLoader.DumpRanges()

    locationFinder = LocationFinder(loader.RespondersCollection)

    xLoc = []
    yLoc = []
    for colection in rangeCollection:
        tempX, tempY = locationFinder.GetLocation(colection)
        xLoc.append(tempX)
        yLoc.append(tempY)

    axs.plot(xLoc, yLoc, 'ro', color="r")
    

    cs = []
    for rangeList in rangeCollection:
        for r in rangeList:
            curResp = locationFinder.ResponderDic[r.Mac]
            axs.plot(curResp[0], curResp[1], 'g^')
            cs.append(CircleUtility.circleFactory(curResp[0], curResp[1], r.Range))
    # CircleUtility.DrawCircles(axs, cs)

    axs.plot(60, 60)
    axs.plot(-20, -20)
    plt.grid(True)
