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
    def __init__(self, fileName, respondersMacList, measurmentInCycle=43, minValidResponderInFix=5, floor=1):
        self.FileName = fileName
        self.RespondersMacList = respondersMacList
        # Collection of type LightRangeInfo
        self.RangeCollection = []
        self.NumOfMeasuInCycle = measurmentInCycle
        self.MinValidRespondersInFix = minValidResponderInFix
        self.Floor = floor

    def LoadRanges(self, detaliedLog=False):
        data = pd.read_csv('data//Floor' + str(self.Floor) + '//' + self.FileName + '.csv')
        internal = []
        filterCnt = 0
        for index, row in data.iterrows():
            if index != 0 and index % self.NumOfMeasuInCycle == 0:
                if detaliedLog:
                    print("{} valid ranges was found".format(len(internal)))
                if len(internal) >= self.MinValidRespondersInFix:
                    self.RangeCollection.append(internal)
                else:
                    filterCnt += 1
                    if detaliedLog:
                        print("Dropping group due to not enough valid measurments")

                internal = []
            status = row['status']
            r = row['range']
            if (status == "SUCCESS") and (r >= 0) and row['mac'] in self.RespondersMacList:
                internal.append(LightRangeInfo(row['mac'], row['time'], row['range']))

        if len(internal) >= self.MinValidRespondersInFix:
            if detaliedLog:
                print("{} valid ranges was found".format(len(internal)))
            self.RangeCollection.append(internal)
        else:
            filterCnt += 1
            if detaliedLog:
                print("Dropping group due to not enough valid measurments")

        print("{} fixes was found. {} was filtered due to not enough valid measurments".format(len(self.RangeCollection),filterCnt))

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
    def __init__(self, responderSet=None):
        self.ResponderDic = None
        self.InitLocationFinder(responderSet)

    def InitLocationFinder(self, responderSet):
        if responderSet is not None:
            self.ResponderDic = LocationFinder.GetResponderDic(responderSet)

    @staticmethod
    def GetResponderDic(responderSet):
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


class LocationEngine:
    def __init__(self, configuration):
        assert isinstance(configuration, LocationEngineConfig)
        self.EngineConfig = configuration
        self.ResponderLoader = ResponderLoader(configuration.ReponderFileName, configuration.Floor)
        self.FullFixCollection = []
        self.LocationFinder = LocationFinder()

        ''' Location of fixes'''
        self.xLoc = []
        self.yLoc = []
        self.LatLoc = []
        self.LongLoc = []

        self.TrajCsvData = []

    def InitEngine(self):
        print("Loading responders and initializeing LocationFinder")
        self.ResponderLoader.InitResponders()
        self.LocationFinder.InitLocationFinder(self.ResponderLoader.RespondersCollection)

    def PlotResponders(self, axs, color='b'):
        X = []
        Y = []
        for res in self.ResponderLoader.RespondersCollection:
            X.append(res.X)
            Y.append(res.Y)
        axs.plot(X, Y, "ro", color=color)

    def LoadRanges(self, isDump=False):
        i = 0
        for setName in self.EngineConfig.LocationFileNamesSet:
            print("Loading {}.csv range file".format(setName))
            rangeLoader = self.LoadRangeFile(setName)
            assert isinstance(rangeLoader, RangeLoader)
            if isDump:
                print("Dumping {} range file".format(setName))
                rangeLoader.DumpRanges()

            self.FullFixCollection.append(rangeLoader.RangeCollection)
            i += 1

    def LoadRangeFile(self, fileName):
        mac_respondersList = [r.Mac for r in self.ResponderLoader.RespondersCollection]
        rangeLoader = RangeLoader(fileName, mac_respondersList, self.EngineConfig.NumOfMeasurmentsInCycle, self.EngineConfig.MinRangesInFix, self.EngineConfig.Floor)
        rangeLoader.LoadRanges()
        return rangeLoader

    def ExtractFixLocation(self):
        projector = self.ResponderLoader.getProjector()
        assert isinstance(projector, EquirectangularProjector)

        print("Starting Extract Location")
        i = 0
        trajIndex = 0
        for fileCol in self.FullFixCollection:
            LongTemp = []
            LatTemp = []
            for colection in fileCol:
                d1 = datetime.datetime.now()
                tempX, tempY = self.LocationFinder.GetLocation(colection)
                tempLat, tempLong = projector.XYToLatLong(tempX, tempY)
                d2 = datetime.datetime.now()
                # print("Time {} for fix index:{}".format(d2 - d1, i))
                self.xLoc.append(tempX)
                self.yLoc.append(tempY)
                LatTemp.append(tempLat)
                LongTemp.append(tempLong)
                i += 1

            self.LatLoc += LatTemp
            self.LongLoc += LongTemp
            self.AddTrajectoryToTrajFile(LatTemp, LongTemp, trajIndex)
            trajIndex += 1
        assert len(self.xLoc) == len(self.yLoc)

    def AddTrajectoryToTrajFile(self, LatTemp, LongTemp, trajIndex):
        time = 0
        for lat, long in zip(LatTemp, LongTemp):
            self.TrajCsvData.append({"alt": 0, "lat": lat, "long": long, "time": time, "trajIndex": trajIndex})
            time += 1

    def DumpTrajCsv(self, fileName="TauTraj"):
        df = pd.DataFrame(self.TrajCsvData)
        df.to_csv('data//' + fileName + 'Dump.csv')

    def PlotFixes(self, axs, color="r"):
        axs.plot(self.xLoc, self.yLoc, 'ro', color=color)

    def PlotLocationToGoogleMaps(self, fileName = "Range_GMAP.html"):
        googleMapsConfig = GoogleMapsConfig(self.LatLoc,self.LongLoc,fileName)
        GoogleMapFactory.GenerateMap(googleMapsConfig)


class LocationEngineConfig:
    def __init__(self, responderFileName, floor, locationFileNames, numberOfRespondersInCycle, minRangesInFix=3):
        self.ReponderFileName = responderFileName
        self.LocationFileNamesSet = locationFileNames
        self.NumOfMeasurmentsInCycle = numberOfRespondersInCycle
        self.Floor = floor
        self.MinRangesInFix = minRangesInFix


def Run(fileName, axs, LAT, LONG):
    RangeFileName = fileName
    loader = ResponderLoader("Respodenrs", 2)
    loader.InitResponders()
    #   for rep in loader.RespondersCollection:
    # print(rep)

    X = [res.X for res in loader.RespondersCollection]
    Y = [res.Y for res in loader.RespondersCollection]

    axs.plot(X, Y, "ro", color="b")
    # plt.show()
    print("Starting Loading ranges")

    mac_respondersList = [r.Mac for r in loader.RespondersCollection]
    rangeLoader = RangeLoader(RangeFileName, mac_respondersList)
    rangeLoader.LoadRanges()

    rangeCollection = rangeLoader.RangeCollection

    print("Dumping ranges")
    rangeLoader.DumpRanges()

    locationFinder = LocationFinder(loader.RespondersCollection)

    print("Starting Fix Calculation")
    xLoc = []
    yLoc = []
    i = 0
    for colection in rangeCollection:
        d1 = datetime.datetime.now()
        tempX, tempY = locationFinder.GetLocation(colection)
        d2 = datetime.datetime.now()
        print("Time {} for fix index:{}".format(d2 - d1, i))
        xLoc.append(tempX)
        yLoc.append(tempY)
        i += 1

    axs.plot(xLoc, yLoc, 'ro', color="r")
    zipped = zip(xLoc, yLoc)
    assert len(xLoc) == len(yLoc)

    xPos2 = -1000
    yPos2 = -1000

    xPos1, yPos1 = zipped.__next__()
    for i in range(int(len(xLoc)-1)):
        if xPos2 != -1000:
       #     axs.plot([xPos1, xPos2], [yPos1, yPos2])
            xPos1 = xPos2
            yPos1 = yPos2
        xPos2, yPos2 = zipped.__next__()

    print("Starting generate Google maps plot")
    projector = loader.getProjector()
    assert isinstance(projector, EquirectangularProjector)

    for a, b in zip(xLoc, yLoc):
        la, lo = projector.XYToLatLong(a, b)
        LAT.append(la)
        LONG.append(lo)

    '''
    print("starting Cycle drawing")
    cs = []
    for rangeList in rangeCollection:
        for r in rangeList:
            curResp = locationFinder.ResponderDic[r.Mac]
            axs.plot(curResp[0], curResp[1], 'g^')
            cs.append(CircleUtility.circleFactory(curResp[0], curResp[1], r.Range))
    CircleUtility.DrawCircles(axs, cs)
'''
    axs.plot(60, 60)
    axs.plot(-20, -20)


def RunEngine(configuration, axs):
    assert isinstance(configuration,LocationEngineConfig)
    axs.set_title("Min Range in Fix {}".format(configuration.MinRangesInFix))

    Engine = LocationEngine(configuration)
    Engine.InitEngine()
    Engine.PlotResponders(axs)
    Engine.LoadRanges()
    Engine.ExtractFixLocation()
    Engine.PlotFixes(axs)
    Engine.PlotLocationToGoogleMaps()
    Engine.DumpTrajCsv()

    return


class FixFactory:
    @staticmethod
    def Generate():
        main()


def main():
    fig, axs = plt.subplots(1, 2)

    ''', "RangesLobby", "RangesClasses"'''
#    Engine = LocationEngine(LocationEngineConfig("Respodenrs", 2, ["RangesLabs", "RangesLobby", "RangesClasses"], 43))
    # RunEngine(LocationEngineConfig("Respodenrs", 2, ["Labs1", "Labs2"], 16, 3))

    RangeFileList = ["RangeFloor1", "RangeFloor1_1","RangeFloor1_2","RangeFloor1_3", "RangeFloor1_4", "RangeFloor1_5", "RangeFloor1_6"]
    # RangeFileList = ["RangeFloor1_4", "RangeFloor1_5", "RangeFloor1_6"]

    #RunEngine(LocationEngineConfig("RespodenrsFloor1", 1, RangeFileList, 16, 3), axs[0])
    RunEngine(LocationEngineConfig("RespodenrsFloor1", 1, RangeFileList, 16, 3), axs[1])
    axs[0].grid(True)
    axs[1].grid(True)
    plt.show()

if __name__ == "__main__":
    main()
