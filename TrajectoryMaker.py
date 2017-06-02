# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 14:10:18 2017

@author: Eyal
"""
from geopy.distance import vincenty, VincentyDistance
from geopy.distance import Point as gpt
from shapely.geometry import Polygon, MultiPoint, LinearRing, LineString, MultiPolygon
from matplotlib import pyplot as plt
import random
import pandas as pd
from IMCObjects import GPoint, Fix, Trajectory
from Projector import EquirectangularProjector
import gmplot


class Building:
    def __init__(self, outerWallsLinearRing, innerWallsMultiLineString=None, holesPolygonList=[]):
        """
        param outerWallsRing: a LinearRing class object for outer walls. include points for walls connection. <long, lat>!!!!
        param innerWallsList: a list of LineString class objects for inner walls of the building.
        """
        self.outerWallsLinearRing = outerWallsLinearRing
        self.innerWallsMultiLineString = innerWallsMultiLineString
        self.holesPolygonList = holesPolygonList

    def plotOnGoogleMap(self, buildingNum):
        longs, lats = zip(*self.outerWallsLinearRing.coords)
        gmap = gmplot.GoogleMapPlotter(lats[0], longs[0], 20)
        gmap.plot(lats, longs)
        gmap.draw("building_{}_GMAP.html".format(buildingNum))

    def getPolygon(self):
        outerWallsCoords = self.outerWallsLinearRing.coords
        polygon = Polygon(outerWallsCoords)
        if self.holesPolygonList != None:
            polygon = polygon.difference(MultiPolygon(self.holesPolygonList))
        return polygon

    def containsPoint(self, gpoint):
        polygon = self.getPolygon()
        return gpoint.within(polygon)

    def legalTravel(self, pt1, pt2):
        line = LineString([pt1.getCoords(), pt2.getCoords()])
        illegal = line.intersects(self.outerWallsLinearRing) or any(
            line.intersects(geo) for geo in [self.innerWallsMultiLineString] + self.holesPolygonList if geo != None)
        return not illegal

    def buildingBuilder(typenum):
        if typenum == 1:
            return Building.buildBuilding1()
        elif typenum == 2:
            return Building.buildBuilding2()
        else:
            print("No implementation for typenum = {}".format(typenum))
            return None

    def buildBuilding2():

        shiftpoint = (32.114104, 34.803086)
        projector = EquirectangularProjector(None, useSegments=False, minLat_minLong=shiftpoint)

        #        shiftPointInXY = projector.latLongToXY(*shiftpoint)
        #        print("ShiftPoint in XY: {}".format(shiftPointInXY))
        def addTupleElementWise(t1, t2):
            return (t1[0] + t2[0], t1[1] + t2[1])

        #
        #        outer_coords_l = [addTupleElementWise(xy, shiftPointInXY) for xy in [(0, 70), (70, 70), (70, 0), (0, 0), (0, 70)]]
        #        hole1_coords_l = [addTupleElementWise(xy, shiftPointInXY) for xy in [(10, 60), (30, 60), (30, 40), (10, 40), (10,60)]]
        #        hole2_coords_l = [(x + 30, y) for (x,y) in hole1_coords_l] #shift 30 meters on x
        #        hole3_coords_l = [(x, y - 30) for (x,y) in hole2_coords_l] #shift 30 meters down on y of hole 2
        #        hole4_coords_l = [(x, y - 30) for (x,y) in hole1_coords_l] #shift 30 meters down on y of hole 1

        outer_coords_l = [(0, 70), (70, 70), (70, 0), (0, 0), (0, 70)]
        hole1_coords_l = [(10, 60), (30, 60), (30, 40), (10, 40), (10, 60)]
        hole2_coords_l = [(x + 30, y) for (x, y) in hole1_coords_l]  # shift 30 meters on x
        hole3_coords_l = [(x, y - 30) for (x, y) in hole2_coords_l]  # shift 30 meters down on y of hole 2
        hole4_coords_l = [(x, y - 30) for (x, y) in hole1_coords_l]  # shift 30 meters down on y of hole 1

        outer_coords_g = [projector.XYToLatLong(*e)[::-1] for e in
                          outer_coords_l]  # reverse because long is "x" and lat is "y"
        hole1_coords_g = [projector.XYToLatLong(*e)[::-1] for e in hole1_coords_l]
        hole2_coords_g = [projector.XYToLatLong(*e)[::-1] for e in hole2_coords_l]
        hole3_coords_g = [projector.XYToLatLong(*e)[::-1] for e in hole3_coords_l]
        hole4_coords_g = [projector.XYToLatLong(*e)[::-1] for e in hole4_coords_l]

        holesPolygons_list = [Polygon(hole1_coords_g), Polygon(hole2_coords_g), Polygon(hole3_coords_g),
                              Polygon(hole4_coords_g)]
        #        holesMultiPoly = MultiPolygon(holesPolygonsList)
        outerWallsLinRing = LinearRing(outer_coords_g)
        building2 = Building(outerWallsLinRing, holesPolygonList=holesPolygons_list)
        #        fig, axs = plt.subplots()
        #        building2.plot(axs)
        return building2

    def buildBuilding1():
        outerWallsCrds = list()
        meetingInnerOuter = list()
        startPoint = GPoint(32.178074, 34.911721)
        currPoint = GPoint(32.178074, 34.911721)  # Start as north westren point
        outerWallsCrds.append(currPoint.getCoords())
        for i in range(5):
            nextPoint = currPoint.travel(90, 5)
            outerWallsCrds.append(nextPoint.getCoords())
            if not (i == 4):
                meetingInnerOuter.append(nextPoint)
            currPoint = nextPoint

        nextPoint = currPoint.travel(180, 10)
        outerWallsCrds.append(nextPoint.getCoords())
        currPoint = nextPoint

        nextPoint = currPoint.travel(270, 25)
        outerWallsCrds.append(nextPoint.getCoords())
        outerWallsCrds.append(startPoint.getCoords())

        outerWallsLinRing = LinearRing(outerWallsCrds)

        innerWallsCrds = list()
        for i in [0, 2]:
            innerWallsCrds.append(
                [meetingInnerOuter[i].getCoords(),
                 meetingInnerOuter[i].travel(180, 8).getCoords(),
                 meetingInnerOuter[i + 1].travel(180, 8).getCoords(),
                 meetingInnerOuter[i + 1].getCoords()])
        #        innerWallsMultiPolygon = MultiPolygon([Polygon(crds) for crds in innerWallsCrds])

        return Building(outerWallsLinRing, holesPolygonList=[Polygon(crds) for crds in
                                                             innerWallsCrds])  # holesMultiPolygon = innerWallsMultiPolygon)

    def plot(self, axs):
        polygon = self.getPolygon()  # Polygon(self.outerWallsLinearRing)
        xs, ys = polygon.exterior.xy
        axs.fill(xs, ys, alpha=0.5, fc='r', ec='none')

        if self.innerWallsMultiLineString != None:
            for line in self.innerWallsMultiLineString:
                x, y = line.xy
                axs.plot(x, y, color='#6699cc', alpha=0.7, linewidth=3, solid_capstyle='round', zorder=2)

        for phole in self.holesPolygonList:
            xs, ys = phole.exterior.xy
            axs.fill(xs, ys, alpha=1, fc='black', ec='none')


class TrajectoryMaker:
    def __init__(self, building):
        self.building = building
        self.trajectoryCollection = []

    def selectStartingFix(self):
        bnds = self.building.getPolygon().bounds
        while True:
            long = random.uniform(bnds[0], bnds[2])  # easting
            lat = random.uniform(bnds[1], bnds[3])  # northing
            aPnt = Fix(lat, long)
            if aPnt.within(self.building.getPolygon()):
                break
        return aPnt

    def makeFixList(self, dt):
        '''returns a trajectory of Fixes inside the TrajectoryMaker's Polygon'''
        # walk 5km/h = 1.389 m/s
        velocity = 1.389
        # walk for 2 minutes
        T = 120

        t = 0

        currFix = self.selectStartingFix()
        aTraj = Trajectory([currFix])
        currBrng = random.uniform(0, 360)
        while t < T:
            i = 0
            while True:

                bearingNoise = random.uniform(-75, 75)
                nextFix = currFix.travel(currBrng + bearingNoise, velocity * dt).toFix(t + dt)
                if self.building.legalTravel(currFix, nextFix):
                    # if nextFix.within(self.building):
                    #                    print("t is {}".format(t))
                    t += dt
                    aTraj.append(nextFix)
                    currFix = nextFix
                    break
                else:
                    currBrng += random.uniform(-45, 45)
                    i += 1
                    if i > 500:
                        fig, axs = plt.subplots()
                        axs.set_title("failed 500 times to make a fix list")
                        print("failed 500 times to make a fix list")
                        self.building.plot(axs)
                        axs.scatter([currFix.getLong()], [currFix.getLat()])
                        return None
        self.trajectoryCollection.append(aTraj)
        return aTraj

    def makeDataSet(self, filename, numberOfTrajectories, secPerFix):
        data = []
        for i in range(numberOfTrajectories):
            print("makeDataSet: making trajectory {}".format(i))
            self.makeFixList(secPerFix)

        for i in range(len(self.trajectoryCollection)):
            traj = self.trajectoryCollection[i]
            for fix in traj.FixList:
                data.append({"trajIndex": i,
                             "time": fix.time,
                             "lat": fix.geopyPoint[0],
                             "long": fix.geopyPoint[1],
                             "alt": fix.geopyPoint[2]})
        df = pd.DataFrame(data)
        df.to_csv('Data//' + filename + '.csv')


class TrajectoryCollectionCSVLoader:
    def __init__(self):
        '''each element in the dict is of a Trajectory class object'''
        self.trajectoryCollection = {}

    def loadTrajectoryCollectionFromCSV(self, filename):
        data = pd.read_csv('data//' + filename + '.csv')
        for index, row in data.iterrows():
            trajectoryIndex = row['trajIndex']
            if not trajectoryIndex in self.trajectoryCollection:
                self.trajectoryCollection[trajectoryIndex] = Trajectory(list())
            self.trajectoryCollection[trajectoryIndex].append(Fix(row['lat'], row['long'], row['alt'], row['time']))
        return self.trajectoryCollection


'''unused'''


def plotPoly(p, plotColor='#6699cc'):
    x, y = p.exterior.xy
    fig = plt.figure(1, figsize=(5, 5), dpi=90)
    ax = fig.add_subplot(111)
    ax.plot(x, y, color=plotColor, alpha=0.7,
            linewidth=3, solid_capstyle='round', zorder=2)
    ax.set_title('Polygon')


def test1():
    bldng1 = [(32.178028, 34.912227),
              (32.178022, 34.912530),
              (32.177977, 34.912527),
              (32.177972, 34.912886),
              (32.177927, 34.912889),
              (32.177929, 34.912444),
              (32.177976, 34.912443),
              (32.177985, 34.912226)]

    ks_p1 = (32.177916, 34.911915)
    ks_p2 = (32.177925, 34.911857)
    d12 = vincenty(ks_p1, ks_p2).meters
    ptk_p1 = (32.101759, 34.850336)
    ks_p3 = VincentyDistance(200 / 1000).destination(ks_p1, 90)

    # pol1 = Polygon(bldng1)
    pol1 = Polygon([(lon, lat) for (lat, lon) in bldng1])
    pol2 = MultiPoint(bldng1).convex_hull

    # plotPoly(pol1)
    # plot(pol2)

    tm = TrajectoryMaker(pol1)

    fig, axs = plt.subplots()
    xs, ys = pol1.exterior.xy
    axs.fill(xs, ys, alpha=0.5, fc='r', ec='none')

    traj = tm.makeFixList(2)
    x, y = traj.scatterXY()
    axs.plot(x, y)

    tm.makeDataSet('BarIlan', 20, 2)

# bd2 = Building.buildBuilding2()
# tm = TrajectoryMaker(bd2)
# tm.makeDataSet('PTK', 20, 1)


# traj = tm.makeFixList(1)

# fig, axs = plt.subplots()
# bd2.plot(axs)
# traj.plot(axs)
# pts = [tm.selectStartingFix() for i in range(10)]
# longs = [pt.getLong() for pt in pts]
# lats = [pt.getLat() for pt in pts]
# axs.scatter(longs, lats)
