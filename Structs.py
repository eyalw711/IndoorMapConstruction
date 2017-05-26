from IMCObjects import GPoint
from Utility import CircleUtility
import math
from sympy.geometry import Circle, Point


class Responder(GPoint):
    def __init__(self, lat, long, x, y, mac, floor):
        GPoint.__init__(self, lat, long)
        self.X = x
        self.Y = y
        self.Mac = mac
        self.Floor = floor

    def __repr__(self):
        return "Responder (lat: {}, long: {}, x: {}, y: {}, mac: {}, floor {} )".format(self.getLat(),
                                                                                        self.getLong(), self.X, self.Y,
                                                                                        self.Mac, self.Floor)


class Range:
    def __init__(self, responder, time, distance):
        assert (isinstance(responder, Responder))
        self.Responder = responder
        self.Time = time
        self.Range = distance

    def __repr__(self):
        return "Range (Responder: {}, Time: {}, Range in meters: {}".format(self.Responder, self.Time, self.Range)


class RangeSet:
    def __init__(self, listOfRanges):
        assert (isinstance(listOfRanges, list))
        for r in listOfRanges:
            assert (isinstance(r, Range))
        self.RangeSet = listOfRanges

    def CreateCyclesList(self):
        cycles = []
        for r in self.RangeSet:
            c = CircleUtility.circleFactory(r.Responder.X, r.Responder.Y, r.Range)
            cycles.append(c)
        return cycles

    def CalculateIntersectionMatrix(self):
        dic = {}
        cycles = self.CreateCyclesList()
        for c1 in cycles:
            for c2 in cycles:
                if c1 == c2:
                    continue

                l = CircleUtility.circle_intersection_sympy(c1, c2)
                assert (isinstance(l, list))
                assert(all(isinstance(elem, Point) for elem in l))
                for p in l:
                    if (p.x, p.y) not in dic:
                        dic[(p.x,p.y)] = 1
                    else:
                        dic[(p.x, p.y)] += 1
        return dic


if __name__ == "__main__":
    r1 = Responder(0, 0, 0, 0, "a", 2)
    #r2 = Responder(0, 0, 0, 2, "b", 2)
    r3 = Responder(0, 0, 2, 0, "c", 2)

    range2 = Range(r1,0 ,1)
    range3 = Range(r3, 0, 1)
    #range1 = Range(r1, 0, math.sqrt(8))

    rangeSet = [range2, range3]
    rs = RangeSet(rangeSet)
    print(rs.CreateCyclesList())
    print(rs.CalculateIntersectionMatrix())

