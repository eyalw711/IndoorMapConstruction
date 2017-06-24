import numpy as np
from matplotlib import pyplot as plt


# code converted from https://github.com/subpos/subpos_receiver/blob/master/Trilateration.cpp
# RespLocDic - Dictionary with the location of each responder in cartesian coordinate
# ResultDic - Dictionary with ranges of user from each responder (at least 3 valid ranges)
def trilateration(RespLocDic, ResultDic):
    trilObjectList = TrilObject.TriAdapter(RespLocDic, ResultDic)
    return ImprovedTrilateration(trilObjectList)


def ImprovedTrilateration(trilObjectList):
    A = []
    b = []

    if len(trilObjectList) < 3:
        raise Exception("ResultDic should contain at least 3 values")

    cnt = len(trilObjectList)
    rows = (cnt * (cnt-1)) / 2

    for i1 in range(cnt):
        beacon1 = trilObjectList[i1]
        assert isinstance(beacon1, TrilObject)
        for i2 in range (i1+1,cnt):
            beacon2 = trilObjectList[i2]
            assert isinstance(beacon2, TrilObject)
            x1 = beacon1.X
            y1 = beacon1.Y
            r1 = beacon1.Distance

            x2 = beacon2.X
            y2 = beacon2.Y
            r2 = beacon2.Distance

            b1 = ((pow(x1, 2) - pow(x2, 2)) +
                  (pow(y1, 2) - pow(y2, 2)) -
                  (pow(r1, 2) - pow(r2, 2))) / 2
            b.append([b1])
            A1 = [x1 - x2, y1 - y2]
            A.append(A1)
    e = np.linalg.lstsq(A, b)
    xLoc, yLoc = e[0]
    return xLoc,yLoc


class TrilObject:
    def __init__(self, x, y, distance):
        self.X = x
        self.Y = y
        self.Distance = distance

    def __repr__(self):
        return "TrilObj: X:{} Y:{} Distance{}".format(self.X, self.Y, self.Distance)

    ''' Return List of TrilObject'''
    @staticmethod
    def TriAdapter(RespLocDic, ResultDic):
        trilObjectList = []
        for key, value in ResultDic.items():
            d = value
            responder = RespLocDic[key]
            xTemp = responder[0]
            yTemp = responder[1]
            trilObjectList.append(TrilObject(xTemp, yTemp, d))
        return trilObjectList


if __name__ == "__main__":
    respLocDic = {"r1": [0, 0, 0], "r2": [2, 0, 0], "r3": [1, 1, 0]}
    resultDic = {"r1": 1, "r2": 1, "r3": 1}
    x, y = trilateration(respLocDic, resultDic)

    fig, axs = plt.subplots(1, 1)
    axs.plot(x,y, 'ro')
    plt.show()