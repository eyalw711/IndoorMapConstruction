import pandas as pd

from IMCObjects import GPoint
from Projector import EquirectangularProjector
from Structs import Responder


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
            Coords.append(Responder(row['lat'], row['long'],0,0,row['mac'],0))
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

if __name__ == "__main__":
    loader = ResponderLoader("Respodenrs", 2)
    loader.InitResponders()
    for rep in loader.RespondersCollection:
        print(rep)