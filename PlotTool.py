from matplotlib import pyplot as plt

from Projector import EquirectangularProjector


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

