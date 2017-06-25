import matplotlib.pyplot as plt
from scipy.misc import imread
import matplotlib.image as mpimg
from shapely.geometry import Point
import sys

if __name__ == "__main__":
    fig, axs = plt.subplots()
    image = mpimg.imread("arim.png")
    plt.imshow(image)
    plt.show()


class NavigationTool:
    def __init__(self, X, Y):
        assert isinstance(X, list)
        assert isinstance(Y, list)
        self._points = list()
        for x,y in zip(X,Y):
            self._points.append(Point(x,y))
        assert len(self._points) > 0

    def GetClosestPoint(self, p):
        """ Returns closest point to current"""

        assert len(self._points) > 0
        assert isinstance(p, Point)
        minLenght = sys.float_info.max
        closestPoint = self._points[0]

        for pTemp in self._points:
            ''' Calculate distance to current p point'''
            assert isinstance(pTemp, Point)
            dt = p.distance(pTemp)

            ''' Replace closest point in case found point closer'''
            if dt <= minLenght:
                minLenght = dt
                closestPoint = pTemp

        return closestPoint
