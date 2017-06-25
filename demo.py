import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from mapper import IndoorMapper
from descartes import PolygonPatch
import networkx as nx
from shapely.geometry import Polygon
from ResponderLoader import FixFactory
import pickle
import time
import matplotlib.image as mpimg
from ImagePlot import NavigationTool
from shapely.geometry import Point


class Router(object):

    def __init__(self, map_polygon, routing_graph, figure, axs):
        self.is_select_source = True
        self.source = None
        self.destination = None
        self.map_poly = map_polygon
        self.routing_graph = routing_graph
        self.figure = figure
        self.axs = axs
        self.path = None

        ''' Support routing in any place in the map'''
        self.externalSource = None
        self.externalDest = None

    def select_source(self, event):
        print("selecting source")
        self.is_select_source = True

    def select_destination(self, event):
        print("selecting destination")
        self.is_select_source = False

    def route(self, event):
        print("Routing...")
        if self.source is None or self.destination is None:
            print("Select source and destination, then route...")
            return
        if self.source == self.destination:
            print("Cannot route from point to itself!")
            return
        else:
            try:
                path = nx.dijkstra_path(self.routing_graph, self.source, self.destination)
            except nx.NetworkXNoPath:
                print("Source and Destination are not connected!")
                return
            xs, ys = zip(*path)
            xsL = list(xs)
            ysL = list(ys)

            ''' Add external source and dest point to routing graph'''
            if self.externalSource is not None and self.externalDest is not None:
                xsL.insert(0, self.externalSource.x)
                xsL.append(self.externalDest.x)
                ysL.insert(0, self.externalSource.y)
                ysL.append(self.externalDest.y)

            line, = self.axs.plot(xsL, ysL, color='red')
            self.path = line
            self.figure.canvas.draw()

    def clear(self, event):
        self.path.set_data([],[])
        self.figure.canvas.draw()
        self.source = None
        self.destination = None

    def onpick(event, router_instance, axs):

        assert isinstance(router_instance, Router)
        print("router onpick")

        ''' Delete external points'''
        router_instance.externalSource = None
        router_instance.externalDest = None

        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        points = tuple(zip(xdata[ind], ydata[ind]))
        point = points[0]
        Router._onPick(point, router_instance)

    @staticmethod
    def _onPick(point, router_instance):
        if router_instance.is_select_source:
            router_instance.source = point
        else:
            router_instance.destination = point
        print("source: {}, dest: {}".format(router_instance.source, router_instance.destination))

    @staticmethod
    def externalPick(externalPoint, inGraphP, router_instance):
        assert isinstance(externalPoint, Point)
        assert isinstance(inGraphP, Point)
        assert isinstance(router_instance, Router)

        if router_instance.is_select_source:
            router_instance.externalSource = Point(externalPoint.x, externalPoint.y)
        else:
            router_instance.externalDest = Point(externalPoint.x, externalPoint.y)

        print("source: {}, dest: {}".format(router_instance.externalSource, router_instance.externalDest))

        Router._onPick((inGraphP.x, inGraphP.y), router_instance)

class Demo:
    def __init__(self, map_poly, routing_graph):
        fig, axs = plt.subplots()
        self.axs = axs
        self.figure = fig
        self.router = Router(map_poly, routing_graph, fig, axs)

        plt.subplots_adjust(bottom=0.2)
        axsrc = plt.axes([0.48, 0.05, 0.1, 0.075])
        axdst = plt.axes([0.59, 0.05, 0.1, 0.075])
        axsroute = plt.axes([0.7, 0.05, 0.1, 0.075])
        axsclear = plt.axes([0.81, 0.05, 0.1, 0.075])

        callback = self.router

        self.bsrc = Button(axsrc, 'Source', color='0.95')
        self.bsrc.on_clicked(self.select_source)

        self.bdst = Button(axdst, 'Dest', color='0.85')
        self.bdst.on_clicked(self.select_destination)

        self.broute = Button(axsroute, 'Route')
        self.broute.on_clicked(callback.route)

        self.bclear = Button(axsclear, 'Clear')
        self.bclear.on_clicked(callback.clear)

        self.buttonList = [self.bsrc, self.bdst, self.broute, self.bclear]

    def select_source(self, event):
        self.bsrc.color = '0.95'
        self.bdst.color = '0.85'
        self.router.select_source(event)

    def select_destination(self, event):
        self.bsrc.color = '0.85'
        self.bdst.color = '0.95'
        self.router.select_destination(event)

    def run(self):
        if type(self.router.map_poly) is not Polygon:
            for polygon in self.router.map_poly:
                x, y = polygon.exterior.xy
                self.axs.plot(x, y)
                patch = PolygonPatch(polygon, alpha=0.5, zorder=2)
                self.axs.add_patch(patch)
        else:
            x, y = self.router.map_poly.exterior.xy
            self.axs.plot(x, y)
            patch = PolygonPatch(self.router.map_poly, alpha=0.5, zorder=2)
            self.axs.add_patch(patch)

        image = mpimg.imread("arim.png")
        self.axs.imshow(image, extent=(-3, 112, -3, 82))

        xs = list(map(lambda nd: nd[0], self.router.routing_graph.nodes()))
        ys = list(map(lambda nd: nd[1], self.router.routing_graph.nodes()))
        line, = self.axs.plot(xs, ys, '.', picker=5)

        callback = self.router

        self.figure.canvas.mpl_connect('pick_event', lambda event: Router.onpick(event, callback, self.axs))

        self.NavTool = NavigationTool(xs, ys)
        self.figure.canvas.mpl_connect('button_release_event', lambda event: self.brelease(event, callback, self.buttonList))
        plt.show()

    def brelease(self, event, routing_instance, list_of_buttons):
        x = event.xdata
        y = event.ydata

        if x is None or y is None:
            print("Selected is None data in button_release")
            return

        selected_point = Point(x, y)

        ''' Work around to ignore button release events'''
        for bt in list_of_buttons:
            assert isinstance(bt, Button)
            if bt.ax == event.inaxes:
                print("Selected point inside button")
                return

        print("X: {} Y: {} released".format(x, y))
        self.axs.plot(x, y, '*')
        closestP = self.NavTool.GetClosestPoint(Point(x,y))
        assert isinstance(closestP, Point)
        self.axs.plot(closestP.x, closestP.y, 'o')
        self.figure.canvas.draw()

        Router.externalPick(selected_point, Point(closestP.x, closestP.y), routing_instance)

    def createDemo():
        # ------------------ Don't Delete! -----------------#
        # Good for TauTrajDump:
        # mapper = IndoorMapper('TauTrajDump', 12.5, 3)
        # -------------------------------------------------#
        csvName = 'ArimSim'
        mapper = IndoorMapper(csvName, 18, 3)
        map_poly, graph = mapper.run(None, False)
        d = Demo(map_poly, graph)
        try:
            with open("Demos/" + csvName + "_demo_" + time.strftime("%d-%m-%Y_%H-%M-%S", time.gmtime()), "wb") as pickleFile:
                pickle.dump((map_poly, graph), pickleFile)
        except Exception:
            print("fail dumping demo")
        d.run()

    def runDemo(pickeName):
        try:
            with open(pickeName, "rb") as pickleFile:
                map_poly, graph = pickle.load(pickleFile)
                d = Demo(map_poly, graph)
                d.run()
        except Exception:
            print("fail loading demo")

# Demo.createDemo()
if __name__ == "__main__":
    Demo.runDemo("Demos/ArimSim_demo")
# Demo.runDemo("TauTrajDump")
