import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from mapper import IndoorMapper
from descartes import PolygonPatch
import networkx as nx
from shapely.geometry import Polygon
from ResponderLoader import FixFactory

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
            line, = self.axs.plot(xs, ys)
            self.path = line
            self.figure.canvas.draw()

    def clear(self, event):
        self.path.set_data([],[])
        self.figure.canvas.draw()
        self.source = None
        self.destination = None

    def onpick(event, router_instance, axs):
        print("router onpick")
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ind = event.ind
        points = tuple(zip(xdata[ind], ydata[ind]))
        point = points[0]
        if router_instance.is_select_source:
            router_instance.source = point
        else:
            router_instance.destination = point
        print("source: {}, dest: {}".format(router_instance.source, router_instance.destination))


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

        xs = list(map(lambda nd: nd[0], self.router.routing_graph.nodes()))
        ys = list(map(lambda nd: nd[1], self.router.routing_graph.nodes()))
        line, = self.axs.plot(xs, ys, 'o', picker=5)

        callback = self.router
        self.figure.canvas.mpl_connect('pick_event', lambda event: Router.onpick(event, callback, self.axs))

        plt.show()


def main():
    print("Genrating Fixes")
    #FixFactory.Generate()

    print("Starting IndoorMappingAlgorithm")
    mapper = IndoorMapper('TauTrajDump', 5, 3)
    map_poly, graph = mapper.run(None, False)
    d = Demo(map_poly, graph)
    d.run()


if __name__ == "__main__":
    main()
