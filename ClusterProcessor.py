# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 00:38:50 2017

@author: Eyal
"""

from descartes import PolygonPatch
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.cm as cm
from networkx import Graph

from shapely.geometry import LineString, GeometryCollection, Polygon, MultiPoint, MultiPolygon
from shapely.geometry import Point as shapelyPoint

from IMCObjects import LLine
from Projector import RotationalProjector
import copy
import math
import statistics
import itertools



class ClusterProcessor:
    '''
    works each time on a different single cluster
    works in xy-local
    '''

    def __init__(self, MinLns):
        self.localLLinesList = None
        self.MinLns = MinLns
        self.rotationalProjector = RotationalProjector()


    def plot_clusters_and_reprs(self, axs, polygonsAndReprTrajs):
        """ input is a zip of <polygon, traj> """
        if axs == None:
            fig, axs = plt.subplots()
            # -------- debug - clusters plotting -------------------
            fig.suptitle('Cluster Shapes and Representative Trajectories')
            plt.xlabel('local x coordinates [m]')
            plt.ylabel('local y coordinates [m]')
            # -------------------------------------------------------

        polygonsAndReprTrajs = list(polygonsAndReprTrajs)
        colors = iter(cm.rainbow(np.linspace(0, 1, len(polygonsAndReprTrajs))))
        for e in polygonsAndReprTrajs:
            ccolor = next(colors)

            # plot polygon
            polygon = e[0]
            xs, ys = polygon.exterior.xy
            axs.fill(xs, ys, alpha=0.5, fc=ccolor, ec='none')

            # plot representative trajectory
            reprTraj_xist_ylist_tuple = e[1]
            axs.plot(*reprTraj_xist_ylist_tuple, color=ccolor, linewidth=5)
        plt.show()

    def process(self, clusterList):
        """ returns a zip of <polygon, traj> """
        polygons = []  # polygon objects
        trajs = []  # tuple of (xlist,ylist)
        for i, cluster in enumerate(clusterList):
            self.loadLocalCluster_lsegList(cluster)
            reprAndWalls = self.calcRepresentativeLocal(0.4)
            rxs = [raw[0][0] for raw in reprAndWalls]
            rys = [raw[0][1] for raw in reprAndWalls]

            if len(rxs) > 1:
                # evaluated shape:
                wxs = [raw[1][0] for raw in reprAndWalls] + [raw[2][0] for raw in reprAndWalls][::-1]
                wxs += [rxs[0]]
                wys = [raw[1][1] for raw in reprAndWalls] + [raw[2][1] for raw in reprAndWalls][::-1]
                wys += [rys[0]]
                poly_bySigma = Polygon(list(zip(wxs, wys)))
                clusterpoly = poly_bySigma

                # # make convex hull:
                # hull = MultiPoint([shapelyPoint(seg.end1[0], seg.end1[1]) for seg in cluster] + \
                #                   [shapelyPoint(seg.end2[0], seg.end2[1]) for seg in cluster]).convex_hull
                #
                # if type(hull) == GeometryCollection or type(hull) == LineString:
                #     print("cluster {} had no good convex hull".format(i))
                # else:
                #     hxs, hys = hull.exterior.xy
                #     poly_byConvexHull = Polygon(list(zip(hxs, hys)))
                #     clusterpoly = poly_bySigma.intersection(poly_byConvexHull)
                #     if type(clusterpoly) != Polygon:
                #         clusterpoly = poly_bySigma
                trajs.append((rxs, rys,))
                polygons.append(clusterpoly)

            else:
                print("cluster {} had no good representative".format(i))
                continue

        return zip(polygons, trajs)

    def make_map_from_zip_polygons_trajs(self, axs, polys_repr_trajs):
        """
        makes a zip of polys and reprs into a tuple of
        (map polygon, underlying routing graph)
        :param polys_repr_trajs: a zip of polys and reprs
        :return: (map polygon, underlying routing graph)
        """

        map_poly = MultiPolygon()
        polys_repr_trajs = list(polys_repr_trajs)
        polys_repr_trajs.sort(key=lambda tup: tup[0].area, reverse=True)

        polys, trajs = zip(*polys_repr_trajs)  # traj is a (xlist, ylist)

        grp = Graph()
        bridges = []

        # building bridges inside this loop
        for poly_traj_tup_a, poly_traj_tup_b in itertools.combinations(polys_repr_trajs, 2):
            poly_a, poly_b = poly_traj_tup_a[0], poly_traj_tup_b[0]
            traj_a, traj_b = poly_traj_tup_a[1], poly_traj_tup_b[1]
            is_intersect = poly_a.intersects(poly_b)
            is_a_contained = poly_b.contains(poly_a)
            is_b_contained = poly_a.contains(poly_b)
            if is_a_contained or is_b_contained:
                print("special containment case - deal with later")
                continue
            if is_intersect:
                intersection_of_clusters = poly_a.intersection(poly_b)

                if type(intersection_of_clusters) is Polygon:  # Partial intersection
                    bridge = ClusterProcessor.bridge(traj_a, traj_b, intersection_of_clusters)
                    bridges.append(bridge)

                elif type(intersection_of_clusters) is MultiPolygon:
                    for polyg in intersection_of_clusters:
                        bridge = ClusterProcessor.bridge(traj_a, traj_b, polyg)
                        bridges.append(bridge)

                else:
                    print("the type of the intersection was {}".format(type(intersection_of_clusters)))

        # building map:
        for poly_traj_tup in polys_repr_trajs:
            poly = poly_traj_tup[0]
            traj = poly_traj_tup[1]

            # adding poly to map
            map_poly = map_poly.union(poly)

            # adding trajectory to map graph
            xytups = list(zip(*traj))
            grp.add_nodes_from(xytups)

            # add all edges of cluster to the map
            for i, tup in enumerate(xytups):
                if i == len(xytups) - 1:
                    break
                grp.add_edge(tup, xytups[i+1], weight=ClusterProcessor.dist(*tup, *xytups[i+1]))

        for a, m, b in bridges:
            grp.add_node(m)
            grp.add_edge(a, m, weight=ClusterProcessor.dist(*a, *m))
            grp.add_edge(m, b, wieght=ClusterProcessor.dist(*m, *b))

        # optional plotting:
        fig, axs = plt.subplots()
        fig.suptitle('Map and Underlying Routing Graph')
        plt.xlabel('local x coordinates [m]')
        plt.ylabel('local y coordinates [m]')
        if isinstance(map_poly, MultiPolygon):
            for polygon in map_poly:
                x, y = polygon.exterior.xy
                axs.plot(x, y)
                patch = PolygonPatch(polygon, alpha=0.5, zorder=2)
                axs.add_patch(patch)
            for edge in grp.edges():
                xs, ys = list(zip(*edge))
                axs.plot(xs, ys, color="b")
            plt.show()
        else:
            x, y = map_poly.exterior.xy
            axs.plot(x, y)
            patch = PolygonPatch(map_poly, alpha=0.5, zorder=2)
            axs.add_patch(patch)
            for edge in grp.edges():
                xs, ys = list(zip(*edge))
                axs.plot(xs, ys, color="b")
            plt.show()
        return map_poly, grp

    def dist(x1, y1, x2, y2):
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

    def bridge(traj_a, traj_b, polyg):
        """
        gets two trajectories from two clusters,
        and returns 3 points of a short bridge between them.
        :param traj_a: trajectory of cluster a
        :param traj_b: trajectory of cluster b
        :param polyg: the intersection of poly_a and poly_b
        :return: (pta, meanp, ptb)
        """
        ptsa = list(zip(*traj_a))
        ptsb = list(zip(*traj_b))

        x, y = polyg.exterior.xy
        mean_point = sum(x) / len(x), sum(y) / len(y)

        # find closest point in traj_a
        pts_and_distances = [(xytup, ClusterProcessor.dist(*xytup, *mean_point)) for xytup in ptsa]
        min_pt_and_dist_a = min(pts_and_distances, key=lambda pt: pt[1])
        # find closest point in traj_b
        pts_and_distances = [(xytup, ClusterProcessor.dist(*xytup, *mean_point)) for xytup in ptsb]
        min_pt_and_dist_b = min(pts_and_distances, key=lambda pt: pt[1])

        return (min_pt_and_dist_a[0], mean_point, min_pt_and_dist_b[0])

    def loadLocalCluster_lsegList(self, localClusterLlineList):
        self.localLLinesList = localClusterLlineList

    def calcAverageDirection(self):
        """
        assume all LLines are with bearing 0 <= bearing < 180
        :return: xytup_start xytup_end
        """

        if len(self.localLLinesList) < 1:
            raise ValueError("calcMainDirection: cluster has to have segments!")
        origin = self.localLLinesList[0].end1  # (x,y) tuple

        currPosition = origin
        for lline in self.localLLinesList:
            nextPosition_a = lline.llineTravel(currPosition)
            nextPosition_b = lline.llineTravel(currPosition, reverse=True)

            dist_a = math.sqrt((nextPosition_a[0] - origin[0]) ** 2 +
                               (nextPosition_a[1] - origin[1]) ** 2)
            dist_b = math.sqrt((nextPosition_b[0] - origin[0]) ** 2 +
                               (nextPosition_b[1] - origin[1]) ** 2)
            currPosition = nextPosition_a if dist_a > dist_b else nextPosition_b

        totalWay = LLine(origin, currPosition)

        originCopy = copy.deepcopy(origin)
        avgDirEnd2 = totalWay.llineTravel(originCopy, 2 / len(self.localLLinesList))

        avgDir = LLine(originCopy, avgDirEnd2)
        return avgDir

    def avgDirPhi(avgDirEnd1, avgDirEnd2):
        dx = avgDirEnd2[0] - avgDirEnd1[0]
        dy = avgDirEnd2[1] - avgDirEnd1[1]
        if dx == 0 and dy == 0:
            raise Exception("calcRepresentativeLocal: main dir of cluster is a dot!")
        elif dx == 0 and dy > 0:
            phi = 90
        elif dx == 0 and dy < 0:
            phi = -90
        elif dy == 0:
            phi = 0
        elif dy > 0:
            phi = math.degrees(math.atan(dy / dx))
        else:
            phi = - math.degrees(math.atan(-dy / dx))
        return phi

    def rotateLLineCoords(self, lline):
        e1x, e1y = lline.end1
        e1xytag = self.rotationalProjector.xyToRotatedXY(e1x, e1y)

        e2x, e2y = lline.end2
        e2xytag = self.rotationalProjector.xyToRotatedXY(e2x, e2y)
        return LLine(e1xytag, e2xytag)

    def calcRepresentativeLocal(self, gamma):
        """
        calculates the representative trajectory for the cluster
        currently loaded to the processor
        """

        def pltLLines(axs, llinesList, color='r'):
            for lline in llinesList:
                xs = [lline.end1[0], lline.end2[0]]
                ys = [lline.end1[1], lline.end2[1]]
                axs.plot(xs, ys, color=color)

        avgDir = self.calcAverageDirection()
        avgDirEnd1, avgDirEnd2 = avgDir[0], avgDir[1]

        phi = ClusterProcessor.avgDirPhi(avgDirEnd1, avgDirEnd2)
        self.rotationalProjector.rotate(phi)

        rotLLinesList = [self.rotateLLineCoords(lline) for lline in self.localLLinesList]
        rotLLinesDxYy = [(lline.end2[0] - lline.end1[0], lline.end2[1] - lline.end1[1]) for lline in rotLLinesList]

        # sorting llines once by starting points and once by ending points
        startingPoint = lambda l: l.end1[0]
        endingPoint = lambda l: l.end2[0]

        sortedByStartings = sorted(rotLLinesList, key=startingPoint)
        sortedByEndings = sorted(rotLLinesList, key=endingPoint)

        calculationsPool = []  # list of elements of format [x-val, [y-vals]]
        xval = float("-inf")
        while len(sortedByStartings) > 0 or len(sortedByEndings) > 0:
            if len(sortedByStartings) == 0 and len(sortedByEndings) > 0:
                a_starting_point = False
            elif len(sortedByEndings) == 0 and len(sortedByStartings) > 0:
                a_starting_point = True
            else:
                if sortedByStartings[0].end1[0] == sortedByEndings[0].end2[0]:
                    sortedByEndings = sortedByEndings[1:]
                    a_starting_point = True
                else:
                    a_starting_point = sortedByStartings[0].end1[0] <= sortedByEndings[0].end2[0]

            if a_starting_point:
                currLLine = sortedByStartings[0]
                sortedByStartings = sortedByStartings[1:]
                currPt = currLLine.end1
            else:
                currLLine = sortedByEndings[0]
                sortedByEndings = sortedByEndings[1:]
                currPt = currLLine.end2

            if currPt[0] < xval + gamma:
                continue

            xval = currPt[0]
            yvals = []
            # add to not look at too many segments
            for i, lline in enumerate(rotLLinesList):
                if not (lline.end1[0] <= xval and xval <= lline.end2[0]):
                    continue

                dx, dy = rotLLinesDxYy[i]
                if dx == 0 and lline.end1[0] == xval:
                    # edge case exactly perpendicular to axis
                    yvals.append((lline.end1[1] + lline.end2[1]) / 2)
                elif dx != 0:
                    # need to find the intersection with the line.
                    m = dy / dx
                    u, v = lline.end2
                    y = m * (xval - u) + v
                    yvals.append(y)

            if len(yvals) < 2:  # self.MinLns:
                continue

            calculationsPool.append([xval, statistics.mean(yvals), (statistics.stdev(yvals))])

        representative_and_walls = [(self.rotationalProjector.rotatedXYToXY(pe[0], pe[1]), \
                                     self.rotationalProjector.rotatedXYToXY(pe[0], pe[1] +   math.sqrt(3) * pe[2]), \
                                     self.rotationalProjector.rotatedXYToXY(pe[0], pe[1] -   math.sqrt(3) * pe[2]))
                                    for pe in calculationsPool]

        if len(representative_and_walls) < 2:
            return []
        return representative_and_walls
