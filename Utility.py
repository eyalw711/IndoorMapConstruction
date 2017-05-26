from sympy.geometry import Circle, Point


class CircleUtility:
    @staticmethod
    def circle_intersection_sympy(circle1, circle2):
        c1 = circle1
        c2 = circle2

        # Return list[]
        return c1.intersection(c2)

    @staticmethod
    def circleFactory(x, y, r):
        return Circle(Point(x, y), r)

if __name__ == "__main__":
    circleUtility = CircleUtility()
    points = circleUtility.circle_intersection_sympy(CircleUtility.circleFactory(0, 0, 1),
                                                     CircleUtility.circleFactory(0, 1, 1))

    for p in points:
        assert(isinstance(p, Point))
    print(points)
