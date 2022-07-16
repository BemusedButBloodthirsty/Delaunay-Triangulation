from matplotlib import pyplot as plt
import numpy as np


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


class Triangle:
    def __init__(self, p1: Point, p2: Point, p3: Point):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.orient_counterclockwise()
        self.edges = self.determine_edges()
        return


    def points(self):
        return [self.p1, self.p2, self.p3]


    def determine_edges(self):
        """
        Simply calculates the midpoint between the 3 pairs of points:
        
        Edge 1 => (P1, P2)
        Edge 2 => (P2, P3)
        Edge 3 => (P3, P1)
        """
        def midpoint(p1: Point, p2: Point):
            return Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)

        edges = []
        edges.append(midpoint(self.p1, self.p2))
        edges.append(midpoint(self.p2, self.p3))
        edges.append(midpoint(self.p3, self.p1))
        return edges


    def plot(self, points=None, x_limits=(0,10), y_limits=(0,10)):
        plt.figure()
        plt.xlim(x_limits)
        plt.ylim(y_limits)

        plt.plot([self.p1.x, self.p2.x, self.p3.x, self.p1.x],[self.p1.y, self.p2.y, self.p3.y, self.p1.y], '*-')
        if points:
            plt.plot([p.x for p in points], [p.y for p in points], '.k')
        
        plt.show()


    def is_linear(self):
        return abs((self.p2.y - self.p1.y)*(self.p3.x - self.p2.x) - (self.p3.y - self.p2.y)*(self.p2.x - self.p1.x)) <= 1e-3 


    def orient_counterclockwise(self):
        val = (self.p2.y - self.p1.y)*(self.p3.x - self.p2.x) - (self.p3.y - self.p2.y)*(self.p2.x - self.p1.x)
        if val > 0:
            # Orient the tri counterclockwise by swapping p2 and p3.
            temp = self.p2
            self.p2 = self.p3
            self.p3 = temp
            print(f"Orientation of {self} has been corrected to be counterclockwise.")
        return 


    def circumcircle_contains_point(self, p: Point):
        Ax = self.p1.x
        Bx = self.p2.x
        Cx = self.p3.x
        Dx = p.x

        Ay = self.p1.y
        By = self.p2.y
        Cy = self.p3.y
        Dy = p.y

        # Alternate derivation:
        # mat = np.array([
        #     [Ax, Ay, Ax**2 + Ay**2, 1],
        #     [Bx, By, Bx**2 + By**2, 1],
        #     [Cx, Cy, Cx**2 + Cy**2, 1],
        #     [Dx, Dy, Dx**2 + Dy**2, 1]
        # ])

        mat = np.array([
            [Ax - Dx, Ay - Dy, (Ax**2 - Dx**2) + (Ay**2 - Dy**2)],
            [Bx - Dx, By - Dy, (Bx**2 - Dx**2) + (By**2 - Dy**2)],
            [Cx - Dx, Cy - Dy, (Cx**2 - Dx**2) + (Cy**2 - Dy**2)]
        ])

        return np.linalg.det(mat) > 0


class DelaunayTriangulation:
    def __init__(self, points):
        self.points = points
        self.triangles = []
        self.triangulate()


    def triangulate(self):
        """
        Determines the Delaunay Triangulation for a set of points according to the Bowyer-Watson Algorithm.
        
        Ref: https://en.wikipedia.org/wiki/Bowyer%E2%80%93Watson_algorithm
        """


        def points_overlap(p1: Point, p2: Point):
            return ((p1.x - p2.x)**2 + (p1.y - p2.y)**2)**0.5 <= 1e-6


        def find_solitary_edges(bad_triangles):
            solitary_edges = []

            # Brute-force method to determine conflicts:
            for tri_i_index, tri_i in enumerate(np.array(self.triangles)[bad_triangles]):
                for edge_midpoint_i_index, edge_midpoint_i in enumerate(tri_i.edges):
                    no_conflicts = True
                    
                    for tri_j_index, tri_j in enumerate(np.array(self.triangles)[bad_triangles]):
                        if tri_i_index != tri_j_index: # Making sure we don't compare the same triangle.
                            for edge_midpoint_j_index, edge_midpoint_j in enumerate(tri_j.edges):
                                
                                if points_overlap(edge_midpoint_i, edge_midpoint_j):
                                    no_conflicts = False

                    if no_conflicts:
                        p1, p2, p3 = tri_i.p1, tri_i.p2, tri_i.p3

                        if edge_midpoint_i_index == 0:
                            solitary_edges.append((p1, p2))
                        elif edge_midpoint_i_index == 1:
                            solitary_edges.append((p2, p3))    
                        elif edge_midpoint_i_index == 2:
                            solitary_edges.append((p3, p1))
                
            return solitary_edges


        def remove_bad_triangles(bad_triangles):
            bad_triangles.sort(reverse=True)
            for bad_tri in bad_triangles:
                del self.triangles[bad_tri]
            return


        # Generating the "super triangle":
        x_max = max([p.x for p in self.points])
        y_max = max([p.y for p in self.points])
        x_min = min([p.x for p in self.points])
        y_min = min([p.y for p in self.points])

        s1 = Point((x_max + x_min) / 2, 2*y_max)
        s2 = Point(((x_max + x_min) / 2) + ((x_max - x_min) / 2)*((2*y_max - y_min) / y_max) + abs(x_max - x_min), y_min/2)
        s3 = Point(((x_max + x_min) / 2) - ((x_max - x_min) / 2)*((2*y_max - y_min) / y_max) - abs(x_max - x_min), y_min/2)
        super_triangle = Triangle(s1, s2, s3)
        self.triangles.append(super_triangle)
        # Debug plot => super_triangle.plot(self.points, x_limits=(s3.x,s2.x), y_limits=(s2.y, s1.y))

        # Inserting each point for triangulation:
        for p in self.points:
            
            # Determine the bad triangles:
            bad_triangles = []
            for i, tri in enumerate(self.triangles):

                if tri.circumcircle_contains_point(p):
                    # This is a triangle that needs to be deleted.
                    bad_triangles.append(i)
            
            # Determine the bad polygon:
            polygon = find_solitary_edges(bad_triangles) 

            # Remove all bad triangles:
            remove_bad_triangles(bad_triangles)

            # Create new triangles from the solitary edges stored in polygon:
            for solitary_edge in polygon:
                new_tri = Triangle(p, solitary_edge[0], solitary_edge[1])
                self.triangles.append(new_tri)
                # Debug plotting => new_tri.plot(x_limits=(s3.x,s2.x), y_limits=(s2.y, s1.y))

        # Debug plot => self.plot()
        
        # Remove all triangles that contain any of the vertices of the "super triangle":
        bad_triangles = []
        for i, tri in enumerate(self.triangles):
            for point in tri.points():

                if points_overlap(point, s1) or points_overlap(point, s2) or points_overlap(point, s3):
                    bad_triangles.append(i)
                    break

        remove_bad_triangles(bad_triangles)

        return 


    def plot(self):
        plt.figure()

        for tri in self.triangles:
            plt.plot([tri.p1.x, tri.p2.x, tri.p3.x, tri.p1.x], [tri.p1.y, tri.p2.y, tri.p3.y, tri.p1.y], '-k')
        
        plt.plot([p.x for p in self.points], [p.y for p in self.points], '.r')

        plt.show()
        return

