from delaunay_bowyer_watson import *


def main():

    # Generate a random group of points:
    number_of_points = 200
    x_values = [np.random.random()*100 for x in range(number_of_points)]
    y_values = [np.random.random()*100 for y in range(number_of_points)]
    points = [Point(x, y) for x, y in zip(x_values, y_values)] 

    # Compute the Delaunay triangulation using the Bowyer-Watson algorithm:
    delaunay = DelaunayTriangulation(points)
    delaunay.plot()
    
    return

if __name__ == "__main__":
    main()