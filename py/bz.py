import sympy as sy

import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import configparser as cnf

home = "../"
data = "../dat/"
png = "../fig/png/"
svg = "../fig/svg/"

config = cnf.ConfigParser()
config.read(home+'config.ini', encoding='utf-8')

physics = config['physics']
a  = float(physics.get('a'))
c  = float(physics.get('c'))
g0 = float(physics.get('g0'))
bands = int(physics.get('bands'))
bandsT = int(physics.get('bandsT'))
bandsL = int(physics.get('bandsL'))
lowest_T = int(physics.get('lowest_band_T'))
lowest_L = int(physics.get('lowest_band_L'))

colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

b = np.array([[float(physics.get('b1x')), float(physics.get('b1y')), float(physics.get('b1z'))], \
     [float(physics.get('b2x')), float(physics.get('b2y')), float(physics.get('b2z'))], \
     [float(physics.get('b3x')), float(physics.get('b3y')), float(physics.get('b3z'))]])

b = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

figure_size_unit = 6

x, y = np.linspace(-2,2,100), np.linspace(-2,2,100)
kx, ky = np.meshgrid(x, y)
px, pz = np.meshgrid(x, y)
qy, qz = np.meshgrid(x, y)

def get_plane(vec, plane): # {{{
    if vec[2] == 0:
        if vec[1] == 0:
            p, q, r = sy.symbols("p q r")
            eq = plane.equation(p, q, r)
            qx = sy.solve(eq, p)
            args = (q, r)
            qx = sy.lambdify(args, qx, 'numpy')
            xx = qx(qy, qz)
            yy = qy
            zz = qz
        else:
            p, q, r = sy.symbols("p q r")
            eq = plane.equation(p, q, r)
            py = sy.solve(eq, q)
            args = (p, r)
            plane = sy.lambdify(args, py, 'numpy')
            xx = px
            yy = plane(px, pz)[0]
            zz = pz
    else:
        p, q, r = sy.symbols("p q r")
        eq = plane.equation(p, q, r)
        kz = sy.solve(eq, r)
        args = (p, q)
        plane = sy.lambdify(args, kz, 'numpy')
        xx = kx
        yy = ky
        zz = plane(kx, ky)[0]

    return xx, yy, zz
# }}}

def get_line(line): # {{{
    xx = []
    yy = []
    zz = []

    t = sy.symbols("t")
    c = line.arbitrary_point(t)
    lx = sy.lambdify(t, c.x, 'numpy')
    ly = sy.lambdify(t, c.y, 'numpy')
    lz = sy.lambdify(t, c.z, 'numpy')
    ti = np.linspace(-5,5,100)
    for i in ti:
        xx.append(lx(i))
        yy.append(ly(i))
        zz.append(lz(i))

    return xx, yy, zz
# }}}

def get_segment(segment): # {{{
    xx = []
    yy = []
    zz = []

    t = sy.symbols("t")
    c = segment.arbitrary_point(t)

    lx = sy.lambdify(t, c.x, 'numpy')
    ly = sy.lambdify(t, c.y, 'numpy')
    lz = sy.lambdify(t, c.z, 'numpy')
    ti = np.linspace(0,1,20)
    for i in ti:
        xx.append(lx(i))
        yy.append(ly(i))
        zz.append(lz(i))

    return xx, yy, zz
# }}}

def get_point(point): # {{{
    xx, yy, zz = point.x, point.y, point.z
    return xx, yy, zz
# }}}

def perpendicular_plane(vec): # {{{
    p = vec / 2e0
    plane = sy.Plane(sy.Point3D(p), normal_vector=p)
#    print(plane)
    return plane
# }}}

def intersection_planes(planes): # {{{
    lines = []

    BZ_A = planes

    BZ_B = BZ_A.copy()
    for planeA in BZ_A:
        print(len(BZ_B))
        BZ_B.remove(planeA)
        for planeB in BZ_B:
            if planeA.is_parallel(planeB):
                continue

            line = planeA.intersection(planeB)[0]
            lines.append(line)
#            print(planeA, planeB, line)

    return lines
# }}}

def intersection_lines(lines): # {{{
    points = []

    tmp = lines.copy()
    for lineA in lines:
        print(len(tmp))
        tmp.remove(lineA)
        for lineB in tmp:
            if lineA.is_parallel(lineB):
                continue

            point = lineA.intersection(lineB)

            if len(point) == 0:
                continue
            points.append(point[0])

    return points
# }}}

def get_segments(points, lines): # {{{
    segments = []

    for line in lines:
        segment_points = []
        for point in points:
            if line.contains(point):
                segment_points.append(point)

        tmp = segment_points.copy()
        for pointA in segment_points:
            tmp.remove(pointA)
            for pointB in tmp:
                if pointB.equals(pointA):
                    continue
                s = sy.Segment3D(pointA, pointB)
                segments.append(s)

#        if len(segment_points) > 2:
#            origin = sy.Point3D(0, 0, 0)
#            

    return segments
# }}}

def get_BZ(): # {{{
    BZ = []

    n  = np.arange(0,2)
#    n  = np.arange(-1,2)
    for i in n:
        for j in n:
            for k in n:
                if [i,j,k] == [0,0,0]:
                    continue
                p = i*b[0] + j*b[1] + k*b[2]
                print(i, j, k, p)
                plane = perpendicular_plane(p)
                BZ.append([p, plane])

    planes = []
    for row in BZ:
        planes.append(row[1])

    lines = intersection_planes(planes)

    points = intersection_lines(lines)

    segments = get_segments(points, lines)

#    BZ = [points, segments, BZ]
    BZ = [points, lines, BZ]

    return BZ
# }}}

def plot_BZ(BZ): # {{{
    fig = plt.figure(figsize=(figure_size_unit, figure_size_unit))
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    ax.set_xlabel("$k_x$")
    ax.set_ylabel("$k_y$")
    ax.set_zlabel("$k_z$")
    ax.set_xlim(-2,2)
    ax.set_ylim(-2,2)
    ax.set_zlim(-2,2)

#    for point in BZ[0]:
#        xx, yy, zz = get_point(point)
#        ax.scatter(xx, yy, zz)

#    for segment in BZ[1]:
#        xx, yy, zz = get_segment(segment)
#        ax.scatter(xx, yy, zz)

#    for [p, plane] in BZ[2]:
#        print(plane)
#        xx, yy, zz = get_plane(p, plane)
#        ax.plot_surface(xx, yy, zz)

    plt.show()
# }}}

def main():
    BZ = get_BZ()
    plot_BZ(BZ)

if __name__ == "__main__":
    main()
