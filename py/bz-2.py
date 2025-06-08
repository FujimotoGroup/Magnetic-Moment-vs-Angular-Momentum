from scipy.spatial import Voronoi, Delaunay, ConvexHull
import numpy as np
import mpl_toolkits.mplot3d as a3
from matplotlib.colors import LightSource

import matplotlib as mpl
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
pdf = "../fig/pdf/"

mpl.use('Agg')
#mpl.use('pgf')
plt.rcParams.update({
    "lines.linewidth": 2.0,
    "font.size": 24,
    "text.usetex": True,
#    "font.family": "Times New Roman",
    "font.family": "Helvetica",
    "mathtext.fontset": 'cm',
#    "pgf.texsystem": "lualatex",
#    "pgf.rcfonts": False,    # don't setup fonts from rc parameters
#    "pgf.preamble": "\n".join([
#        r"\usepackage{url}",            # load additional packages
#        r"\usepackage{unicode-math}",   # unicode math setup
#        r"\usepackage{newcomputermodern}",   # unicode math setup
#        r"\setmainfont{NewCM10-Book}",  # serif font via preamble
#     ])
})
plt.rc('text.latex', preamble=r'\usepackage{bm}')

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

figure_size_unit = 10

b = np.array([[float(physics.get('b1x')), float(physics.get('b1y')), float(physics.get('b1z'))], \
     [float(physics.get('b2x')), float(physics.get('b2y')), float(physics.get('b2z'))], \
     [float(physics.get('b3x')), float(physics.get('b3y')), float(physics.get('b3z'))]])

b_abs = [np.sqrt(np.sum(b[i]*b[i])) for i in np.arange(3)]

def main():
    points = []
    n  = np.arange(-1,2)
    for i in n:
        for j in n:
            for k in n:
                p = i*b[0] + j*b[1] + k*b[2]
#                print(i, j, k, p)
                points.append((p[0], p[1], p[2]))

    vor = Voronoi(points) # Voronoi diagram計算

#    # -------------------------------------------------------------------
#
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.set_xlabel("$k_x$")
#    ax.set_ylabel("$k_y$")
#    ax.set_zlabel("$k_z$")
#
#    ax.set_title('k-points of interest')
#    x, y, z = vor.points[:, 0], vor.points[:, 1], vor.points[:, 2]
#    ax.plot(x, y, z, "o", ms=5, mew=0.5)
#    plt.show()
#
#    # -------------------------------------------------------------------
#
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.set_xlabel("$k_x$")
#    ax.set_ylabel("$k_y$")
#    ax.set_zlabel("$k_z$")
#    ax.set_title('Vertices of Voronoi cells')
#
#    x, y, z = vor.vertices[:, 0], vor.vertices[:, 1], vor.vertices[:, 2]
#    # excluding vertices too far from ROI
#    roi = (np.abs(x) < 20) & (np.abs(y) < 20)
#    ax.plot(x[roi], y[roi], z[roi], "o", color="#00cccc", ms=4, mew=0.5)
#    plt.show()
#
#    # -------------------------------------------------------------------
#
#    fig = plt.figure()
#    ax = fig.add_subplot(111, projection='3d')
#    ax.set_xlabel("$k_x$")
#    ax.set_ylabel("$k_y$")
#    ax.set_zlabel("$k_z$")
#    ax.set_title('Pairs of atoms between which each Voronoi ridge plane lies')
#
#    for pair in vor.points[vor.ridge_points]:
#        ax.plot(pair[:, 0], pair[:, 1], pair[:, 2], color='C0')
#    plt.show()
#
#    # -------------------------------------------------------------------

    fig = plt.figure(figsize=(figure_size_unit,figure_size_unit))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_axis_off()
#    ax.set_xlabel("$k_x$")
#    ax.set_ylabel("$k_y$")
#    ax.set_zlabel("$k_z$")
#    ax.set_aspect('equal')

    window = [-1.3, 1.3]
    ax.set_xlim(window)
    ax.set_ylim(window)
    ax.set_zlim(window)

    a = 1.2
    m = 0.1
    ax.plot([-a, a], [ 0, 0], [ 0, 0], color=colors[11])
    ax.plot([ 0, 0], [-a, a], [ 0, 0], color=colors[11])
    ax.plot([ 0, 0], [ 0, 0], [-a, a], color=colors[11])
    ax.text(0, 0, a, "Trigonal",  ha='left', va='bottom')
#    ax.text(0, 0, a, "Trigonal ($k_z$)",  ha='left', va='bottom')
    ax.text(a, 0, 0, "Binary",  ha='right', va='bottom')
#    ax.text(a, 0, 0, "($k_x$)", ha='right', va='top')
    ax.text(-0.01, a+m, 0, "Bisectrix", ha='left', va='bottom')
#    ax.text(-0.01, a+m, 0, "($k_y$)",   ha='left', va='top', fontsize="small")

    elev = 15
    azim = 59

    ls = LightSource(azdeg=-azim, altdeg=90)

    ridgenum = 0
    for i in np.array(vor.ridge_vertices):
        # exclude ridge planes with points outside ROI
        if -1 not in i:
            vertices = vor.vertices[i]
            kx, ky, kz = vertices[:, 0], vertices[:,1], vertices[:, 2]
            c = 1.2
            if (np.abs(kx).max() <= b_abs[0]/c) & (np.abs(ky).max() <= b_abs[1]/c) & (np.abs(kz).max() <= b_abs[2]/c):
                ridgenum += 1
#                poly=a3.art3d.Poly3DCollection([vertices], alpha=0.1, edgecolors=colors[11], facecolors=colors[4], shade=True, lightsource=ls)
                poly=a3.art3d.Poly3DCollection([vertices], alpha=0.1, edgecolors=colors[11], facecolors=colors[4])
                ax.add_collection3d(poly)

    name = data+'triangle_T-mu0e0-rough'
    readfile = name+'_vertex.csv'
    Tf = pd.read_csv(readfile,header=None)
    Tf = Tf.groupby((Tf.isnull().all(axis=1)).cumsum())
    for index, g in Tf:
        g = g.dropna()
        if len(g) > 0:
            g = g.values[:,:-1]
            ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[7])
    ax.text(float(physics.get('kTx')), float(physics.get('kTy'))+0.05, float(physics.get('kTz')), '$h$', ha='left', va='bottom')

    offset = [[  0, 0, 0.1 ], \
              [  0, 0, 0.05], \
              [0.0, 0, 0.05]]
    label_L = ["$e_2$" , "$e_3$", "$e_1$"]
    for i in np.arange(1,4):
        name = data+'triangle_L'+str(i)+'-mu0e0-rough'
        readfile = name+'_face.csv'
        Lf = pd.read_csv(readfile,header=None)
        Lf = Lf.groupby((Lf.isnull().all(axis=1)).cumsum())
        for index, g in Lf:
            g = g.dropna()
            if len(g) > 0:
                g = g.values
                ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[4])

        ax.text(float(physics.get('kL'+str(i)+'x'))+offset[i-1][0], \
                float(physics.get('kL'+str(i)+'y'))+offset[i-1][1], \
                float(physics.get('kL'+str(i)+'z'))+offset[i-1][2], label_L[i-1], ha='left', va='bottom')

    fig.tight_layout()
    ax.view_init(elev=elev, azim=azim)

    plt.savefig(pdf+"fermi_surface.pdf", dpi=300)
#    plt.rc("svg", fonttype="none")
#    plt.savefig(svg+"Fermi_surface.svg")

    plt.show()
#    print('No. of all ridge planes:', len(vor.ridge_vertices))
#    print('No. of ridge planes displayed:', ridgenum)
    plt.close()

# T {{{
    name = data+'triangle_T-mu0e0-rough'
    readfile = name+'_vertex.csv'
    T = pd.read_csv(readfile,header=None).values
    T_kx_max, T_kx_min = np.max(T[:,0]), np.min(T[:,0])
    T_ky_max, T_ky_min = np.max(T[:,1]), np.min(T[:,1])
    T_kz_max, T_kz_min = np.max(T[:,2]), np.min(T[:,2])
    T_len = max(T_kx_max - T_kx_min, T_ky_max - T_ky_min, T_kz_max - T_kz_min)
    T_center_x = np.mean(T[:,0])
    T_center_y = np.mean(T[:,1])
    T_center_z = np.mean(T[:,2])
    T_window = [[T_center_x - T_len/2, T_center_x + T_len/2],
                [T_center_y - T_len/2, T_center_y + T_len/2],
                [T_center_z - T_len/2, T_center_z + T_len/2]]

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_axis_off()
#    ax.set_title("T")
#    ax.set_xlabel("kx")
#    ax.set_ylabel("ky")
#    ax.set_zlabel("kz")
    ax.set_xlim(T_window[0])
    ax.set_ylim(T_window[1])
    ax.set_zlim(T_window[2])
    readfile = name+'_face.csv'
    T = pd.read_csv(readfile,header=None)
    T = T.groupby((T.isnull().all(axis=1)).cumsum())
    for index, g in T:
        g = g.dropna()
        if len(g) > 0:
            g = g.values
            ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[7], linewidth=1)

    fn = index
    str_vertex_num = "face num = "+str(fn)
    ax_pos = ax.get_position()
    fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_vertex_num)

    plt.savefig(pdf+"fermi_surface_T.pdf", bbox_inches = 'tight', dpi=300)
#    plt.show()
    plt.close()
# }}}


if __name__ == "__main__":
    main()
