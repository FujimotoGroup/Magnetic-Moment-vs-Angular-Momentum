import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.spatial import Delaunay
import open3d as o3d

if __name__ == "__main__":
    for i in np.arange(5):
        readfile = '../dat/L0_band0k'+str(i)+'.csv'
        readfile = '../dat/T_band4k'+str(i)+'.csv'
        position = pd.read_csv(readfile, header=None).values
        x = position[:,0]
        y = position[:,1]
        z = position[:,2]

        fig = plt.figure(figsize=plt.figaspect(0.5))
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        ax.scatter(x, y, z, label="e = "+str(position[0,6]))
        ax.legend()
        plt.show()

        point_cloud =  pd.read_csv(readfile, header=None).values
        pcd = o3d.geometry.PointCloud()
        pcd.points = o3d.utility.Vector3dVector(point_cloud[:,0:3])
        pcd.normals = o3d.utility.Vector3dVector(point_cloud[:,3:6])
#        o3d.visualization.draw_geometries([pcd])

        convex_mesh = pcd.compute_convex_hull()[0]
        print("convex_mesh is water-tight: ", convex_mesh.is_watertight())
        o3d.visualization.draw_geometries([pcd, convex_mesh])

        radii = [0.005, 0.01, 0.02, 0.04]
        rec_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting(
                    pcd, o3d.utility.DoubleVector(radii))
        print("rec_mesh is water-tight: ", rec_mesh.is_watertight())
        poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
        print("poisson_mesh is water-tight: ", poisson_mesh.is_watertight())
        o3d.visualization.draw_geometries([pcd, rec_mesh])
        o3d.visualization.draw_geometries([pcd, poisson_mesh])
        print(poisson_mesh)

#for i in np.arange(5):
#    readfile = '../dat/L0_band0k'+str(i)+'.csv'
#    readfile = '../dat/T_band4k'+str(i)+'.csv'
#    position = pd.read_csv(readfile, header=None).values
#    x = position[:,0]
#    y = position[:,1]
#    z = position[:,2]
#
#    plt.rcParams["axes.facecolor"] = 'white'
#    fig = plt.figure(figsize=plt.figaspect(0.5))
#
#    ax = fig.add_subplot(1, 2, 1, projection='3d')
#    ax.scatter(x, y, z, label="e = "+str(position[0,6]))
#    ax.legend()
#    plt.show()
