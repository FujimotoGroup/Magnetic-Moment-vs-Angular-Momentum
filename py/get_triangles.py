import sys
import numpy as np
import pandas as pd
import open3d as o3d

args = sys.argv

if __name__ == "__main__":
    readfile = args[1]
    outputfile = args[2]
    point_cloud =  pd.read_csv(readfile, header=None).values
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(point_cloud[:,0:3])
    pcd.normals = o3d.utility.Vector3dVector(point_cloud[:,3:6])
    poisson_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
    o3d.io.write_triangle_mesh(filename=outputfile, mesh=poisson_mesh, write_vertex_colors=False, write_triangle_uvs=False)
