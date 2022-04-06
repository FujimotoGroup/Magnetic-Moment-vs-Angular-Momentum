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
    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=False)[0]
#    mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(pcd, depth=8, width=0, scale=1.1, linear_fit=True)[0]
#    print("poisson_mesh is water-tight: ", mesh.is_watertight())
#    mesh = pcd.compute_convex_hull()[0]
#    print("convex_mesh is water-tight: ", mesh.is_watertight())
    print("surface area = ", mesh.get_surface_area())
#    with open("surface_area.dat", mode='a') as f:
#        f.write(str(mesh.get_surface_area())+"\n")
#    print("is intersectioning = ", mesh.is_self_intersecting())
    o3d.io.write_triangle_mesh(filename=outputfile, mesh=mesh, write_vertex_colors=False, write_triangle_uvs=False)
