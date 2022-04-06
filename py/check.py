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

numeric = config['numeric']

colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

readfile = '../triangle_L1_vertex.csv'
#T = pd.read_csv(readfile,header=None).values
#T_kx_max, T_kx_min = np.max(T[:,0]), np.min(T[:,0])
#T_ky_max, T_ky_min = np.max(T[:,1]), np.min(T[:,1])
#T_kz_max, T_kz_min = np.max(T[:,2]), np.min(T[:,2])
#T_len = max(T_kx_max - T_kx_min, T_ky_max - T_ky_min, T_kz_max - T_kz_min)
#T_center_x = np.mean(T[:,0])
#T_center_y = np.mean(T[:,1])
#T_center_z = np.mean(T[:,2])
#T_window = [[T_center_x - T_len/2, T_center_x + T_len/2],
#            [T_center_y - T_len/2, T_center_y + T_len/2],
#            [T_center_z - T_len/2, T_center_z + T_len/2]]

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title("T")
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("kz")
#ax.set_xlim(T_window[0])
#ax.set_ylim(T_window[1])
#ax.set_zlim(T_window[2])
readfile = '../triangle_L1_face.csv'
#readfile = '../check.csv'
T = pd.read_csv(readfile,header=None)
T = T.groupby((T.isnull().all(axis=1)).cumsum())
for index, g in T:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax.plot(g[:,0],  g[:,1],  g[:,2])
plt.show()
