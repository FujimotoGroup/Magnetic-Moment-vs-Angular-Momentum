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

# {{{

# T {{{
name = data+'triangle_T-mu0e0'
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
ax.set_title("T")
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("kz")
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
        ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[9])

fn = index
str_vertex_num = "face num = "+str(fn)
ax_pos = ax.get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_vertex_num)

plt.savefig(png+"Fermi_surface_T.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"Fermi_surface_T.svg")
#plt.show()
plt.close()
# }}}

# L1 {{{
name = data+'triangle_L1-mu0e0'
readfile = name+'_vertex.csv'
L1 = pd.read_csv(readfile,header=None).values
L1_kx_max, L1_kx_min = np.max(L1[:,0]), np.min(L1[:,0])
L1_ky_max, L1_ky_min = np.max(L1[:,1]), np.min(L1[:,1])
L1_kz_max, L1_kz_min = np.max(L1[:,2]), np.min(L1[:,2])
L1_len = max(L1_kx_max - L1_kx_min, L1_ky_max - L1_ky_min, L1_kz_max - L1_kz_min)
L1_center_x = np.mean(L1[:,0])
L1_center_y = np.mean(L1[:,1])
L1_center_z = np.mean(L1[:,2])
L1_window = [[L1_center_x - L1_len/2, L1_center_x + L1_len/2],
             [L1_center_y - L1_len/2, L1_center_y + L1_len/2],
             [L1_center_z - L1_len/2, L1_center_z + L1_len/2]]

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title("L1")
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("kz")
ax.set_xlim(L1_window[0])
ax.set_ylim(L1_window[1])
ax.set_zlim(L1_window[2])
readfile = name+'_face.csv'
L1 = pd.read_csv(readfile,header=None)
L1 = L1.groupby((L1.isnull().all(axis=1)).cumsum())
for index, g in L1:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[4])

fn = index
str_vertex_num = "face num = "+str(fn)
ax_pos = ax.get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_vertex_num)

plt.savefig(png+"Fermi_surface_L1.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"Fermi_surface_L1.svg")
#plt.show()
plt.close()
# }}}

# L2 {{{
name = data+'triangle_L2-mu0e0'
readfile = name+'_vertex.csv'
L2 = pd.read_csv(readfile,header=None).values
L2_kx_max, L2_kx_min = np.max(L2[:,0]), np.min(L2[:,0])
L2_ky_max, L2_ky_min = np.max(L2[:,1]), np.min(L2[:,1])
L2_kz_max, L2_kz_min = np.max(L2[:,2]), np.min(L2[:,2])
L2_len = max(L2_kx_max - L2_kx_min, L2_ky_max - L2_ky_min, L2_kz_max - L2_kz_min)
L2_center_x = np.mean(L2[:,0])
L2_center_y = np.mean(L2[:,1])
L2_center_z = np.mean(L2[:,2])
L2_window = [[L2_center_x - L2_len/2, L2_center_x + L2_len/2],
             [L2_center_y - L2_len/2, L2_center_y + L2_len/2],
             [L2_center_z - L2_len/2, L2_center_z + L2_len/2]]

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title("L2")
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("kz")
ax.set_xlim(L2_window[0])
ax.set_ylim(L2_window[1])
ax.set_zlim(L2_window[2])
readfile = name+'_face.csv'
L2 = pd.read_csv(readfile,header=None)
L2 = L2.groupby((L2.isnull().all(axis=1)).cumsum())
for index, g in L2:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[7])

fn = index
str_vertex_num = "face num = "+str(fn)
ax_pos = ax.get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_vertex_num)

plt.savefig(png+"Fermi_surface_L2.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"Fermi_surface_L2.svg")
#plt.show()
plt.close()
# }}}

# L3 {{{
name = data+'triangle_L3-mu0e0'
readfile = name+'_vertex.csv'
L3 = pd.read_csv(readfile,header=None).values
L3_kx_max, L3_kx_min = np.max(L3[:,0]), np.min(L3[:,0])
L3_ky_max, L3_ky_min = np.max(L3[:,1]), np.min(L3[:,1])
L3_kz_max, L3_kz_min = np.max(L3[:,2]), np.min(L3[:,2])
L3_len = max(L3_kx_max - L3_kx_min, L3_ky_max - L3_ky_min, L3_kz_max - L3_kz_min)
L3_center_x = np.mean(L3[:,0])
L3_center_y = np.mean(L3[:,1])
L3_center_z = np.mean(L3[:,2])
L3_window = [[L3_center_x - L3_len/2, L3_center_x + L3_len/2],
             [L3_center_y - L3_len/2, L3_center_y + L3_len/2],
             [L3_center_z - L3_len/2, L3_center_z + L3_len/2]]

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, projection='3d')
ax.set_title("L3")
ax.set_xlabel("kx")
ax.set_ylabel("ky")
ax.set_zlabel("kz")
ax.set_xlim(L3_window[0])
ax.set_ylim(L3_window[1])
ax.set_zlim(L3_window[2])
readfile = name+'_face.csv'
L3 = pd.read_csv(readfile,header=None)
L3 = L3.groupby((L3.isnull().all(axis=1)).cumsum())
for index, g in L3:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax.plot(g[:,0],  g[:,1],  g[:,2], color=colors[6])

fn = index
str_vertex_num = "face num = "+str(fn)
ax_pos = ax.get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_vertex_num)

plt.savefig(png+"Fermi_surface_L3.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"Fermi_surface_L3.svg")
#plt.show()
plt.close()
# }}}

# Total {{{
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1, 2, 1, projection='3d')

ax1.set_xlabel("kx")
ax1.set_ylabel("ky")
ax1.set_zlabel("kz")
ax1.set_xlim(-1,1)
ax1.set_ylim(-1,1)
ax1.set_zlim(-1,1)
for index, g in T:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax1.plot(g[:,0],  g[:,1],  g[:,2], color=colors[9])
for index, g in L1:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax1.plot(g[:,0],  g[:,1],  g[:,2], color=colors[4])
for index, g in L2:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax1.plot(g[:,0],  g[:,1],  g[:,2], color=colors[7])
for index, g in L3:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        ax1.plot(g[:,0],  g[:,1],  g[:,2], color=colors[6])

axT  = fig.add_subplot(2,4,3, projection='3d')
axT.set_title("T")
axT.set_xlabel("kx")
axT.set_ylabel("ky")
axT.set_zlabel("kz")
axT.set_xlim(T_window[0])
axT.set_ylim(T_window[1])
axT.set_zlim(T_window[2])
for index, g in T:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        axT.plot(g[:,0],  g[:,1],  g[:,2], color=colors[9])

axL1 = fig.add_subplot(2,4,4, projection='3d')
axL1.set_title("L1")
axL1.set_xlabel("kx")
axL1.set_ylabel("ky")
axL1.set_zlabel("kz")
axL1.set_xlim(L1_window[0])
axL1.set_ylim(L1_window[1])
axL1.set_zlim(L1_window[2])
for index, g in L1:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        axL1.plot(g[:,0],  g[:,1],  g[:,2], color=colors[4])
axL2 = fig.add_subplot(2,4,7, projection='3d')
axL2.set_title("L2")
axL2.set_xlabel("kx")
axL2.set_ylabel("ky")
axL2.set_zlabel("kz")
axL2.set_xlim(L2_window[0])
axL2.set_ylim(L2_window[1])
axL2.set_zlim(L2_window[2])
for index, g in L2:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        axL2.plot(g[:,0],  g[:,1],  g[:,2], color=colors[7])
axL3 = fig.add_subplot(2,4,8, projection='3d')
axL3.set_title("L3")
axL3.set_xlabel("kx")
axL3.set_ylabel("ky")
axL3.set_zlabel("kz")
axL3.set_xlim(L3_window[0])
axL3.set_ylim(L3_window[1])
axL3.set_zlim(L3_window[2])
for index, g in L3:
    g = g.dropna()
    if len(g) > 0:
        g = g.values
        axL3.plot(g[:,0],  g[:,1],  g[:,2], color=colors[6])

plt.savefig(png+"Fermi_surface.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"Fermi_surface.svg")
#plt.show()
plt.close()
# }}}

# }}}
