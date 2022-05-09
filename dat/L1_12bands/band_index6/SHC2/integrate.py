import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import configparser as cnf
import glob

home = "./"
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

markers = ["o", ",", "D", "v", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

SHC_files1 = sorted(glob.glob("./spin_conductivity2_mu-*"),reverse=True)
SHC_files2 = sorted(glob.glob("./spin_conductivity2_mu0*"))
SHC_files = SHC_files1 + SHC_files2
for file in SHC_files:
    print(file)
    SHC = pd.read_csv(file,header=0).values
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_xlim(-1e-3,1e-3)
    ax.set_ylim(-20,20)
    ax.scatter(SHC[:,0], SHC[:,6], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None')
    ax.plot(SHC[:,0], SHC[:,6])
    plt.show()
    plt.close()
