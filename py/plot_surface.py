import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from scipy.spatial import Delaunay

for i in np.arange(5):
    readfile = '../dat/L0_band0k'+str(i)+'.csv'
    readfile = '../dat/T_band4k'+str(i)+'.csv'
    position = pd.read_csv(readfile, header=None).values
    x = position[:,0]
    y = position[:,1]
    z = position[:,2]

    plt.rcParams["axes.facecolor"] = 'white'
    fig = plt.figure(figsize=plt.figaspect(0.5))

    ax = fig.add_subplot(1, 2, 1, projection='3d')
    ax.scatter(x, y, z, label="e = "+str(position[0,6]))
    ax.legend()
    plt.show()
