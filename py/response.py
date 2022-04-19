import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
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

# T conductivity {{{

conductivity = np.array([]);

fig, axes = plt.subplots(2,3,figsize=(12,5))
axes = axes.flatten()
for i in np.arange(6):
    label = "{:.6f}".format(1e-1**(i+1))
    readfile = data+'T'+str(bandsT)+'bands/conductivity_eps'+label+'.csv'
    conductivity = pd.read_csv(readfile,header=0).values
    axes[i].set_title("eps = "+label)
    axes[i].set_xlabel("mu")
    axes[i].plot(conductivity[:,0], conductivity[:,1], color=colors[0])
    axes[i].plot(conductivity[:,0], conductivity[:,5], color=colors[1])
    axes[i].plot(conductivity[:,0], conductivity[:,9], color=colors[2])

#plt.savefig(png+"dispersion_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"dispersion_T-"+str(bandsT)+"bands.svg")
##plt.show()
plt.close()
# }}}
