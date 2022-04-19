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
# 2,3 layout - bare {{{
fig, axes = plt.subplots(2,3,figsize=(12,10))
axes = axes.flatten()
for i in np.arange(6):
    label = "{:.6f}".format(1e-1**(i+1))
    readfile = data+'T'+str(bandsT)+'bands/conductivity_eps'+label+'.csv'
    conductivity = pd.read_csv(readfile,header=0).values
    axes[i].set_title("eps = "+label)
    axes[i].set_xlabel("mu")
    axes[i].plot(conductivity[:,0], conductivity[:,1], color=colors[3], label="xx")
    axes[i].plot(conductivity[:,0], conductivity[:,5], color=colors[2], label="yy")
    axes[i].plot(conductivity[:,0], conductivity[:,9], color=colors[1], label="zz")
    axes[i].legend()

#plt.savefig(png+"conductivity_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"conductivity_T-"+str(bandsT)+"bands.svg")
#plt.show()
plt.close()
# }}}
# 1,3 layout - log {{{
title = ["xx", "yy", "zz"]
index = [1, 5, 9]
fig, axes = plt.subplots(1,3,figsize=(12,10))
axes = axes.flatten()
for i in np.arange(3):
    axes[i].set_title("sigma_"+title[i])
    axes[i].set_xlabel("mu")
    axes[i].set_ylim(-5,25)
    for j in np.arange(6):
        label = "{:.6f}".format(1e-1**(j+1))
        readfile = data+'T'+str(bandsT)+'bands/conductivity_eps'+label+'.csv'
        conductivity = pd.read_csv(readfile,header=0).values
        axes[i].plot(conductivity[:,0], np.log(conductivity[:,index[i]]), color=colors[j], label=label)
    axes[i].legend()

plt.savefig(png+"conductivity_log_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity_log_T-"+str(bandsT)+"bands.svg")
#plt.show()
plt.close()
# }}}
# }}}
