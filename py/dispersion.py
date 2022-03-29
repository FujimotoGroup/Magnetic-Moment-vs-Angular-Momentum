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

# T point {{{
fig = plt.figure(figsize=(20,10))
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)
ax1.set_xlabel("kx")
ax2.set_xlabel("kx")
ax3.set_xlabel("kx")
ax1.set_ylabel("Energy")
ax1.set_ylim(-13,3)
ax2.set_ylim(-13,3)
ax1.set_title("16-full bands")
ax2.set_title(str(bandsT)+"-effective bands")
ax3.set_title("comparison")

read_full = data+'EigenValue_T_x_full.dat'
full = pd.read_csv(read_full,header=None).values
for i in np.arange(1,bands,2):
    ax1.scatter(full[:,0], full[:,i], color=colors[i-1])

read_effective = data+'EigenValue_T_'+str(bandsT)+'bands.dat'
effective = pd.read_csv(read_effective,header=None).values
for i in np.arange(1,bandsT,2):
    ax2.plot(effective[:,0], effective[:,i], color=colors[i+lowest_T-2])

for i in np.arange(lowest_T,lowest_T+bandsT,2):
    ax3.scatter(full[:,0], full[:,i], label="band#"+str(i))
for i in np.arange(1,bandsT,2):
    ax3.plot(effective[:,0], effective[:,i], color=colors[i+lowest_T-2], label="effective band#"+str(i))
ax3.legend()

plt.savefig(png+"dispersion_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"dispersion_T-"+str(bandsT)+"bands.svg")
#plt.show()
plt.close()
# }}}

# L points {{{
for valley in np.arange(1,4):
    fig = plt.figure(figsize=(20,10))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    ax1.set_xlabel("kx")
    ax2.set_xlabel("kx")
    ax3.set_xlabel("kx")
    ax1.set_ylabel("Energy")
    ax1.set_ylim(-13,3)
    ax2.set_ylim(-13,3)
    ax1.set_title("16-full bands")
    ax2.set_title(str(bandsL)+"-effective bands")
    ax3.set_title("comparison")

    read_full = data+"EigenValue_L"+str(valley)+"_x_full.dat"
    full = pd.read_csv(read_full,header=None).values
    for i in np.arange(1,16,2):
        ax1.scatter(full[:,0], full[:,i], color=colors[i-1])

    read_effective = data+"EigenValue_L"+str(valley)+"_"+str(bandsL)+"bands.dat"
    effective = pd.read_csv(read_effective,header=None).values
    for i in np.arange(1,4,2):
        ax2.plot(effective[:,0], effective[:,i], color=colors[i+lowest_L-2])

    for i in np.arange(lowest_L,lowest_L+bandsL,2):
        ax3.scatter(full[:,0], full[:,i], label="band#"+str(i))
    for i in np.arange(1,bandsL,2):
        ax3.plot(effective[:,0], effective[:,i], color=colors[i+lowest_L-2], label="effective band#"+str(i))
    ax3.legend()

    plt.savefig(png+"dispersion_L"+str(valley)+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"dispersion_L"+str(valley)+".svg")
#    plt.show()
    plt.close()
# }}}

