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

# T effective dispersion {{{
fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(1,3,1)
ax2 = fig.add_subplot(1,3,2)
ax3 = fig.add_subplot(1,3,3)
ax1.set_xlabel("kx")
ax2.set_xlabel("ky")
ax3.set_xlabel("kz")
ax1.set_ylabel("Energy")
ax1.set_ylim(-2.0,1.5)
ax2.set_ylim(-2.0,1.5)
ax3.set_ylim(-2.0,1.5)
ax1.set_title("kx-cut")
ax2.set_title("ky-cut")
ax3.set_title("kz-cut")

read_effective_x = data+'EigenValue_T_x_'+str(bandsT)+'bands.dat'
effective_x = pd.read_csv(read_effective_x,header=None).values
read_effective_y = data+'EigenValue_T_y_'+str(bandsT)+'bands.dat'
effective_y = pd.read_csv(read_effective_y,header=None).values
read_effective_z = data+'EigenValue_T_z_'+str(bandsT)+'bands.dat'
effective_z = pd.read_csv(read_effective_z,header=None).values
for i in np.arange(1,bandsT,2):
    ax1.plot(effective_x[:,0], effective_x[:,i], color=colors[i+lowest_T-2])
    ax2.plot(effective_y[:,0], effective_y[:,i], color=colors[i+lowest_T-2])
    ax3.plot(effective_z[:,0], effective_z[:,i], color=colors[i+lowest_T-2])

plt.savefig(png+"dispersion_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"dispersion_T-"+str(bandsT)+"bands.svg")
#plt.show()
plt.close()
# }}}

# L effective dispersion {{{
for valley in np.arange(1,4):
    fig = plt.figure(figsize=(12,5))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    ax1.set_xlabel("kx")
    ax2.set_xlabel("ky")
    ax3.set_xlabel("kz")
    ax1.set_ylabel("Energy")
    ax1.set_ylim(-2.0,1.5)
    ax2.set_ylim(-2.0,1.5)
    ax3.set_ylim(-2.0,1.5)
    ax1.set_title("kx-cut")
    ax2.set_title("ky-cut")
    ax3.set_title("kz-cut")

    read_effective_x = data+"EigenValue_L"+str(valley)+"_x_"+str(bandsL)+"bands.dat"
    effective_x = pd.read_csv(read_effective_x,header=None).values
    read_effective_y = data+"EigenValue_L"+str(valley)+"_y_"+str(bandsL)+"bands.dat"
    effective_y = pd.read_csv(read_effective_y,header=None).values
    read_effective_z = data+"EigenValue_L"+str(valley)+"_z_"+str(bandsL)+"bands.dat"
    effective_z = pd.read_csv(read_effective_z,header=None).values
    for i in np.arange(1,bandsL,2):
        ax1.plot(effective_x[:,0], effective_x[:,i], color=colors[i+lowest_L-2])
        ax2.plot(effective_y[:,0], effective_y[:,i], color=colors[i+lowest_L-2])
        ax3.plot(effective_z[:,0], effective_z[:,i], color=colors[i+lowest_L-2])

    plt.savefig(png+"dispersion_L"+str(valley)+"-"+str(bandsL)+"bands.png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"dispersion_L"+str(valley)+"-"+str(bandsL)+"bands.svg")
    #plt.show()
    plt.close()
# }}}
