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

EL0 = float(physics.get('EL0'))
EL2 = float(physics.get('EL2'))
delta = (EL2 - EL0)/2e0

numeric = config['numeric']

colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

# dos {{{
# 1,1 layout -  {{{
def dos(x):
    if abs(x - delta) > 0e0:
        return np.abs(x)*np.sqrt(x**2 - delta**2)
    else:
        return 0e0

fig, axes = plt.subplots(1,1,figsize=(10,10))
i = 3
label = "{:.6f}".format(1e-1**(i+1))
readfile = data+'L1_'+str(bandsL)+'bands/dos_eps'+label+'.csv'
dos = pd.read_csv(readfile,header=0).values

n = 1000.0
mu_range = np.linspace(dos[:,0].min()*n, dos[:,0].max()*n, 1000)
mu_range = mu_range / n

axes.set_title("dos")
axes.set_xlabel("mu")
axes.set_ylim(-0.00001,0.001)
axes.scatter(dos[:,0], dos[:,1], s=2, label="numeric")

#plt.savefig(png+"conductivity_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"conductivity_T-"+str(bandsT)+"bands.svg")
plt.show()
plt.close()
# }}}
# }}}

# L conductivity {{{
# 1,2 layout - diagonal vs off-diagonal {{{
fig, axes = plt.subplots(1,2,figsize=(10,5))
axes = axes.flatten()
i = 3
label = "{:.6f}".format(1e-1**(i+1))
readfile = data+'L1_'+str(bandsL)+'bands/conductivity_real_eps'+label+'.csv'
conductivity = pd.read_csv(readfile,header=0).values
axes[0].set_title("diagonal")
axes[0].set_xlabel("mu")
axes[0].set_ylabel("conductivity")
axes[0].set_ylim(-0.2,17.5)
axes[0].scatter(conductivity[:,0], conductivity[:,1], color=colors[3], s=1, label="xx")
axes[0].scatter(conductivity[:,0], conductivity[:,5], color=colors[2], s=1, label="yy")
axes[0].scatter(conductivity[:,0], conductivity[:,9], color=colors[1], s=1, label="zz")
axes[0].legend()

axes[1].set_title("off-diagonal")
axes[1].set_xlabel("mu")
axes[1].set_ylim(-0.2,17.5)
axes[1].scatter(conductivity[:,0], conductivity[:,2], color=colors[3], s=1, label="xy")
axes[1].scatter(conductivity[:,0], conductivity[:,3], color=colors[2], s=1, label="xz")
axes[1].scatter(conductivity[:,0], conductivity[:,4], color=colors[1], s=1, label="yx")
axes[1].scatter(conductivity[:,0], conductivity[:,6], color=colors[4], s=1, label="yz")
axes[1].scatter(conductivity[:,0], conductivity[:,7], color=colors[5], s=1, label="zx")
axes[1].scatter(conductivity[:,0], conductivity[:,8], color=colors[6], s=1, label="zy")
axes[1].legend()

axes[1].text(-0.02, 19, "eps = "+label+" [eV]")

plt.savefig(png+"conductivity_L1-"+str(bandsL)+"bands.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity_L1-"+str(bandsL)+"bands.svg")
#plt.show()
plt.close()
# }}}
# }}}
# L spin conductivity {{{
# 1,3 layout -  {{{
#fig, axes = plt.subplots(1,3,figsize=(10,5))
fig, axes = plt.subplots(1,3)
axes = axes.flatten()
i = 3
label = "{:.6f}".format(1e-1**(i+1))
readfile = data+'L1_'+str(bandsL)+'bands/spin_Hall_conductivity1_real_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
conductivity = df.values
window = [-0.015,0.015]

axes[0].set_ylabel("spin conductivity")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu")
axes[0].set_ylim(window)
plot_lists = [6, 16, 20]
for j in plot_lists:
    axes[0].scatter(conductivity[:,0], conductivity[:,j], s=1, label=titles[j])
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(conductivity[:,0], conductivity[:,j], s=1, label=titles[j])
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(conductivity[:,0], conductivity[:,j], s=1, label=titles[j])
axes[2].legend()

axes[2].text(-0.02, 19, "eps = "+label+" [eV]")

#plt.savefig(png+"spin_conductivity_L1-"+str(bandsL)+"bands.png", bbox_inches = 'tight')
plt.savefig(png+"spin_conductivity_L1-"+str(bandsL)+"bands.png")
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_conductivity_L1-"+str(bandsL)+"bands.svg")
#plt.show()
plt.close()
# }}}
# 1,1 layout -  {{{
fig, axes = plt.subplots(1,1,figsize=(10,10))
i = 3
label = "{:.6f}".format(1e-1**(i+1))
readfile = data+'L1_'+str(bandsL)+'bands/spin_Hall_conductivity1_real_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
conductivity = df.values

n = 1000.0
mu_range = np.linspace(conductivity[:,0].min()*n, conductivity[:,0].max()*n, 1000)
mu_range = mu_range / n

axes.set_title("fitting")
axes.set_xlabel("mu")
axes.set_ylim(-0.001,0.015)
plot_lists = [6, 16, 20]
for j in plot_lists:
    axes.scatter(conductivity[:,0], conductivity[:,j], s=2, label=titles[j])
axes.legend()

axes.text(-0.02, 19, "eps = "+label+" [eV]")

#plt.savefig(png+"conductivity_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"conductivity_T-"+str(bandsT)+"bands.svg")
plt.show()
plt.close()
# }}}
# }}}
