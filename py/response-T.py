import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import pandas as pd
import configparser as cnf
from matplotlib.ticker import ScalarFormatter

home = "../"
data0 = "../dat/"
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

angstrom = 1e-10
hbar     = 6.582119569e-16
mass     = 9.10938356e-31
charge   = 1.60217662e-19
muB      = 9.2740100783e-24
muBeV    = 5.7883818060e-5
v0       = 1e6
g_ast    = 2e0*mass*v0**2/(delta*charge)

numeric = config['numeric']

markers = ["o", ",", "D", "v", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

epsilon = 1e-4
label = "{:.6f}".format(epsilon)
str_damping = "Gamma = "+"{:.2f}".format(epsilon*1e3)+" [meV]"

total_dos = []
total_conductivity = []
total_spin_conductivity1 = []
total_spin_conductivity2 = []

# T {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'

# dos T
# 1,1 layout -  {{{
fig, axes = plt.subplots(1,1,figsize=(7,7))
readfile = data+'dos_eps'+label+'.csv'
print(readfile)
dos = pd.read_csv(readfile,header=0).values
total_dos.append(dos)

axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
axes.set_title("dos")
axes.set_xlabel("mu [eV]")
axes.set_ylabel("DOS [/eV m3]")
axes.scatter(dos[:,0], dos[:,1], s=4, label="numeric")
axes.plot(dos[:,0], dos[:,1])
axes.legend()

plt.savefig(png+"dos_T-"+str(bandsT)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"dos_T-"+str(bandsT)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# T conductivity
# 1,2 layout - diagonal vs off-diagonal {{{
fig, axes = plt.subplots(1,2,figsize=(14,7))
axes = axes.flatten()
readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
conductivity = pd.read_csv(readfile,header=0).values
total_conductivity.append(conductivity)
window = [-0.01e0*conductivity[:,1:].max(), conductivity[:,1:].max()*1.1e0]
for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
axes[0].set_title("diagonal")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylabel("conductivity [/Ohm m]")
axes[0].set_ylim(window)
axes[0].scatter(conductivity[:,0], conductivity[:,1], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None', label="xx")
axes[0].scatter(conductivity[:,0], conductivity[:,5], s=3, marker=markers[1], edgecolors=colors[2], facecolor='None', label="yy")
axes[0].scatter(conductivity[:,0], conductivity[:,9], s=3, marker=markers[3], edgecolors=colors[1], facecolor='None', label="zz")
axes[0].plot(conductivity[:,0], conductivity[:,1], color=colors[3])
axes[0].plot(conductivity[:,0], conductivity[:,5], color=colors[2])
axes[0].plot(conductivity[:,0], conductivity[:,9], color=colors[1])
axes[0].legend()

axes[1].set_title("off-diagonal")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
axes[1].scatter(conductivity[:,0], conductivity[:,2], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None', label="xy")
axes[1].scatter(conductivity[:,0], conductivity[:,3], s=3, marker=markers[1], edgecolors=colors[2], facecolor='None', label="xz")
axes[1].scatter(conductivity[:,0], conductivity[:,4], s=3, marker=markers[3], edgecolors=colors[1], facecolor='None', label="yx")
axes[1].scatter(conductivity[:,0], conductivity[:,6], s=3, marker=markers[2], edgecolors=colors[4], facecolor='None', label="yz")
axes[1].scatter(conductivity[:,0], conductivity[:,7], s=3, marker=markers[4], edgecolors=colors[5], facecolor='None', label="zx")
axes[1].scatter(conductivity[:,0], conductivity[:,8], s=3, marker=markers[5], edgecolors=colors[6], facecolor='None', label="zy")
axes[1].plot(conductivity[:,0], conductivity[:,2], color=colors[3])
axes[1].plot(conductivity[:,0], conductivity[:,3], color=colors[2])
axes[1].plot(conductivity[:,0], conductivity[:,4], color=colors[1])
axes[1].plot(conductivity[:,0], conductivity[:,6], color=colors[4])
axes[1].plot(conductivity[:,0], conductivity[:,7], color=colors[5])
axes[1].plot(conductivity[:,0], conductivity[:,8], color=colors[6])
axes[1].legend()

ax_pos = axes[1].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.1, str_damping)

plt.savefig(png+"conductivity_T-"+str(bandsT)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity_T-"+str(bandsT)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

readfile = data+'mu-dependence/spin-conductivity1_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
conductivity1 = df.values
total_spin_conductivity1.append(conductivity1)

readfile = data+'mu-dependence/spin-conductivity2_eps'+label+'.csv'
conductivity2 = pd.read_csv(readfile,header=0).values
total_spin_conductivity2.append(conductivity2)

maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]

# T spin conductivity 1
# 1,3 layout -  {{{
fig, axes = plt.subplots(1,3,figsize=(21,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylim(window)
plot_lists = [6, 16, 20]
i = 0
for j in plot_lists:
    axes[0].scatter(conductivity1[:,0], conductivity1[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(conductivity1[:,0], conductivity1[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(conductivity1[:,0], conductivity1[:,j], s=3, label=titles[j])
axes[2].legend()

ax_pos = axes[2].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# T spin conductivity 2
# 1,3 layout -  {{{
fig, axes = plt.subplots(1,3,figsize=(21,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylim(window)
plot_lists = [6, 16, 20]
i = 0
for j in plot_lists:
    axes[0].scatter(conductivity2[:,0], conductivity2[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(conductivity2[:,0], conductivity2[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(conductivity2[:,0], conductivity2[:,j], s=1, label=titles[j])
axes[2].legend()

ax_pos = axes[2].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# T total spin Hall conductivity
# 1,1 layout -  {{{
conductivity3 = conductivity1[:,6] + conductivity2[:,6]
fig, axes = plt.subplots(1,1,figsize=(7,7))

axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes.set_title("compare")
axes.set_xlabel("mu [eV]")
axes.set_ylabel("SHC [/Ohm m]")
axes.scatter(conductivity1[:,0], conductivity1[:,6], s=2, label="SHC1:"+titles[6])
axes.scatter(conductivity2[:,0], conductivity2[:,6], s=2, label="SHC2:"+titles[6])
axes.scatter(conductivity1[:,0], conductivity3,      s=2, label="sum :"+titles[6])
axes.legend()

ax_pos = axes.get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.1, str_damping)

plt.savefig(png+"spin_conductivity_T-"+str(bandsT)+"bands_compare_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_conductivity_T-"+str(bandsT)+"bands_compare_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}
