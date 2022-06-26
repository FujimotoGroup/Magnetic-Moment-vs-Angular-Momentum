import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import pandas as pd
import configparser as cnf
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

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
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo', '#990099']

epsilon = 5e-4
label = "{:.6f}".format(epsilon)
str_damping = "Gamma = "+"{:.2f}".format(epsilon*1e3)+" [meV]"

total_dos = []
total_conductivity = []
total_spin_conductivity1 = []
total_spin_conductivity2 = []

cutoff = 0.1
param = "cutoff"+str(cutoff)+"eV"
# T {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'

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

plt.savefig(png+"dos_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"dos_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
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

plt.savefig(png+"conductivity_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}

plot_dict1 = {1:6, 2:16, 3:20}
plot_dict1_inset = {2:16, 3:20}
plot_dict2 = {1:12, 2:8, 3:22}
plot_dict2_inset = {2:8, 3:22}
plot_dict3 = {1:1, 2:13, 3:5, 4:11}
plot_dict4a = {0:2, 1:3, 2:4, 3:7, 4:9, 5:10, 6:14, 7:15, 8:17}
plot_dict4b = {0:18, 1:19, 2:21, 3:23, 4:24, 5:25, 6:26, 7:27}

plot_dict = {1:12, 2:8, 3:20, 4:13}

##############################################################################

readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
conductivity1 = df.values
total_spin_conductivity1.append(conductivity1)

readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
conductivity2 = pd.read_csv(readfile,header=0).values
total_spin_conductivity2.append(conductivity2)

maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]

# T spin magnetic conductivity 1
# 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylim(window)
for key, j in plot_dict1.items():
    axes[0].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].legend()

left, bottom, width, height = [0.2, 0.6, 0.35, 0.35]
ax0 = axes[0].inset_axes([left, bottom, width, height])
ax0.set_xlim([-cutoff, -cutoff+0.01])
ax0.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax0.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict1_inset.items():
    ax0.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].indicate_inset_zoom(ax0, edgecolor="black")

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
for key, j in plot_dict2.items():
    axes[1].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].legend()

left, bottom, width, height = [0.2, 0.1, 0.3, 0.3]
ax1 = axes[1].inset_axes([left, bottom, width, height])
ax1.set_xlim([-cutoff, -cutoff+0.01])
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict2_inset.items():
    ax1.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].indicate_inset_zoom(ax1, edgecolor="black")

axes[2].set_title("xxx, etc.")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
for key, j in plot_dict3.items():
    axes[2].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].legend()

left, bottom, width, height = [0.2, 0.1, 0.3, 0.3]
ax2 = axes[2].inset_axes([left, bottom, width, height])
ax2.set_xlim([-cutoff, -cutoff+0.01])
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict3.items():
    ax2.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].indicate_inset_zoom(ax2, edgecolor="black")

axes[3].set_title("others")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
for key, j in plot_dict4a.items():
    axes[3].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
for key, j in plot_dict4b.items():
    axes[3].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], c=colors[key+len(plot_dict4a)], label=titles[j])
axes[3].legend()

left, bottom, width, height = [0.2, 0.1, 0.3, 0.3]
ax3 = axes[3].inset_axes([left, bottom, width, height])
ax3.set_xlim([-cutoff, -cutoff+0.01])
ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax3.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict4a.items():
    ax3.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
for key, j in plot_dict4b.items():
    axes[3].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], c=colors[key+len(plot_dict4a)], label=titles[j])
axes[3].indicate_inset_zoom(ax3, edgecolor="black")

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_magnetic_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}

# T spin magnetic conductivity 2
# 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylim(window)
for key, j in plot_dict1.items():
    axes[0].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].legend()

left, bottom, width, height = [0.4, 0.6, 0.35, 0.35]
ax0 = axes[0].inset_axes([left, bottom, width, height])
ax0.set_xlim([cutoff-0.01, cutoff])
ax0.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax0.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict1_inset.items():
    ax0.scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].indicate_inset_zoom(ax0, edgecolor="black")

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
for key, j in plot_dict2.items():
    axes[1].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].legend()

left, bottom, width, height = [0.4, 0.1, 0.3, 0.3]
ax1 = axes[1].inset_axes([left, bottom, width, height])
ax1.set_xlim([cutoff-0.01, cutoff])
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict2_inset.items():
    ax1.scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].indicate_inset_zoom(ax1, edgecolor="black")

axes[2].set_title("xxx, etc.")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
for key, j in plot_dict3.items():
    axes[2].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].legend()

left, bottom, width, height = [0.4, 0.1, 0.3, 0.3]
ax2 = axes[2].inset_axes([left, bottom, width, height])
ax2.set_xlim([cutoff-0.01, cutoff])
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict3.items():
    ax2.scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].indicate_inset_zoom(ax2, edgecolor="black")

axes[3].set_title("others")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
for key, j in plot_dict4a.items():
    axes[3].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
for key, j in plot_dict4b.items():
    axes[3].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], c=colors[key+len(plot_dict4a)], label=titles[j])
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_magnetic_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}

# T total spin magnetic Hall conductivity
# 1,4 layout -  {{{
conductivity3 = conductivity1 + conductivity2
conductivity3[:,0] = conductivity1[:,0]

maximun = np.abs(conductivity3[:,1:].max()*1.1e0)
window = [-maximun*0.01, maximun]

fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

axes[0].set_ylabel("SHC [/Ohm m]")
for j in np.arange(4):
    axes[j].set_ylim(window)
    axes[j].xaxis.set_major_locator(MultipleLocator(0.05))
    axes[j].xaxis.set_major_formatter(FormatStrFormatter('%0.2f'))
    axes[j].xaxis.set_minor_locator(MultipleLocator(0.01))
    axes[j].yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes[j].ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    axes[j].set_xlabel("mu [eV]")

    axes[j].scatter(conductivity1[:,0], conductivity1[:,plot_dict[j+1]], s=7, marker=markers[0], label="SHC1:"+titles[plot_dict[j+1]])
    axes[j].scatter(conductivity2[:,0], conductivity2[:,plot_dict[j+1]], s=7, marker=markers[1], label="SHC2:"+titles[plot_dict[j+1]])
    axes[j].scatter(conductivity3[:,0], conductivity3[:,plot_dict[j+1]], s=7, marker=markers[2], label="sum :"+titles[plot_dict[j+1]])
    axes[j].legend()

for j in np.arange(1,4):
    left, bottom, width, height = [0.1, 0.3, 0.8, 0.35]
    ax = axes[j].inset_axes([left, bottom, width, height])
    ax.set_xticks(np.arange(-2,3)*0.05)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.scatter(conductivity1[:,0], conductivity1[:,plot_dict[j+1]], s=7, marker=markers[0], label="SHC1:"+titles[plot_dict[j+1]])
    ax.scatter(conductivity2[:,0], conductivity2[:,plot_dict[j+1]], s=7, marker=markers[1], label="SHC2:"+titles[plot_dict[j+1]])
    ax.scatter(conductivity3[:,0], conductivity3[:,plot_dict[j+1]], s=7, marker=markers[2], label="sum :"+titles[plot_dict[j+1]])
    axes[j].indicate_inset_zoom(ax, edgecolor="black")

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.1, ax_pos.y1 + 0.01, str_damping)

plt.savefig(png+"spin_magnetic_conductivity_T-"+str(bandsT)+"bands_compare_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity_T-"+str(bandsT)+"bands_compare_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}

total_spin_conductivity1 = []
total_spin_conductivity2 = []

readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
conductivity1 = df.values
total_spin_conductivity1.append(conductivity1)

readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
conductivity2 = pd.read_csv(readfile,header=0).values
total_spin_conductivity2.append(conductivity2)

maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]

# T spin angular conductivity 1
# 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylim(window)
for key, j in plot_dict1.items():
    axes[0].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].legend()

left, bottom, width, height = [0.2, 0.6, 0.35, 0.35]
ax0 = axes[0].inset_axes([left, bottom, width, height])
ax0.set_xlim([-cutoff, -cutoff+0.01])
ax0.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax0.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict1_inset.items():
    ax0.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].indicate_inset_zoom(ax0, edgecolor="black")

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
for key, j in plot_dict2.items():
    axes[1].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].legend()

left, bottom, width, height = [0.2, 0.1, 0.3, 0.3]
ax1 = axes[1].inset_axes([left, bottom, width, height])
ax1.set_xlim([-cutoff, -cutoff+0.01])
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict2_inset.items():
    ax1.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].indicate_inset_zoom(ax1, edgecolor="black")

axes[2].set_title("yxx, etc.")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
for key, j in plot_dict3.items():
    axes[2].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].legend()

left, bottom, width, height = [0.2, 0.1, 0.3, 0.3]
ax2 = axes[2].inset_axes([left, bottom, width, height])
ax2.set_xlim([-cutoff, -cutoff+0.01])
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict3.items():
    ax2.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].indicate_inset_zoom(ax2, edgecolor="black")

axes[3].set_title("others")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
for key, j in plot_dict4a.items():
    axes[3].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
for key, j in plot_dict4b.items():
    axes[3].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], c=colors[key+len(plot_dict4a)], label=titles[j])
axes[3].legend()

left, bottom, width, height = [0.2, 0.1, 0.3, 0.3]
ax3 = axes[3].inset_axes([left, bottom, width, height])
ax3.set_xlim([-cutoff, -cutoff+0.01])
ax3.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax3.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict4a.items():
    ax3.scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
for key, j in plot_dict4b.items():
    axes[3].scatter(conductivity1[:,0], conductivity1[:,j], s=7, marker=markers[key], c=colors[key+len(plot_dict4a)], label=titles[j])
axes[3].indicate_inset_zoom(ax3, edgecolor="black")

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_angular_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}

# T spin angular conductivity 2
# 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
axes[0].set_title("off-diagonal 1")
axes[0].set_xlabel("mu [eV]")
axes[0].set_ylim(window)
for key, j in plot_dict1.items():
    axes[0].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].legend()

left, bottom, width, height = [0.4, 0.6, 0.35, 0.35]
ax0 = axes[0].inset_axes([left, bottom, width, height])
ax0.set_xlim([cutoff-0.01, cutoff])
ax0.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax0.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict1_inset.items():
    ax0.scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].indicate_inset_zoom(ax0, edgecolor="black")

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
for key, j in plot_dict2.items():
    axes[1].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].legend()

left, bottom, width, height = [0.4, 0.1, 0.3, 0.3]
ax1 = axes[1].inset_axes([left, bottom, width, height])
ax1.set_xlim([cutoff-0.01, cutoff])
ax1.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax1.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict2_inset.items():
    ax1.scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].indicate_inset_zoom(ax1, edgecolor="black")

axes[2].set_title("yxx, etc.")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
for key, j in plot_dict3.items():
    axes[2].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].legend()

left, bottom, width, height = [0.4, 0.1, 0.3, 0.3]
ax2 = axes[2].inset_axes([left, bottom, width, height])
ax2.set_xlim([cutoff-0.01, cutoff])
ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax2.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict3.items():
    ax2.scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].indicate_inset_zoom(ax2, edgecolor="black")

axes[3].set_title("others")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
for key, j in plot_dict4a.items():
    axes[3].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
for key, j in plot_dict4b.items():
    axes[3].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], c=colors[key+len(plot_dict4a)], label=titles[j])
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_angular_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}

# T total spin angular Hall conductivity
# 1,4 layout -  {{{
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
fig.text(ax_pos.x1 - 0.05, ax_pos.y1 + 0.1, str_damping)

plt.savefig(png+"spin_angular_conductivity_T-"+str(bandsT)+"bands_compare_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity_T-"+str(bandsT)+"bands_compare_gamma"+label+"_"+param+".svg")
#plt.show()
plt.close()
# }}}
