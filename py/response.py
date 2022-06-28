import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import pandas as pd
import configparser as cnf
from matplotlib.ticker import ScalarFormatter
from scipy import interpolate

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

markers = ["o", ",", "D", "v", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

epsilon = 5e-4
label = "{:.6f}".format(epsilon)
str_damping = "Gamma = "+"{:.2f}".format(epsilon*1e3)+" [meV]"

labels = {"T,4":0, "L1,4":1, "L1,6":2, "L2,4":3, "L2,6":4, "L3,4":5, "L3,6":6}
valleys = {"T":0, "L1":1, "L2":2, "L3":3}
bands  = 12
bandsT = 12
bandsL = 4
lower_band_L = 0
upper_band_L = 2
#bandsL = 8
#lower_band_L = 2
#upper_band_L = lower_band_L + 2
#bandsL = 12
#lower_band_L = 4
#upper_band_L = 6
#bandsL = 16
#lower_band_L = 8
#upper_band_L = 10

e_max = []
e_min = []
dos = []
dos_i = []
dos_valley = []

conductivity = []
conductivity_i = []
conductivity_i_sum = []
conductivity_valley = []

conductivity1 = []
conductivity2 = []
conductivity1_i = []
conductivity2_i = []
conductivity1_i_sum = []
conductivity2_i_sum = []
conductivity1_valley = []
conductivity2_valley = []

axises = {"x":0, "y":1, "z":2}
e_cutoff = 0.08e0
k_cutoff = 0.27722e0

cutoff = 0.08
param = "cutoff"+str(cutoff)+"eV"

# dispersion T point {{{
fig, ax = plt.subplots(1,3,figsize=(15,7))
fig.tight_layout()
ax = ax.flatten()
ax[0].set_ylabel("Energy")
ax[1].set_title("T point")
ax[1].tick_params(labelleft=False)
ax[2].tick_params(labelleft=False)
for axis, val in axises.items():
    ax[val].set_xlabel("k"+axis)
    ax[val].set_ylim(-0.2,0.2)
    ax[val].axhspan(-e_cutoff, e_cutoff, color="gray", alpha=0.3)
    ax[val].axvspan(-k_cutoff, k_cutoff, color="tab:cyan", alpha=0.2)

    read_full = data0+'EigenValue_T_'+axis+'_full.dat'
    full = pd.read_csv(read_full,header=None).values
    read_effective = data0+'EigenValue_T_'+axis+'_'+str(bandsT)+'bands.dat'
    effective = pd.read_csv(read_effective,header=None).values

    ax[val].plot(full[:,0], full[:,9], color=colors[0], label="Liu-Allen")
    ax[val].plot(effective[:,0], effective[:,5], color=colors[7], label="effective")
    ax[val].legend()
#plt.show()

plt.subplots_adjust(wspace=0.05, hspace=0.35)

plt.savefig(png+"response_dispersion_T-"+str(bandsT)+"bands.png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"response_dispersion_T-"+str(bandsT)+"bands.svg")
plt.close()
# }}}

# dispersion L points {{{
for valley in np.arange(1,4):
    fig, ax = plt.subplots(1,3,figsize=(15,7))
    fig.tight_layout()
    ax = ax.flatten()
    ax[0].set_ylabel("Energy")
    ax[1].set_title("L"+str(valley)+" point")
    ax[1].tick_params(labelleft=False)
    ax[2].tick_params(labelleft=False)
    for axis, val in axises.items():
        ax[val].set_xlabel("k"+axis)
        ax[val].set_ylim(-0.2,0.2)
        ax[val].axhspan(-e_cutoff, e_cutoff, color="gray", alpha=0.3)
#        ax[val].axvspan(-k_cutoff, k_cutoff, color="tab:cyan", alpha=0.2)

        read_full = data0+'EigenValue_L'+str(valley)+'_'+axis+'_full.dat'
        full = pd.read_csv(read_full,header=None).values
        read_effective = data0+'EigenValue_L'+str(valley)+'_'+axis+'_'+str(bandsL)+'bands.dat'
        effective = pd.read_csv(read_effective,header=None).values

        ax[val].plot(full[:,0], full[:,9], color=colors[0], label="Liu-Allen")
        ax[val].plot(full[:,0], full[:,11], color=colors[0])
        ax[val].plot(effective[:,0], effective[:,lower_band_L+1], color=colors[7], label="effective")
        ax[val].plot(effective[:,0], effective[:,upper_band_L+1], color=colors[7])
        ax[val].legend()
    #plt.show()

    plt.savefig(png+"response_dispersion_L"+str(valley)+"-"+str(bandsL)+"bands.png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"response_dispersion_L"+str(valley)+"-"+str(bandsL)+"bands.svg")
    plt.close()
# }}}

# dos {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
readfile = data+'dos_eps'+label+'.csv'
d = pd.read_csv(readfile,header=None).values
dos.append(d)
e_max.append(d[:,0].max())
e_min.append(d[:,0].min())
dos_valley.append(d)

for valley in np.arange(1,4):
    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'
    readfile = data+'dos_eps'+label+'.csv'
    d4 = pd.read_csv(readfile,header=None).values
    dos.append(d4)
    e_max.append(d4[:,0].max())
    e_min.append(d4[:,0].min())

    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'
    readfile = data+'dos_eps'+label+'.csv'
    d6 = pd.read_csv(readfile,header=None).values
    dos.append(d6)
    e_max.append(d6[:,0].max())
    e_min.append(d6[:,0].min())

    dos_valley.append(np.append(d4, d6, axis=0))

e_max = sorted(set(e_max))
e_min = sorted(set(e_min))
print(e_max, e_min)
n = 16000.0
e_range = np.linspace(-8e-2*n, 8e-2*n, int(n)-1)
x = e_range / n
for d in dos_valley:
    ip1d = interpolate.interp1d(d[:,0], d[:,1])
    dos_i.append(ip1d(x))
dos_i = np.array(dos_i)
dos_i_sum = np.sum(dos_i, axis=0)

fig, ax = plt.subplots(1,1,figsize=(7,7))
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
ax.set_title("dos")
ax.set_xlabel("mu [eV]")
ax.set_ylabel("DOS [/eV m3]")
ax.plot(x, dos_i_sum, label="total", c="red")
for d in dos_i:
    ax.plot(x, d, c="black", lw=1)
for key, val in labels.items():
    ax.scatter(dos[val][:,0], dos[val][:,1], s=4, label=key)
ax.legend()
#plt.show()

plt.savefig(png+"dos_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"dos_total_"+str(bands)+"bands_gamma"+label+".svg")
plt.close()
# }}}

# conductivity {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
d = df.values
conductivity.append(d)
conductivity_valley.append(d)

for valley in np.arange(1,4):
    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'
    readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
    d4 = pd.read_csv(readfile,header=0).values
    conductivity.append(d4)

    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'
    readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
    d6 = pd.read_csv(readfile,header=0).values
    conductivity.append(d6)

    conductivity_valley.append(np.append(d4, d6, axis=0))

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

plot_lists = [1, 5, 9]
for s in np.arange(0,3):
    conductivity_tmp = []
    axes[s].set_title(titles[plot_lists[s]])
    for d in conductivity_valley:
        ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
        conductivity_tmp.append(ip1d(x))
    conductivity_tmp = np.array(conductivity_tmp)
    conductivity_i.append(conductivity_tmp)
    conductivity_tmp = np.sum(conductivity_tmp, axis=0)
    conductivity_i_sum.append(conductivity_tmp)

maximun = 0e0
for s in np.arange(0,3):
    maximun = max([maximun, np.abs(conductivity_i_sum[s].max()), np.abs(conductivity_i_sum[s].min())])

window = [-maximun*0.01e0, maximun*1.1e0]

axes[0].set_ylabel("conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].plot(x, conductivity_i_sum[s], label="total", c="red")
    for d in conductivity_i[s]:
        axes[s].plot(x, d, c="black", lw=1)
    for key, val in labels.items():
        axes[s].scatter(conductivity[val][:,0], conductivity[val][:,plot_lists[s]], s=4, label=key)
    axes[s].legend()
#plt.show()

plt.savefig(png+"conductivity_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity_total_"+str(bands)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin magnetic conductivity1 {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
d = df.values
conductivity1.append(d)
conductivity1_valley.append(d)

for valley in np.arange(1,4):
    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'
    readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
    d4 = pd.read_csv(readfile,header=0).values
    conductivity1.append(d4)

    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'
    readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
    d6 = pd.read_csv(readfile,header=0).values
    conductivity1.append(d6)

    conductivity1_valley.append(np.append(d4, d6, axis=0))

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

plot_lists = [6, 16, 20]
for s in np.arange(0,3):
    conductivity1_tmp = []
    axes[s].set_title(titles[plot_lists[s]])
    for d in conductivity1_valley:
        ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
        conductivity1_tmp.append(ip1d(x))
    conductivity1_tmp = np.array(conductivity1_tmp)
    conductivity1_i.append(conductivity1_tmp)
    conductivity1_tmp = np.sum(conductivity1_tmp, axis=0)
    conductivity1_i_sum.append(conductivity1_tmp)

maximun = 0e0
for s in np.arange(0,3):
    maximun = max([maximun, np.abs(conductivity1_i_sum[s].max()), np.abs(conductivity1_i_sum[s].min())])

window = [-maximun*1.1e0, maximun*1.1e0]

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].plot(x, conductivity1_i_sum[s], label="total", c="red")
    for d in conductivity1_i[s]:
        axes[s].plot(x, d, c="black", lw=1)
    for key, val in labels.items():
        axes[s].scatter(conductivity1[val][:,0], conductivity1[val][:,plot_lists[s]], s=4, label=key)
    axes[s].legend()
#plt.show()

plt.savefig(png+"spin-magnetic-conductivity1_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin-magnetic-conductivity1_total_"+str(bands)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin magnetic conductivity2 {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
d = pd.read_csv(readfile,header=0).values
conductivity2.append(d)
conductivity2_valley.append(d)

for valley in np.arange(1,4):
    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'
    readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
    d4 = pd.read_csv(readfile,header=0).values
    conductivity2.append(d4)

    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'
    readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
    d6 = pd.read_csv(readfile,header=0).values

    val = np.copy(d4[-1, :])
    val[0] = 0e0
    d6 = d6 + val

    conductivity2.append(d6)
    conductivity2_valley.append(np.append(d4, d6, axis=0))

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

maximun = 0e0
plot_lists = [6, 16, 20]
for s in np.arange(0,3):
    conductivity2_tmp = []
    axes[s].set_title(titles[plot_lists[s]])
    for d in conductivity2_valley:
        ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
        conductivity2_tmp.append(ip1d(x))
    conductivity2_tmp = np.array(conductivity2_tmp)
    conductivity2_i.append(conductivity2_tmp)
    conductivity2_tmp = np.sum(conductivity2_tmp, axis=0)
    conductivity2_i_sum.append(conductivity2_tmp)

    maximun = max(maximun, np.abs(conductivity2_tmp.min()))

window = [-maximun*1.1e0, maximun*1.1e0]

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].plot(x, conductivity2_i_sum[s], label="total", c="red")
    for d in conductivity2_i[s]:
        axes[s].plot(x, d, c="black", lw=1)
    for key, val in labels.items():
        axes[s].scatter(conductivity2[val][:,0], conductivity2[val][:,plot_lists[s]], s=4, label=key)
    axes[s].legend()

plt.savefig(png+"spin-magnetic-conductivity2_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin-magnetic-conductivity2_total_"+str(bands)+"bands_gamma"+label+".svg")
plt.close()
# }}}

# spin magnetic conductivity 1 + 2 {{{
conductivity3_i_sum = list(range(3))
plot_lists = [6, 16, 20]
for s in np.arange(0,3):
    conductivity3_i_sum[s] = conductivity1_i_sum[s] + conductivity2_i_sum[s]

maximun = 0e0
for s in np.arange(0,3):
    maximun = max([maximun, np.abs(conductivity3_i_sum[s].max()), np.abs(conductivity3_i_sum[s].min())])

window = [-maximun*1.1e0, maximun*1.1e0]

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].set_title(titles[plot_lists[s]])
    axes[s].plot(x, conductivity1_i_sum[s], label="SHC1",  c=colors[4])
    axes[s].plot(x, conductivity2_i_sum[s], label="SHC2",  c=colors[7])
    axes[s].plot(x, conductivity3_i_sum[s], label="total", c=colors[6])
    axes[s].legend()

plt.savefig(png+"spin-magnetic-conductivity3_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin-magnetic-conductivity3_total_"+str(bands)+"bands_gamma"+label+".svg")
plt.close()
# }}}

conductivity1 = []
conductivity2 = []
conductivity1_i = []
conductivity2_i = []
conductivity1_i_sum = []
conductivity2_i_sum = []
conductivity1_valley = []
conductivity2_valley = []

# spin angular conductivity1 {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
d = df.values
conductivity1.append(d)
conductivity1_valley.append(d)

for valley in np.arange(1,4):
    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'
    readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    d4 = pd.read_csv(readfile,header=0).values
    conductivity1.append(d4)

    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'
    readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    d6 = pd.read_csv(readfile,header=0).values
    conductivity1.append(d6)

    conductivity1_valley.append(np.append(d4, d6, axis=0))

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

plot_lists = [6, 16, 20]
for s in np.arange(0,3):
    conductivity1_tmp = []
    axes[s].set_title(titles[plot_lists[s]])
    for d in conductivity1_valley:
        ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
        conductivity1_tmp.append(ip1d(x))
    conductivity1_tmp = np.array(conductivity1_tmp)
    conductivity1_i.append(conductivity1_tmp)
    conductivity1_tmp = np.sum(conductivity1_tmp, axis=0)
    conductivity1_i_sum.append(conductivity1_tmp)

maximun = 0e0
for s in np.arange(0,3):
    maximun = max([maximun, np.abs(conductivity1_i_sum[s].max()), np.abs(conductivity1_i_sum[s].min())])

window = [-maximun*1.1e0, maximun*1.1e0]

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].plot(x, conductivity1_i_sum[s], label="total", c="red")
    for d in conductivity1_i[s]:
        axes[s].plot(x, d, c="black", lw=1)
    for key, val in labels.items():
        axes[s].scatter(conductivity1[val][:,0], conductivity1[val][:,plot_lists[s]], s=4, label=key)
    axes[s].legend()
#plt.show()

plt.savefig(png+"spin-angular-conductivity1_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin-angular-conductivity1_total_"+str(bands)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin angular conductivity2 {{{
data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
d = pd.read_csv(readfile,header=0).values
conductivity2.append(d)
conductivity2_valley.append(d)

for valley in np.arange(1,4):
    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'
    readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    d4 = pd.read_csv(readfile,header=0).values
    conductivity2.append(d4)

    data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'
    readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    d6 = pd.read_csv(readfile,header=0).values

    val = np.copy(d4[-1, :])
    val[0] = 0e0
    d6 = d6 + val

    conductivity2.append(d6)
    conductivity2_valley.append(np.append(d4, d6, axis=0))

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

maximun = 0e0
plot_lists = [6, 16, 20]
for s in np.arange(0,3):
    conductivity2_tmp = []
    axes[s].set_title(titles[plot_lists[s]])
    for d in conductivity2_valley:
        ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
        conductivity2_tmp.append(ip1d(x))
    conductivity2_tmp = np.array(conductivity2_tmp)
    conductivity2_i.append(conductivity2_tmp)
    conductivity2_tmp = np.sum(conductivity2_tmp, axis=0)
    conductivity2_i_sum.append(conductivity2_tmp)

    maximun = max(maximun, np.abs(conductivity2_tmp.min()))

window = [-maximun*1.1e0, maximun*1.1e0]

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].plot(x, conductivity2_i_sum[s], label="total", c="red")
    for d in conductivity2_i[s]:
        axes[s].plot(x, d, c="black", lw=1)
    for key, val in labels.items():
        axes[s].scatter(conductivity2[val][:,0], conductivity2[val][:,plot_lists[s]], s=4, label=key)
    axes[s].legend()

plt.savefig(png+"spin-angular-conductivity2_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin-angular-conductivity2_total_"+str(bands)+"bands_gamma"+label+".svg")
plt.close()
# }}}

# spin angular conductivity 1 + 2 {{{
conductivity3_i_sum = list(range(3))
plot_lists = [6, 16, 20]
for s in np.arange(0,3):
    conductivity3_i_sum[s] = conductivity1_i_sum[s] + conductivity2_i_sum[s]

maximun = 0e0
for s in np.arange(0,3):
    maximun = max([maximun, np.abs(conductivity3_i_sum[s].max()), np.abs(conductivity3_i_sum[s].min())])

window = [-maximun*1.1e0, maximun*1.1e0]

fig, axes = plt.subplots(1,3,figsize=(15,7))
axes = axes.flatten()

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for ax in axes:
    ax.set_xlabel("mu [eV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlim(-0.08, 0.08)
    ax.set_ylim(window)

for s in np.arange(0,3):
    axes[s].set_title(titles[plot_lists[s]])
    axes[s].plot(x, conductivity1_i_sum[s], label="SHC1",  c=colors[4])
    axes[s].plot(x, conductivity2_i_sum[s], label="SHC2",  c=colors[7])
    axes[s].plot(x, conductivity3_i_sum[s], label="total", c=colors[6])
    axes[s].legend()

plt.savefig(png+"spin-angular-conductivity3_total_"+str(bands)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin-angular-conductivity3_total_"+str(bands)+"bands_gamma"+label+".svg")
plt.close()
# }}}
