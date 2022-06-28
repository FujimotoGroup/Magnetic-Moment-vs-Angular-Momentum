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

import re

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

cutoff = 0.09
param = "cutoff"+str(cutoff)+"eV"

markers = ["o", ",", "D", "v", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo', '#990099']

n = 1000.0
mu_range1 = np.linspace((-delta-5e-2)*n, -delta*n, int(n))
x1 = mu_range1 / n
n = 100.0
mu_range2 = np.linspace(-delta*n, 0e0, int(n))
x2 = mu_range2 / n
x = np.append(x1, x2)

epsilon = 5e-4
label = "{:.6f}".format(epsilon)
str_damping = "Gamma = "+"{:.2f}".format(epsilon*1e3)+" [meV]"

total_dos = []
total_conductivity = []
total_spin_magnetic_conductivity1 = []
total_spin_magnetic_conductivity2 = []
total_spin_angular_conductivity1 = []
total_spin_angular_conductivity2 = []

for valley in np.arange(1,4):
    data1 = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index0/'+param+'/'
    data2 = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index2/'+param+'/'

    readfile = data1+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    conductivity1 = df.values
    readfile = data2+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    conductivity1 = np.append(conductivity1, df.values, axis=0)

    readfile = data1+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    conductivity2 = pd.read_csv(readfile,header=0).values
    readfile = data2+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    val = np.copy(conductivity2[-1, :])
    val[0] = 0e0
    conductivity2_tmp = pd.read_csv(readfile,header=0).values + val
    conductivity2 = np.append(conductivity2, conductivity2_tmp, axis=0)

    axis = ["x", "y", "z"]
    i = 1
    name = {}
    for external in axis:
        for response in axis:
            for spin in axis:
                name[i] = external+response+spin
                i = i + 1

    def R(theta):
        return np.matrix([[ np.cos(theta), -np.sin(theta), 0e0],
                         [ np.sin(theta),  np.cos(theta), 0e0],
                         [ 0e0          ,            0e0, 1e0]])

    for a in axis:
        component = {}
        tmp = []
        for key, val in name.items():
            if re.match('.*'+a+'$', val):
                component[key] = val
                tmp.append(key)
        matrix1 = conductivity1[:,tmp]
        matrix2 = conductivity2[:,tmp]

        matrix1a = matrix1.reshape(matrix1.shape[0],3,3)
        matrix2a = matrix2.reshape(matrix2.shape[0],3,3)
        matrix1b = np.empty(matrix1a.shape)
        matrix2b = np.empty(matrix2a.shape)
        matrix1c = np.empty(matrix1a.shape)
        matrix2c = np.empty(matrix2a.shape)

        theta = 2e0*np.pi/3e0
        for m in np.arange(matrix1.shape[0]):
            matrix1b[m] = np.matmul(np.matmul(R(theta)[:,:], matrix1a[m,:,:]), R(-theta)[:,:])
            matrix2b[m] = np.matmul(np.matmul(R(theta)[:,:], matrix2a[m,:,:]), R(-theta)[:,:])

        theta =-2e0*np.pi/3e0
        for m in np.arange(matrix1.shape[0]):
            matrix1c[m] = np.matmul(np.matmul(R(theta)[:,:], matrix1a[m,:,:]), R(-theta)[:,:])
            matrix2c[m] = np.matmul(np.matmul(R(theta)[:,:], matrix2a[m,:,:]), R(-theta)[:,:])

        matrix1 = matrix1a + matrix1b + matrix1c
        matrix2 = matrix2a + matrix2b + matrix2c

        matrix1 = matrix1.reshape(matrix1.shape[0],9)
        matrix2 = matrix2.reshape(matrix2.shape[0],9)

        for j in np.arange(9):
            conductivity1[:,tmp[j]] = matrix1[:,j]
            conductivity2[:,tmp[j]] = matrix2[:,j]

    total_spin_angular_conductivity1.append(conductivity1)
    total_spin_angular_conductivity2.append(conductivity2)

conductivity1 = (total_spin_angular_conductivity1[0] + total_spin_angular_conductivity1[1] + total_spin_angular_conductivity1[2])/3e0
conductivity2 = (total_spin_angular_conductivity2[0] + total_spin_angular_conductivity2[1] + total_spin_angular_conductivity2[2])/3e0
conductivity1[:,0] = total_spin_angular_conductivity1[0][:,0]
conductivity2[:,0] = total_spin_angular_conductivity2[0][:,0]

maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]

plot_dict1 = {1:6, 2:16, 3:20}
plot_dict1_inset = {2:16, 3:20}
plot_dict2 = {1:12, 2:8, 3:22}
plot_dict2_inset = {2:8, 3:22}
plot_dict3 = {1:1, 2:13, 3:5, 4:11}
plot_dict4a = {0:2, 1:3, 2:4, 3:7, 4:9, 5:10, 6:14, 7:15, 8:17}
plot_dict4b = {0:18, 1:19, 2:21, 3:23, 4:24, 5:25, 6:26, 7:27}

plot_dict = {1:12, 2:8, 3:20, 4:13}

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

left, bottom, width, height = [0.2, 0.1, 0.35, 0.35]
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

left, bottom, width, height = [0.2, 0.05, 0.3, 0.3]

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

#plt.savefig(png+"spin_angular_conductivity1_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"spin_angular_conductivity1_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+"_"+param+".svg")
plt.show()
plt.close()
# }}}

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

left, bottom, width, height = [0.2, 0.6, 0.35, 0.35]
ax0 = axes[0].inset_axes([left, bottom, width, height])
ax0.set_xlim([-cutoff, -cutoff+0.01])
ax0.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax0.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
for key, j in plot_dict1_inset.items():
    ax0.scatter(conductivity2[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[0].indicate_inset_zoom(ax0, edgecolor="black")

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
for key, j in plot_dict2.items():
    axes[1].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[1].legend()

axes[2].set_title("xxx, etc.")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
for key, j in plot_dict3.items():
    axes[2].scatter(conductivity1[:,0], conductivity2[:,j], s=7, marker=markers[key], edgecolors=colors[key], facecolor='None', label=titles[j])
axes[2].legend()

left, bottom, width, height = [0.4, 0.1, 0.3, 0.3]
ax2 = axes[2].inset_axes([left, bottom, width, height])
ax2.set_xlim([cutoff-0.01, cutoff])
mask = conductivity1[:,0] > cutoff-0.01
masked = conductivity2[mask,:]
maximun2 = np.abs(masked[:,list(plot_dict3.values())].max()*1.1e0)
ax2.set_ylim([-maximun2, maximun2])
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

#plt.savefig(png+"spin_angular_conductivity2_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"spin_angular_conductivity2_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+"_"+param+".svg")
plt.show()
plt.close()
# }}}


total_spin_angular_conductivity1 = []
total_spin_angular_conductivity2 = []

data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'

readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
df = pd.read_csv(readfile,header=0)
titles = df.columns.values
conductivity1 = df.values

readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
conductivity2 = pd.read_csv(readfile,header=0).values

maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]

fig, axes = plt.subplots(1,3,figsize=(21,7))
axes = axes.flatten()
i = 0
for a in axis:
    component = {}
    tmp = []
    for key, val in name.items():
        if re.match('.*'+a+'$', val):
            component[key] = val
            tmp.append(key)
    matrix1 = conductivity1[:,tmp]
    matrix2 = conductivity2[:,tmp]

    matrix1a = matrix1.reshape(matrix1.shape[0],3,3)
    matrix2a = matrix2.reshape(matrix2.shape[0],3,3)
    matrix1b = np.empty(matrix1a.shape)
    matrix2b = np.empty(matrix2a.shape)
    matrix1c = np.empty(matrix1a.shape)
    matrix2c = np.empty(matrix2a.shape)

    theta = 2e0*np.pi/3e0
    for m in np.arange(matrix1.shape[0]):
        matrix1b[m] = np.matmul(np.matmul(R(theta)[:,:], matrix1a[m,:,:]), R(-theta)[:,:])
        matrix2b[m] = np.matmul(np.matmul(R(theta)[:,:], matrix2a[m,:,:]), R(-theta)[:,:])

    theta =-2e0*np.pi/3e0
    for m in np.arange(matrix1.shape[0]):
        matrix1c[m] = np.matmul(np.matmul(R(theta)[:,:], matrix1a[m,:,:]), R(-theta)[:,:])
        matrix2c[m] = np.matmul(np.matmul(R(theta)[:,:], matrix2a[m,:,:]), R(-theta)[:,:])

    matrix1 = (matrix1a + matrix1b + matrix1c)/3e0
    matrix2 = (matrix2a + matrix2b + matrix2c)/3e0

    matrix1 = matrix1.reshape(matrix1.shape[0],9)
    matrix2 = matrix2.reshape(matrix2.shape[0],9)

    for j in np.arange(9):
        conductivity1[:,tmp[j]] = matrix1[:,j]
        conductivity2[:,tmp[j]] = matrix2[:,j]

        axes[i].set_ylim(window)
        axes[i].scatter(conductivity1[:,0], matrix2[:,j], s=7, marker=markers[j], edgecolors=colors[j], facecolor='None', label=titles[tmp[j]])
    axes[i].legend()

    i = i + 1
plt.show()

total_spin_angular_conductivity1.append(conductivity1)
total_spin_angular_conductivity2.append(conductivity2)

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

#plt.savefig(png+"spin_angular_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"spin_angular_conductivity1_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
plt.show()
plt.close()
# }}}

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

#plt.savefig(png+"spin_angular_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".png", bbox_inches = 'tight', dpi=300)
#plt.rc("svg", fonttype="none")
#plt.savefig(svg+"spin_angular_conductivity2_T-"+str(bandsT)+"bands_gamma"+label+"_"+param+".svg")
plt.show()
plt.close()
# }}}
