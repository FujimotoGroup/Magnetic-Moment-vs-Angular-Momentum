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

# L {{{
for valley in np.arange(1,4):
    data1 = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index0/'
    data2 = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index2/'

# dos L {{{
# 1,1 layout -  {{{
    fig, axes = plt.subplots(1,1,figsize=(7,7))
    readfile = data1+'dos_eps'+label+'.csv'
    dos = pd.read_csv(readfile,header=0).values
    readfile = data2+'dos_eps'+label+'.csv'
    dos = np.append(dos, pd.read_csv(readfile,header=0).values, axis=0)
    total_dos.append(dos)

    axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    axes.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    axes.set_title("dos")
    axes.set_xlabel("mu [eV]")
    axes.set_ylabel("DOS [/eV m3]")
    axes.scatter(dos[:,0], dos[:,1], s=4, label="numeric")
    axes.plot(dos[:,0], dos[:,1])
    axes.legend()

    plt.savefig(png+"dos_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"dos_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}

# L conductivity {{{
# 1,2 layout - diagonal vs off-diagonal {{{
    fig, axes = plt.subplots(1,2,figsize=(14,7))
    axes = axes.flatten()
    readfile = data1+'mu-dependence/conductivity_eps'+label+'.csv'
    conductivity = pd.read_csv(readfile,header=0).values
    readfile = data2+'mu-dependence/conductivity_eps'+label+'.csv'
    conductivity = np.append(conductivity, pd.read_csv(readfile,header=0).values, axis=0)
    total_conductivity.append(conductivity)
    window = [-1.1e0*conductivity[:,1:].max(), conductivity[:,1:].max()*1.1e0]
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

    plt.savefig(png+"conductivity_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"conductivity_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}

    readfile = data1+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    conductivity1 = df.values
    readfile = data2+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    conductivity1 = np.append(conductivity1, df.values, axis=0)
    total_spin_magnetic_conductivity1.append(conductivity1)

    readfile = data1+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
    conductivity2 = pd.read_csv(readfile,header=0).values
    readfile = data2+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
    val = np.copy(conductivity2[-1, :])
    val[0] = 0e0
    conductivity2_tmp = pd.read_csv(readfile,header=0).values + val
    conductivity2 = np.append(conductivity2, conductivity2_tmp, axis=0)
    total_spin_magnetic_conductivity2.append(conductivity2)

    maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
    window = [-maximun, maximun]

# L spin magnetic conductivity 1 {{{
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

    plt.savefig(png+"spin_magnetic_conductivity1_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"spin_magnetic_conductivity1_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}
# L spin magnetic conductivity 2 {{{
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

    plt.savefig(png+"spin_magnetic_conductivity2_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"spin_magnetic_conductivity2_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}
# L spin magnetic Hall conductivity {{{
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

    plt.savefig(png+"spin_magnetic_conductivity_L"+str(valley)+"-"+str(bandsL)+"bands_compare_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"spin_magnetic_conductivity_L"+str(valley)+"-"+str(bandsL)+"bands_compare_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}

    readfile = data1+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    conductivity1 = df.values
    readfile = data2+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    conductivity1 = np.append(conductivity1, df.values, axis=0)
    total_spin_angular_conductivity1.append(conductivity1)

    readfile = data1+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    conductivity2 = pd.read_csv(readfile,header=0).values
    readfile = data2+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    val = np.copy(conductivity2[-1, :])
    val[0] = 0e0
    conductivity2_tmp = pd.read_csv(readfile,header=0).values + val
    conductivity2 = np.append(conductivity2, conductivity2_tmp, axis=0)
    total_spin_angular_conductivity2.append(conductivity2)

    maximun = max(np.abs(conductivity1[:,1:].max()*1.1e0), np.abs(conductivity2[:,1:].max()*1.1e0))
    window = [-maximun, maximun]

# L spin angular conductivity 1 {{{
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

    plt.savefig(png+"spin_angular_conductivity1_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"spin_angular_conductivity1_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}
# L spin angular conductivity 2 {{{
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

    plt.savefig(png+"spin_angular_conductivity2_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"spin_angular_conductivity2_L"+str(valley)+"-"+str(bandsL)+"bands_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}
# L spin angular Hall conductivity {{{
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

    plt.savefig(png+"spin_angular_conductivity_L"+str(valley)+"-"+str(bandsL)+"bands_compare_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
    plt.rc("svg", fonttype="none")
    plt.savefig(svg+"spin_angular_conductivity_L"+str(valley)+"-"+str(bandsL)+"bands_compare_gamma"+label+".svg")
    #plt.show()
    plt.close()
# }}}
# }}}

dos = total_dos[0] + total_dos[1] + total_dos[2]
conductivity = total_conductivity[0] + total_conductivity[1] + total_conductivity[2]
magnetic_conductivity1 = total_spin_magnetic_conductivity1[0] + total_spin_magnetic_conductivity1[1] + total_spin_magnetic_conductivity1[2]
magnetic_conductivity2 = total_spin_magnetic_conductivity2[0] + total_spin_magnetic_conductivity2[1] + total_spin_magnetic_conductivity2[2]
angular_conductivity1 = total_spin_angular_conductivity1[0] + total_spin_angular_conductivity1[1] + total_spin_angular_conductivity1[2]
angular_conductivity2 = total_spin_angular_conductivity2[0] + total_spin_angular_conductivity2[1] + total_spin_angular_conductivity2[2]
dos[:,0] = total_dos[0][:,0]
conductivity[:,0] = total_conductivity[0][:,0]
magnetic_conductivity1[:,0] = total_spin_magnetic_conductivity1[0][:,0]
magnetic_conductivity2[:,0] = total_spin_magnetic_conductivity2[0][:,0]
angular_conductivity1[:,0] = total_spin_angular_conductivity1[0][:,0]
angular_conductivity2[:,0] = total_spin_angular_conductivity2[0][:,0]

#dos
# 1,1 layout -  {{{
fig, axes = plt.subplots(1,1,figsize=(7,7))

axes.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
axes.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
axes.set_title("dos")
axes.set_xlabel("mu [eV]")
axes.set_ylabel("DOS [/eV m3]")
axes.scatter(dos[:,0], dos[:,1], s=4, label="numeric")
axes.plot(dos[:,0], dos[:,1])
axes.legend()

plt.savefig(png+"dos_L-total-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"dos_L-total-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# conductivity
# 1,2 layout - diagonal vs off-diagonal {{{
fig, axes = plt.subplots(1,2,figsize=(14,7))
axes = axes.flatten()
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

plt.savefig(png+"conductivity_L-total-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity_L-total-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

maximun = max(np.abs(magnetic_conductivity1[:,1:].max()*1.1e0), np.abs(magnetic_conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]
# spin magnetic conductivity1
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
    axes[0].scatter(magnetic_conductivity1[:,0], magnetic_conductivity1[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(magnetic_conductivity1[:,0], magnetic_conductivity1[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(magnetic_conductivity1[:,0], magnetic_conductivity1[:,j], s=3, label=titles[j])
axes[2].legend()

ax_pos = axes[2].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_magnetic_conductivity1_L-total-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity1_L-total-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin magnetic conductivity2
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
    axes[0].scatter(magnetic_conductivity2[:,0], magnetic_conductivity2[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(magnetic_conductivity2[:,0], magnetic_conductivity2[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(magnetic_conductivity2[:,0], magnetic_conductivity2[:,j], s=1, label=titles[j])
axes[2].legend()

ax_pos = axes[2].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_magnetic_conductivity2_L-total-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity2_L-total-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

maximun = max(np.abs(angular_conductivity1[:,1:].max()*1.1e0), np.abs(angular_conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]
# spin angular conductivity1
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
    axes[0].scatter(angular_conductivity1[:,0], angular_conductivity1[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(angular_conductivity1[:,0], angular_conductivity1[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(angular_conductivity1[:,0], angular_conductivity1[:,j], s=3, label=titles[j])
axes[2].legend()

ax_pos = axes[2].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_angular_conductivity1_L-total-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity1_L-total-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin angular conductivity2
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
    axes[0].scatter(angular_conductivity2[:,0], angular_conductivity2[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[0].legend()

axes[1].set_title("off-diagonal 2")
axes[1].set_xlabel("mu [eV]")
axes[1].set_ylim(window)
plot_lists = [8, 12, 22]
for j in plot_lists:
    axes[1].scatter(angular_conductivity2[:,0], angular_conductivity2[:,j], s=3, marker=markers[i], edgecolors=colors[i], facecolor='None', label=titles[j])
    i = i + 1
axes[1].legend()

axes[2].set_title("others")
axes[2].set_xlabel("mu [eV]")
axes[2].set_ylim(window)
plot_lists = [1, 2, 3, 4, 5, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 23, 24, 25, 26, 27]
for j in plot_lists:
    axes[2].scatter(angular_conductivity2[:,0], angular_conductivity2[:,j], s=1, label=titles[j])
axes[2].legend()

ax_pos = axes[2].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_angular_conductivity2_L-total-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity2_L-total-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# comparison
# conductivity 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()
window = [-0.01e0*conductivity[:,1:].max(), conductivity[:,1:].max()*1.1e0]

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("conductivity [/Ohm m]")
for i in np.arange(3):
    axes[i].set_title("L"+str(i+1))
    axes[i].set_xlabel("mu [eV]")
    axes[i].set_ylim(window)
    axes[i].scatter(total_conductivity[i][:,0], total_conductivity[i][:,1], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None', label="xx")
    axes[i].scatter(total_conductivity[i][:,0], total_conductivity[i][:,5], s=3, marker=markers[1], edgecolors=colors[2], facecolor='None', label="yy")
    axes[i].scatter(total_conductivity[i][:,0], total_conductivity[i][:,9], s=3, marker=markers[3], edgecolors=colors[1], facecolor='None', label="zz")
    axes[i].plot(total_conductivity[i][:,0], total_conductivity[i][:,1], color=colors[3])
    axes[i].plot(total_conductivity[i][:,0], total_conductivity[i][:,5], color=colors[2])
    axes[i].plot(total_conductivity[i][:,0], total_conductivity[i][:,9], color=colors[1])
    axes[i].legend()

axes[3].set_title("sum")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
axes[3].scatter(conductivity[:,0], conductivity[:,1], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None', label="xx")
axes[3].scatter(conductivity[:,0], conductivity[:,5], s=3, marker=markers[1], edgecolors=colors[2], facecolor='None', label="yy")
axes[3].scatter(conductivity[:,0], conductivity[:,9], s=3, marker=markers[3], edgecolors=colors[1], facecolor='None', label="zz")
axes[3].plot(conductivity[:,0], conductivity[:,1], color=colors[3])
axes[3].plot(conductivity[:,0], conductivity[:,5], color=colors[2])
axes[3].plot(conductivity[:,0], conductivity[:,9], color=colors[1])
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"conductivity1_L-comparison-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"conductivity1_L-comparison-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

maximun = max(np.abs(magnetic_conductivity1[:,1:].max()*1.1e0), np.abs(magnetic_conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]
# spin magnetic conductivity1 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for i in np.arange(3):
    axes[i].set_title("L"+str(i+1))
    axes[i].set_xlabel("mu [eV]")
    axes[i].set_ylim(window)
    k = 0
    plot_lists = [6, 16, 20]
    for j in plot_lists:
        axes[i].scatter(total_spin_magnetic_conductivity1[i][:,0], total_spin_magnetic_conductivity1[i][:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
        axes[i].plot(total_spin_magnetic_conductivity1[i][:,0], total_spin_magnetic_conductivity1[i][:,j], color=colors[k])
        k = k + 1
    axes[i].legend()

axes[3].set_title("sum")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
k = 0
plot_lists = [6, 16, 20]
for j in plot_lists:
    axes[3].scatter(magnetic_conductivity1[:,0], magnetic_conductivity1[:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
    axes[3].plot(magnetic_conductivity1[:,0], magnetic_conductivity1[:,j], color=colors[k])
    k = k + 1
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_magnetic_conductivity1_L-comparison-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity1_L-comparison-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin magnetic conductivity2 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for i in np.arange(3):
    axes[i].set_title("L"+str(i+1))
    axes[i].set_xlabel("mu [eV]")
    axes[i].set_ylim(window)
    k = 0
    plot_lists = [6, 16, 20]
    for j in plot_lists:
        axes[i].scatter(total_spin_magnetic_conductivity2[i][:,0], total_spin_magnetic_conductivity2[i][:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
        axes[i].plot(total_spin_magnetic_conductivity2[i][:,0], total_spin_magnetic_conductivity2[i][:,j], color=colors[k])
        k = k + 1
    axes[i].legend()

axes[3].set_title("sum")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
k = 0
plot_lists = [6, 16, 20]
for j in plot_lists:
    axes[3].scatter(magnetic_conductivity2[:,0], magnetic_conductivity2[:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
    axes[3].plot(magnetic_conductivity2[:,0], magnetic_conductivity2[:,j], color=colors[k])
    k = k + 1
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_magnetic_conductivity2_L-comparison-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_magnetic_conductivity2_L-comparison-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

maximun = max(np.abs(angular_conductivity1[:,1:].max()*1.1e0), np.abs(angular_conductivity2[:,1:].max()*1.1e0))
window = [-maximun, maximun]
# spin angular conductivity1 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for i in np.arange(3):
    axes[i].set_title("L"+str(i+1))
    axes[i].set_xlabel("mu [eV]")
    axes[i].set_ylim(window)
    k = 0
    plot_lists = [6, 16, 20]
    for j in plot_lists:
        axes[i].scatter(total_spin_angular_conductivity1[i][:,0], total_spin_angular_conductivity1[i][:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
        axes[i].plot(total_spin_angular_conductivity1[i][:,0], total_spin_angular_conductivity1[i][:,j], color=colors[k])
        k = k + 1
    axes[i].legend()

axes[3].set_title("sum")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
k = 0
plot_lists = [6, 16, 20]
for j in plot_lists:
    axes[3].scatter(angular_conductivity1[:,0], angular_conductivity1[:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
    axes[3].plot(angular_conductivity1[:,0], angular_conductivity1[:,j], color=colors[k])
    k = k + 1
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_angular_conductivity1_L-comparison-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity1_L-comparison-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}

# spin angular conductivity2 1,4 layout -  {{{
fig, axes = plt.subplots(1,4,figsize=(28,7))
axes = axes.flatten()

for ax in axes:
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))

axes[0].set_ylabel("spin conductivity [/Ohm m]")
for i in np.arange(3):
    axes[i].set_title("L"+str(i+1))
    axes[i].set_xlabel("mu [eV]")
    axes[i].set_ylim(window)
    k = 0
    plot_lists = [6, 16, 20]
    for j in plot_lists:
        axes[i].scatter(total_spin_angular_conductivity2[i][:,0], total_spin_angular_conductivity2[i][:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
        axes[i].plot(total_spin_angular_conductivity2[i][:,0], total_spin_angular_conductivity2[i][:,j], color=colors[k])
        k = k + 1
    axes[i].legend()

axes[3].set_title("sum")
axes[3].set_xlabel("mu [eV]")
axes[3].set_ylim(window)
k = 0
plot_lists = [6, 16, 20]
for j in plot_lists:
    axes[3].scatter(conductivity2[:,0], conductivity2[:,j], s=3, marker=markers[k], edgecolors=colors[k], facecolor='None', label=titles[j])
    axes[3].plot(conductivity2[:,0], conductivity2[:,j], color=colors[k])
    k = k + 1
axes[3].legend()

ax_pos = axes[3].get_position()
fig.text(ax_pos.x1 - 0.15, ax_pos.y1 + 0.05, str_damping)

plt.savefig(png+"spin_angular_conductivity2_L-comparison-"+str(bandsL)+"bands_gamma"+label+".png", bbox_inches = 'tight', dpi=300)
plt.rc("svg", fonttype="none")
plt.savefig(svg+"spin_angular_conductivity2_L-comparison-"+str(bandsL)+"bands_gamma"+label+".svg")
#plt.show()
plt.close()
# }}}
# }}}
