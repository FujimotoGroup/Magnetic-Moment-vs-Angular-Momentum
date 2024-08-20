import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import pandas as pd
import configparser as cnf
from matplotlib.ticker import ScalarFormatter
from scipy import interpolate

matplotlib.rcParams["font.family"] = 'Times New Roman'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams["font.size"] = 16
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

home = "../"
data0 = "../dat/"
pdf = "../fig/pdf/"

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

labels = {"T,4":0, "L1,4":1, "L1,6":2, "L2,4":3, "L2,6":4, "L3,4":5, "L3,6":6}
valleys = ["T", "L1", "L2", "L3", "Total", "L sum"]
bands  = 12
bandsT = 12
bandsL = 4
lower_band_L = 0
upper_band_L = 2

axises = {"x":0, "y":1, "z":2}

def set_conductivity(nbands): # {{{
    conductivity = []

    data = data0+'T'+str(bandsT)+'bands/damping-dependence'
    readfile = data+'/conductivity_T.csv'
    df = pd.read_csv(readfile,header=0)
    d = df.values
    conductivity.append(d)

    for valley in np.arange(1,4):
        readfile = "../dat/L"+str(nbands)+"bands/damping-dependence/conductivity_L"+str(valley)+".csv"
        d = pd.read_csv(readfile,header=0).values
        conductivity.append(d)

    conductivity.append(conductivity[0] + conductivity[1] + conductivity[2] + conductivity[3])
    conductivity[4][:,0] = conductivity[0][:,0]

    conductivity.append(conductivity[1] + conductivity[2] + conductivity[3])
    conductivity[5][:,0] = conductivity[0][:,0]

    return conductivity
# }}}
def set_conductivity_constant_gamma(nbands): # {{{
    conductivity = []

    data = data0+'T'+str(bandsT)+'bands/damping-dependence'
    readfile = data+'/conductivity_constant-gamma_T.csv'
    df = pd.read_csv(readfile,header=0)
    d = df.values
    conductivity.append(d)

    for valley in np.arange(1,4):
        readfile = "../dat/L"+str(nbands)+"bands/damping-dependence/conductivity_constant-gamma_L"+str(valley)+".csv"
        d = pd.read_csv(readfile,header=0).values
        conductivity.append(d)

    conductivity.append(conductivity[0] + conductivity[1] + conductivity[2] + conductivity[3])
    conductivity[4][:,0] = conductivity[0][:,0]

    conductivity.append(conductivity[1] + conductivity[2] + conductivity[3])
    conductivity[5][:,0] = conductivity[0][:,0]

    return conductivity
# }}}

def plot_diagonal(conductivity, label): # {{{
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    maximun = 0e0
    minimun = 1e15
    for s in np.arange(0,5):
        for i in plot_list_conductivity_diagonal.keys():
            maximun = max([maximun, np.abs(conductivity[s][:,i].max()), np.abs(conductivity[s][:,i].min())])
            minimun = min([minimun, np.abs(conductivity[s][:,i].min())])
    window = [minimun*0.9e0, maximun*1.1e0]

    fig, axes = plt.subplots(1,5,figsize=(20,7))
    axes = axes.flatten()

    axes[0].set_ylabel("conductivity [/Ohm m]")
    for ax in axes:
        ax.set_xlabel("damping constant [meV]")
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.set_yscale('log')
        ax.set_ylim(window)

    for s in np.arange(0,5):
        axes[s].set_title(valleys[s])
        i = 0
        for key, val in plot_list_conductivity_diagonal.items():
            axes[s].scatter(conductivity[s][:,0]*1e3, conductivity[s][:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
            i = i + 1
        axes[s].legend()
    #plt.show()
    plt.savefig(pdf+"conductivity_damping_dependence_"+label+"_diagonal.pdf")
    plt.close()

    plot_list_conductivity_offdiagonal = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}
    maximun = 0e0
    minimun = 1e9
    for s in np.arange(0,5):
        for key, val in plot_list_conductivity_offdiagonal.items():
            maximun = max([maximun, np.abs(conductivity[s][:,key].max()), np.abs(conductivity[s][:,i].min())])
    window = [-maximun*1.1e0, maximun*1.1e0]
# }}}
def plot_diagonal_total(conductivity, label): #  {{{
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_ylabel("conductivity [/Ohm m]")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_title(valleys[-2])
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(conductivity[-1][:,0]*1e3, conductivity[-1][:,key], edgecolors=colors[i], facecolor='None', label=val)
        i = i + 1
    ax.legend()
    #plt.show()
    plt.savefig(pdf+"conductivity_damping_dependence_"+label+"_diagonal_total.pdf")
    plt.close()
# }}}

def plot_off_diagonal(conductivity, label): # {{{
    plot_list_conductivity_offdiagonal = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}

    maximun = 0e0
    for s in np.arange(0,5):
        for i in plot_list_conductivity_offdiagonal.keys():
            maximun = max([maximun, np.abs(conductivity[s][:,i].max()), np.abs(conductivity[s][:,i].min())])
    window = [-maximun*1.1e0, maximun*1.1e0]

    fig, axes = plt.subplots(1,5,figsize=(20,7))
    axes[0].set_ylabel("conductivity [/Ohm m]")
    for ax in axes:
        ax.set_xlabel("damping constant [meV]")
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    #    ax.set_yscale('log')
    #    ax.set_xlim(-0.08, 0.08)
        ax.set_ylim(window)

    for s in np.arange(0,5):
        axes[s].set_title(valleys[s])
        for key, val in plot_list_conductivity_offdiagonal.items():
            axes[s].scatter(conductivity[s][:,0]*1e3, conductivity[s][:,key], label=val)
        axes[s].legend()
    #plt.show()
    plt.savefig(pdf+"conductivity_damping_dependence_"+label+"_off-diagonal.pdf")
    plt.close()
# }}}

def plot_T_ratio(conductivity, conductivity0, label1, label2, label): # {{{
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    ratio  = [conductivity[0][:,0]]
    ratio0 = [conductivity0[0][:,0]]
    for i in np.arange(1,10):
        ratio.append(conductivity[0][:,i] / conductivity[0][:,1])
        ratio0.append(conductivity0[0][:,i] / conductivity0[0][:,1])
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_title(valleys[0])
    ax.set_ylabel("ratio")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_ylim([0,1.1])
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio[0]*1e3,  ratio[key],  marker=markers[i], s=30, label=label1+" "+val)
        i = i + 1
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio0[0]*1e3, ratio0[key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=label2+" "+val)
        i = i + 1
    data = [1, 2, 0.159]
    padding = [0.01, 0.01, 0.01]
    for i in range(len(data)):
        ax.plot([ratio[0][0], ratio[0][-1]*1e3], [data[i], data[i]])
        ax.text(ratio[0][0], data[i]+padding[i], str(data[i]))
    ax.legend()
#    plt.show()
    plt.savefig(pdf+"conductivity_T_damping_dependence_"+label+"_diagonal_ratio.pdf")
    plt.close()
# }}}
def plot_L3_ratio(conductivity, conductivity0, label1, label2, label): # {{{
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    ratio  = [conductivity[3][:,0]]
    ratio0 = [conductivity0[3][:,0]]
    for i in np.arange(1,10):
        ratio.append(conductivity[3][:,i] / conductivity[3][:,1])
        ratio0.append(conductivity0[3][:,i] / conductivity0[3][:,1])
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_title(valleys[3])
    ax.set_ylabel("ratio")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio[0]*1e3,  ratio[key],  marker=markers[i], s=30, label=label1+" "+val)
        i = i + 1
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio0[0]*1e3, ratio0[key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=label2+" "+val)
        i = i + 1
#    data = [1, 0.02, 0.62]
#    padding = [0.01, 0.01, -0.04]
    data = [1, 0.007, 0.67]
    padding = [0.01, 0.01, -0.04]
    for i in range(len(data)):
        ax.plot([ratio[0][0], ratio[0][-1]*1e3], [data[i], data[i]])
        ax.text(ratio[0][0], data[i]+padding[i], str(data[i]))
    ax.legend()
#    plt.show()
    plt.savefig(pdf+"conductivity_L3_damping_dependence_"+label+"_diagonal_ratio.pdf")
    plt.close()
# }}}
def plot_Lsum_ratio(conductivity, conductivity0, label1, label2, label): # {{{
    j = 5
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    ratio  = [conductivity[j][:,0]]
    ratio0 = [conductivity0[j][:,0]]
    for i in np.arange(1,10):
        ratio.append(conductivity[j][:,i] / conductivity[j][:,1])
        ratio0.append(conductivity0[j][:,i] / conductivity0[j][:,1])
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_title(valleys[j])
    ax.set_ylabel("ratio")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio[0]*1e3,  ratio[key],  marker=markers[i], s=30, label=label1+" "+val)
        i = i + 1
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio0[0]*1e3, ratio0[key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=label2+" "+val)
        i = i + 1
    ax.legend()
#    plt.show()
    plt.savefig(pdf+"conductivity_Lsum_damping_dependence_"+label+"_diagonal_ratio.pdf")
    plt.close()
# }}}
def plot_total_ratio(conductivity, conductivity0, label1, label2, label): # {{{
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    ratio  = [conductivity[0][:,0]]
    ratio0 = [conductivity0[0][:,0]]
    for i in np.arange(1,10):
        ratio.append(conductivity[4][:,i] / conductivity[4][:,1])
        ratio0.append(conductivity0[4][:,i] / conductivity0[4][:,1])
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_title(valleys[4])
    ax.set_ylabel("ratio")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio[0]*1e3,  ratio[key],  marker=markers[i], s=30, label=label1+" "+val)
        i = i + 1
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio0[0]*1e3, ratio0[key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=label2+" "+val)
        i = i + 1
    ax.legend()
#    plt.show()
    plt.savefig(pdf+"conductivity_total_damping_dependence_"+label+"_diagonal_ratio.pdf")
    plt.close()
# }}}

def plot_total_ratio_to_fit_experiment(conductivity, label): # {{{
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    ratio  = [conductivity[0][:,0]]
    amp = conductivity[5][:,1].mean() / conductivity[0][:,1].mean() * (22/56.5)
    for i in np.arange(1,10):
        ratio.append((conductivity[0][:,i] * amp + conductivity[5][:,i]) / (conductivity[0][:,1] * amp + conductivity[5][:,1]))
    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.set_title(valleys[4])
    ax.set_ylabel("ratio")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_ylim([0, 1.1])
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(ratio[0]*1e3,  ratio[key],  marker=markers[i], s=30, label=label+" "+val)
        i = i + 1
    ax.legend()
    print(amp)
    ax.text(ratio[0][0]*1e3, 0.5, "$\gamma_L / \gamma_T = {:.3f}$".format(amp))

    data = 70.5 / 78.5
    padding = - 0.04
    ax.plot([ratio[0][0], ratio[0][-1]*1e3], [data, data])
    ax.text(ratio[0][0], data+padding, "{:.3f}".format(data))

#    plt.show()
    plt.savefig(pdf+"conductivity_total_damping_dependence_"+label+"_amp_diagonal_ratio.pdf")
    plt.close()
# }}}

def plot_diagonal_single(j, conductivity, conductivity0, label1, label2, label): # {{{
    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    maximun = 0e0
    minimun = 1e15
    for s in np.arange(0,5):
        for i in plot_list_conductivity_diagonal.keys():
            maximun = max([maximun, np.abs(conductivity[s][:,i].max()),
                                    np.abs(conductivity[s][:,i].min()),
                                    np.abs(conductivity0[s][:,i].max()),
                                    np.abs(conductivity0[s][:,i].min())])
            minimun = min([minimun, np.abs(conductivity[s][:,i].min()), np.abs(conductivity0[s][:,i].min())])
    window = [minimun*0.9e0, maximun*1.1e0]

    fig, ax = plt.subplots(1,1,figsize=(7,7))

    ax.set_ylabel("conductivity [/Ohm m]")
    ax.set_xlabel("damping constant [meV]")
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_yscale('log')
    ax.set_ylim(window)

    ax.set_title(valleys[j])
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(conductivity[j][:,0]*1e3, conductivity[j][:,key], marker=markers[i], s=30, label=label1+" "+val)
        i = i + 1
    i = 0
    for key, val in plot_list_conductivity_diagonal.items():
        ax.scatter(conductivity0[j][:,0]*1e3, conductivity0[j][:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=label2+" "+val)
        i = i + 1
    ax.legend()
    #plt.show()
    plt.savefig(pdf+"conductivity_"+str(valleys[j])+"_damping_dependence_"+label+"_diagonal.pdf")
    plt.close()
# }}}

def main():
    conductivity  = set_conductivity(12)
    conductivity0 = set_conductivity_constant_gamma(12)
    plot_diagonal(conductivity,  "12bands-V1V2")
    plot_diagonal(conductivity0, "12bands-constant-gamma")
    for i in range(6):
        plot_diagonal_single(i, conductivity, conductivity0, "Born", "const. $\gamma$", "12bands")
    plot_off_diagonal(conductivity,  "12bands-V1V2")
    plot_off_diagonal(conductivity0, "12bands-constant-gamma")
    plot_L3_ratio(conductivity, conductivity0, "Born", "const. $\gamma$", "12bands")
    plot_Lsum_ratio(conductivity, conductivity0, "Born", "const. $\gamma$", "12bands")
    plot_T_ratio(conductivity, conductivity0, "Born", "const. $\gamma$", "12bands")
    plot_total_ratio(conductivity, conductivity0, "Born", "const. $\gamma$", "12bands")

    conductivity_L4bands  = set_conductivity(4)
    conductivity0_L4bands = set_conductivity_constant_gamma(4)
    for i in range(6):
        plot_diagonal_single(i, conductivity_L4bands, conductivity0_L4bands, "Born", "const. $\gamma$", "4bands")
    plot_diagonal(conductivity_L4bands,  "4bands-V1V2")
    plot_diagonal(conductivity0_L4bands, "4bands-constant-gamma")
    plot_L3_ratio(conductivity_L4bands, conductivity0_L4bands, "Born", "const. $\gamma$", "4bands")
    plot_Lsum_ratio(conductivity_L4bands, conductivity0_L4bands, "Born", "const. $\gamma$", "4bands")
    plot_T_ratio(conductivity_L4bands, conductivity0_L4bands, "Born", "const. $\gamma$", "4bands")
    plot_total_ratio(conductivity_L4bands, conductivity0_L4bands, "Born", "const. $\gamma$", "4bands")

    plot_L3_ratio(conductivity, conductivity_L4bands, "12bands", "4bands", "12bands_vs_4bands")

    plot_total_ratio_to_fit_experiment(conductivity0, "12bands_constant-gamma")
    plot_total_ratio_to_fit_experiment(conductivity0_L4bands, "4bands_constant-gamma")

if __name__ == "__main__":
    main()
