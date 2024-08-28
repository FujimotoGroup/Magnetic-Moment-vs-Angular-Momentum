import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import pandas as pd
import configparser as cnf
from matplotlib.ticker import ScalarFormatter
from scipy import interpolate
from decimal import Decimal, Context

mpl.rcParams["font.family"] = 'Times New Roman'
mpl.rcParams['mathtext.fontset'] = 'cm'
#mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['lines.linewidth'] = 2.0
mpl.rcParams["font.size"] = 16
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc('text.latex', preamble=r'\usepackage{bm}')

figure_size_unit = 8
bbox = {
    "facecolor" : "white",
    "edgecolor" : "black",
    "boxstyle" : "round, pad=0.5",
    "linewidth" : 2
}

home = "../"
data0 = "../dat/"
png = "../fig/png/"
svg = "../fig/svg/"
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

markers = ["o", "v", "D", ",", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

epsilon = 5e-4
label = "{:.6f}".format(epsilon)
str_damping = "Gamma = "+"{:.2f}".format(epsilon*1e3)+" [meV]"

labels = {"T,4":0, "L1,0":1, "L1,2":2, "L2,0":3, "L2,2":4, "L3,0":5, "L3,2":6}
valleys = {"T":0, "L1":1, "L2":2, "L3":3}
bands  = 12
bandsT = 12
bandsL = 4
lower_band_L = 0
upper_band_L = 2
#bandsL = 12
#lower_band_L = 4
#upper_band_L = 6

axises = {"x":0, "y":1, "z":2}
k_cutoff = 0.27722e0

cutoffs = ["0.08", "0.10"]
cutoff = cutoffs[1]
params = ["cutoff"+cutoff+"eV" for cutoff in cutoffs]
param = "cutoff"+cutoff
n = 16000.0
e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
x = e_range / n

def decimal_normalize(f, r='0.001'):
    """数値fの誤差をrの位で四捨五入する"""
    f = Decimal(str(f)).quantize(Decimal(r))
    """数値fの小数点以下を正規化し、文字列で返す"""
    def _remove_exponent(d):
        return d.quantize(Decimal(1)) if d == d.to_integral() else d.normalize()
    a = Decimal.normalize(Decimal(str(f)))
    b = _remove_exponent(a)
    return str(b)

xticks = np.linspace(-float(cutoff), float(cutoff), 4+1)
xticklabels = ["$"+decimal_normalize(tick)+"$" for tick in xticks]

def plot_dispersion(point, full, effective, output): # dispersion {{{
    fig, axes = plt.subplots(1,3,figsize=(15,7))
    fig.tight_layout()
    axes = axes.flatten()
    axes[0].set_ylabel("Energy")
    axes[1].set_title(point+" point")
    axes[1].tick_params(labelleft=False)
    axes[2].tick_params(labelleft=False)
    for axis, val in axises.items():
        axes[val].set_xlabel("$k_"+axis+"$")
        axes[val].set_ylim(-0.2,0.2)
#        axes[val].axhspan(-float(cutoff), float(cutoff), color="gray", alpha=0.3)
        axes[val].axvspan(-k_cutoff, k_cutoff, color="tab:cyan", alpha=0.2)

    for axis, val in full.items():
        for i in np.arange(1, val.shape[1]):
            (p1,) = axes[axises[axis]].plot(val[:,0], val[:,i], color=colors[0])
    for axis, val in full.items():
        for i in np.arange(1, val.shape[1]):
            p2 = axes[axises[axis]].scatter(val[:,0], val[:,i], color=colors[7])

    for ax in axes:
        ax.legend([p1, p2], ["Liu-Allen", "effective"])
    #plt.show()

    plt.subplots_adjust(wspace=0.05, hspace=0.35)

    plt.savefig(png+output+".png", bbox_inches='tight', dpi=300)
    plt.savefig(pdf+output+".pdf")
    plt.close()
# }}}

def plot_dos(e_max, e_min, dos, output): # dos {{{
    e_max = sorted(set(e_max))
    e_min = sorted(set(e_min))
    tmp = []
    for d in dos:
        ip1d = interpolate.interp1d(d[:,0], d[:,1])
        tmp.append(ip1d(x))
    dos_sum = np.sum(tmp, axis=0)

    fig, ax = plt.subplots(1,1,figsize=(7,7))
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.set_xlabel("$\mu~\mathrm{[eV]}$")
    ax.set_ylabel("$\mathrm{DOS~[/eV m^3]}$")

    ax.plot(x, dos_sum, label="total", c="red")
    for d in tmp:
        ax.plot(x, d, c="black", lw=1)
    for key, val in valleys.items():
        ax.scatter(dos[val][:,0], dos[val][:,1], s=4, label=key)
    ax.legend()
    #plt.show()

    plt.savefig(png+output+".png", bbox_inches = 'tight', dpi=300)
    plt.savefig(pdf+output+".pdf")
    plt.close()
# }}}

def plot_conductivity_valleys(titles, conductivity, conductivity_valley, output): # conductivity {{{
    conductivity_i = []
    conductivity_i_sum = []

    fig, axes = plt.subplots(1,3,figsize=(15,7))
    axes = axes.flatten()

    plot_lists = [1, 5, 9]
    for s in np.arange(0,3):
        conductivity_tmp = []
        axes[s].set_title("$\sigma_{"+titles[plot_lists[s]]+"}$")
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

    axes[0].set_ylabel("$\mathrm{conductivity~[/\Omega\ m]}$")
    for ax in axes:
        ax.set_xlabel("$\mu~\mathrm{[eV]}$")
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

    axes[-1].legend()
    #plt.show()

    plt.savefig(png+output+".png", bbox_inches = 'tight', dpi=300)
    plt.savefig(pdf+output+".pdf")
    #plt.show()
    plt.close()
# }}}

def plot_conductivity(conductivity, conductivity_valley, output): # individual + total electric conductivity {{{
    n = 100e0
    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
    x = e_range / n

    label_cutoff = 'cutoff'+cutoff+'eV'
    label_eps = "{:.6f}".format(epsilon)
    str_damping = "$\gamma = "+decimal_normalize(epsilon*1e3)+"~\mathrm{meV}$"

    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    plot_list_conductivity_offdiagonal = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}

    sigma_T = conductivity_valley[0]

    e = conductivity_valley[1][:,0]
    conductivity = conductivity_valley[1] + conductivity_valley[2] + conductivity_valley[3]
    conductivity[:,0] = e
    sigma_L = conductivity

    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(x) for i in np.arange(1,10)])
    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(x) for i in np.arange(1,10)])

    conductivity = np.array([is_T, is_L])
    conductivity_total = np.sum(conductivity, axis=0)
    is_T = np.concatenate([np.array([x]), is_T])
    is_L = np.concatenate([np.array([x]), is_L])
    conductivity_total = np.concatenate([np.array([x]), conductivity_total])
    conductivity = np.array([is_T, is_L, conductivity_total])

    window_diagonal = [0, 0]
    window_offdiagonal = [0, 0]
    for key in plot_list_conductivity_diagonal.keys():
        window_diagonal[0] = min(window_diagonal[0], is_T[key].min(), is_L[key].min())
        window_diagonal[1] = max(window_diagonal[1], is_T[key].max(), is_L[key].max())

    for key in plot_list_conductivity_offdiagonal.keys():
        window_offdiagonal[0] = min(window_offdiagonal[0], is_T[key].min(), is_L[key].min())
        window_offdiagonal[1] = max(window_offdiagonal[1], is_T[key].max(), is_L[key].max())

    window_diagonal = np.array([-window_diagonal[1]*0.02, window_diagonal[1]*1.1])
    window_offdiagonal = np.array(window_offdiagonal)*10

    col_num = 2
    row_num = 3
    fig, axes = plt.subplots(col_num,row_num,figsize=(figure_size_unit*row_num,figure_size_unit*col_num))
    axes = axes.T

    for ax_col in axes:
        for ax in ax_col:
            ax.set_xlim([-float(cutoff), float(cutoff)])
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.grid()
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
    for ax in axes.T[1]:
        ax.set_xlabel("$\mu~[\mathrm{eV}]$")

    axes[0][0].set_ylabel("$\mathrm{longitudinal\ conductivity~[/\Omega \, m]}$")
    axes[0][1].set_ylabel("$\mathrm{transverse\ conductivity  ~[/\Omega \, m]}$")

    labels = ["T", "L", "Total"]
    for n in np.arange(0,3):
        axes[n][0].text(-float(cutoff)*0.9, window_diagonal[1]*0.93, labels[n], bbox=bbox, fontsize="large")
        axes[n][1].text(-float(cutoff)*0.9, window_offdiagonal[1]*0.89, labels[n], bbox=bbox, fontsize="large")
        i = 1
        for key, val in plot_list_conductivity_diagonal.items():
            axes[n][0].set_ylim(window_diagonal)
            axes[n][0].scatter(x, conductivity[n][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
            i = i + 1

        for key, val in plot_list_conductivity_offdiagonal.items():
            axes[n][1].set_ylim(window_offdiagonal)
            axes[n][1].scatter(x, conductivity[n][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
            i = i + 1

    axes[0][0].legend()
    axes[0][1].legend()

#    plt.show()
    fig.tight_layout()
    plt.savefig(pdf+output+".pdf", bbox_inches='tight')
    plt.close()
# }}}

def plot_spin_conductivity_valleys(sign, titles, conductivity, conductivity_valley, output): # spin conductivity {{{
    conductivity_i = []
    conductivity_i_sum = []

    fig, axes = plt.subplots(1,3,figsize=(15,7))
    axes = axes.flatten()

    plot_lists = [6, 16, 20]
    for s in np.arange(0,3):
        conductivity_tmp = []
        axes[s].set_title(titles[plot_lists[s]])
        force = titles[plot_lists[s]][1]
        flow  = titles[plot_lists[s]][2]
        spin  = titles[plot_lists[s]][3]
        axes[s].set_title("$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
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

    window = [-maximun*1.1e0, maximun*1.1e0]

    axes[0].set_ylabel("$\mathrm{spin\ conductivity~[/\Omega\ m]}$")
    for ax in axes:
        ax.set_xlabel("$\mu~\mathrm{[eV]}$")
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.set_xlim(-float(cutoff), float(cutoff))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_ylim(window)

    for s in np.arange(0,3):
#        axes[s].plot(x, conductivity_i_sum[s], label="total", c="red")
#        for d in conductivity_i[s]:
#            axes[s].plot(x, d, c="black", lw=1)
        for key, val in labels.items():
            axes[s].scatter(conductivity[val][:,0], conductivity[val][:,plot_lists[s]], s=4, label=key)
    axes[-1].legend()
    #plt.show()

    plt.savefig(png+output+".png", bbox_inches='tight', dpi=300)
    plt.savefig(pdf+output+".pdf")
    #plt.show()
    plt.close()
# }}}

def plot_spin_conductivity(sign, conductivity1_valley, conductivity2_valley, output): # individual spin conductivity {{{
    n = 100e0
    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
    x = e_range / n

    label_cutoff = 'cutoff'+cutoff+'eV'
    label_eps = "{:.6f}".format(epsilon)
    str_damping = "$\gamma = "+decimal_normalize(epsilon*1e3)+"~\mathrm{meV}$"

    plot_list1  = { 6:"xyz", 16:"yzx", 20:"zxy"}
    plot_list2  = {12:"yxz",  8:"xzy", 22:"zyx"}
    plot_list3  = { 1:"xxx", 13:"yyx",  5:"xyy", 11:"yxy"}
    plot_list4a = { 2:"xxy",  3:"xxz",  4:"xyx",  7:"xzx",  9:"xzz", 10:"yxx", 14:"yyy", 15:"yyz", 17:"yzy"}
    plot_list4b = {18:"yzz", 19:"zxx", 21:"zxz", 23:"zyy", 24:"zyz", 25:"zzx", 26:"zzy", 27:"zzz"}
    plot_lists = [ plot_list1, plot_list2, plot_list3, plot_list4a, plot_list4b ]

    sigma_T1 = conductivity1_valley[0]
    sigma_T2 = conductivity2_valley[0]

    e = conductivity1_valley[1][:,0]
    conductivity = conductivity1_valley[1] + conductivity1_valley[2] + conductivity1_valley[3]
    conductivity[:,0] = e
    sigma_L1 = conductivity
    e = conductivity2_valley[1][:,0]
    conductivity = conductivity2_valley[1] + conductivity2_valley[2] + conductivity2_valley[3]
    conductivity[:,0] = e
    sigma_L2 = conductivity

    is_T1 = np.array([interpolate.interp1d(sigma_T1[:,0], sigma_T1[:,i])(x) for i in np.arange(1,28)])
    is_T2 = np.array([interpolate.interp1d(sigma_T2[:,0], sigma_T2[:,i])(x) for i in np.arange(1,28)])
    is_L1 = np.array([interpolate.interp1d(sigma_L1[:,0], sigma_L1[:,i])(x) for i in np.arange(1,28)])
    is_L2 = np.array([interpolate.interp1d(sigma_L2[:,0], sigma_L2[:,i])(x) for i in np.arange(1,28)])

    conductivity_T = np.sum(np.array([is_T1, is_T2]), axis=0)
    conductivity_L = np.sum(np.array([is_L1, is_L2]), axis=0)
    conductivity_total = np.sum(np.array([conductivity_T, conductivity_L]), axis=0)
    is_T = np.concatenate([np.array([x]), conductivity_T])
    is_L = np.concatenate([np.array([x]), conductivity_L])
    conductivity_total = np.concatenate([np.array([x]), conductivity_total])
    conductivity = np.array([is_T, is_L, conductivity_total])

    window = []
    for j, plot_list in enumerate(plot_lists):
        tmp = [0, 0]
        for n in np.arange(0,3):
            for key in plot_list.keys():
                tmp[0] = min(tmp[0], conductivity[n][key].min())
                tmp[1] = max(tmp[1], conductivity[n][key].max())
        tmp = np.array(tmp) * 1.1e0
        window.append(tmp)

    window[3] = window[3]*10e0
    window[4] = window[4]*10e0

    col_num = 3
    row_num = 5
    fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
    axes = axes.T

    for ax_list in axes:
        for ax in ax_list:
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.grid()
            ax.set_xlim(-float(cutoff), float(cutoff))
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
    for ax in axes.T[4]:
        ax.set_xlabel("$\mu~[\mathrm{eV}]$")
    for ax in axes[0]:
        ax.set_ylabel("$\mathrm{spin\ conductivity~[/\Omega \, m]}$")

    labels = ["T", "L", "Total"]
    for n in np.arange(0,3):
        k = 0
        for j, plot_list in enumerate(plot_lists):
            i = 0
            axes[n][k].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][k].transAxes)
            for key, val in plot_list.items():
                axes[n][j].set_ylim(window[j])
                force = val[0]
                flow  = val[1]
                spin  = val[2]
                axes[n][j].scatter(x, conductivity[n][key], marker=markers[i], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
                i += 1
            axes[0][j].legend()
            k += 1

#    plt.show()
    fig.tight_layout()
    fig.savefig(pdf+output+".pdf", bbox_inches='tight')
    plt.close()
# }}}

def plot_total_spin_conductivity(sign, titles, conductivity1_valley, conductivity2_valley, output): # {{{
    conductivity1_i = []
    conductivity1_i_sum = []

    plot_lists = [6, 16, 20]

    for s in np.arange(0,3):
        conductivity1_tmp = []
        for d in conductivity1_valley:
            ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
            conductivity1_tmp.append(ip1d(x))
        conductivity1_tmp = np.array(conductivity1_tmp)
        conductivity1_i.append(conductivity1_tmp)
        conductivity1_tmp = np.sum(conductivity1_tmp, axis=0)
        conductivity1_i_sum.append(conductivity1_tmp)

    conductivity2_i = []
    conductivity2_i_sum = []

    for s in np.arange(0,3):
        conductivity2_tmp = []
        for d in conductivity2_valley:
            ip1d = interpolate.interp1d(d[:,0], d[:,plot_lists[s]])
            conductivity2_tmp.append(ip1d(x))
        conductivity2_tmp = np.array(conductivity2_tmp)
        conductivity2_i.append(conductivity2_tmp)
        conductivity2_tmp = np.sum(conductivity2_tmp, axis=0)
        conductivity2_i_sum.append(conductivity2_tmp)

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

    axes[0].set_ylabel("$\mathrm{spin\ conductivity~[/\Omega\ m]}$")
    for ax in axes:
        ax.set_xlabel("$\mu~\mathrm{[eV]}$")
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.set_xlim(-float(cutoff), float(cutoff))
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_ylim(window)

    for s in np.arange(0,3):
        force = titles[plot_lists[s]][1]
        flow  = titles[plot_lists[s]][2]
        spin  = titles[plot_lists[s]][3]
        axes[s].set_title("$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
        axes[s].plot(x, conductivity1_i_sum[s], label="SHC1",  c=colors[4])
        axes[s].plot(x, conductivity2_i_sum[s], label="SHC2",  c=colors[6])
        axes[s].plot(x, conductivity3_i_sum[s], label="total", c="red")
    axes[-1].legend()

    plt.savefig(png+output+".png", bbox_inches='tight', dpi=300)
    plt.savefig(pdf+output+".pdf")
    plt.close()

# }}}

def main():
    parameter = "T"+str(bandsT)+"bands_L"+str(bandsL)+"bands_gamma"+label+"_"+param+"eV"

# dispersion T point {{{
    full = {}
    effective = {}
    for axis, val in axises.items():
        read_full = data0+'EigenValue_T_'+axis+'_full.dat'
        full[axis] = pd.read_csv(read_full,header=None).values
        read_effective = data0+'EigenValue_T_'+axis+'_'+str(bandsT)+'bands.dat'
        effective[axis] = pd.read_csv(read_effective,header=None).values
    plot_dispersion("T", full, effective, "dispersion_T-"+str(bandsT)+"bands")
# }}}

# dispersion L points {{{
    for valley in np.arange(1,4):
        full = {}
        effective = {}
        for axis, val in axises.items():
            read_full = data0+'EigenValue_L'+str(valley)+'_'+axis+'_full.dat'
            full[axis] = pd.read_csv(read_full,header=None).values
            read_effective = data0+'EigenValue_L'+str(valley)+'_'+axis+'_'+str(bandsL)+'bands.dat'
            effective[axis] = pd.read_csv(read_effective,header=None).values
        plot_dispersion("L"+str(valley), full, effective, "dispersion_L"+str(valley)+"-"+str(bandsL)+"bands")
# }}}

# dos {{{
    e_max = []
    e_min = []
    dos = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    readfile = data+'dos_eps'+label+'.csv'
    d = pd.read_csv(readfile,header=None).values
    dos.append(d)
    e_max.append(d[:,0].max())
    e_min.append(d[:,0].min())

    for valley in np.arange(1,4):
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        readfile = data+'dos_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=None).values
        e = d4[:,0]
        d4 = - d4
        d4[:,0] = e
        e_max.append(d4[:,0].max())
        e_min.append(d4[:,0].min())

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        readfile = data+'dos_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=None).values
        e_max.append(d6[:,0].max())
        e_min.append(d6[:,0].min())

        dos.append(np.append(d4, d6, axis=0))

    plot_dos(e_max,e_min,dos,"dos_"+parameter)
# }}}

# conductivity {{{
    conductivity = []
    conductivity_valley = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity.append(d)
    conductivity_valley.append(d)

    for valley in np.arange(1,4):
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        e = d4[:,0]
        d4 = - d4
        d4[:,0] = e
        conductivity.append(d4)

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values
        conductivity.append(d6)

        conductivity_valley.append(np.append(d4, d6, axis=0))

    plot_conductivity_valleys(titles,conductivity,conductivity_valley,"conductivity"+parameter)
    plot_conductivity(conductivity,conductivity_valley,"conductivity_tensor"+parameter)
# }}}

# spin magnetic conductivity1 {{{
    conductivity1 = []
    conductivity1_valley = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
    print(readfile)
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity1.append(d)
    conductivity1_valley.append(d)

    for valley in np.arange(1,4):
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
        print(readfile)
        d4 = pd.read_csv(readfile,header=0).values
        e = d4[:,0]
        d4 = - d4
        d4[:,0] = e
        conductivity1.append(d4)

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
        print(readfile)
        d6 = pd.read_csv(readfile,header=0).values
        conductivity1.append(d6)

        conductivity1_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("m", titles,conductivity1,conductivity1_valley,"spin_magnetic_conductivity1_"+parameter)
# }}}

# spin magnetic conductivity2 {{{
    conductivity2 = []
    conductivity2_valley = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
    d = pd.read_csv(readfile,header=0).values
    conductivity2.append(d)
    conductivity2_valley.append(d)

    for valley in np.arange(1,4):
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        e = d4[:,0]
        d4 = - d4
        d4[:,0] = e
        conductivity2.append(d4)

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values

        val = np.copy(d4[-1, :])
        val[0] = 0e0
        d6 = d6 + val

        conductivity2.append(d6)
        conductivity2_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("m", titles,conductivity2,conductivity2_valley,"spin_magnetic_conductivity2_"+parameter)
# }}}

# spin magnetic conductivity 1 + 2 {{{
    plot_spin_conductivity("m", conductivity1_valley, conductivity2_valley, "spin_magnetic_conductivity_tensor_"+parameter)
    plot_total_spin_conductivity("m", titles, conductivity1_valley, conductivity2_valley, "total_spin_magnetic_conductivity_"+parameter)
# }}}

# spin angular conductivity1 {{{
    conductivity1 = []
    conductivity1_valley = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity1.append(d)
    conductivity1_valley.append(d)

    for valley in np.arange(1,4):
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        e = d4[:,0]
        d4 = - d4
        d4[:,0] = e
        conductivity1.append(d4)

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values
        conductivity1.append(d6)

        conductivity1_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("a", titles,conductivity1,conductivity1_valley,"spin_angular_conductivity1_"+parameter)
# }}}

# spin angular conductivity2 {{{
    conductivity2 = []
    conductivity2_valley = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    d = pd.read_csv(readfile,header=0).values
    conductivity2.append(d)
    conductivity2_valley.append(d)

    for valley in np.arange(1,4):
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        e = d4[:,0]
        d4 = - d4
        d4[:,0] = e
        conductivity2.append(d4)

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values

        val = np.copy(d4[-1, :])
        val[0] = 0e0
        d6 = d6 + val

        conductivity2.append(d6)
        conductivity2_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("a", titles,conductivity2,conductivity2_valley,"spin_angular_conductivity2_"+parameter)
# }}}

# spin angular conductivity 1 + 2 {{{
    plot_spin_conductivity("a", conductivity1_valley, conductivity2_valley, "spin_angular_conductivity_tensor_"+parameter)
    plot_total_spin_conductivity("a", titles, conductivity1_valley, conductivity2_valley, "total_spin_angular_conductivity_"+parameter)
# }}}

if __name__ == "__main__":
    main()
