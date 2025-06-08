import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches
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
mpl.rcParams["font.size"] = 28
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
plt.rc('text.latex', preamble=r'\usepackage{bm}')

figure_size_unit = 5
bbox = {
    "facecolor" : "white",
    "edgecolor" : "black",
    "boxstyle" : "round, pad=0.5",
    "linewidth" : 1
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

epsilon_T = 0.7e-4
epsilon_L = 1.6e-4

cutoffs = ["0.08", "0.10"]
cutoff = cutoffs[0]
params = ["cutoff"+cutoff+"eV" for cutoff in cutoffs]
param = "cutoff"+cutoff
n = 16000.0
e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
x = e_range / n

key_L   = [2 , 3, 1]

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
xticklabels = ["$"+decimal_normalize(tick*1e3)+"$" for tick in xticks]

def interpolate_conductivity(conductivity1_valley, conductivity2_valley): # {{{
    sigma_T1 = conductivity1_valley[0]
    sigma_T2 = conductivity2_valley[0]

    n = 100e0
    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
    x = e_range / n

    e = conductivity1_valley[1][:,0]
    conductivity = conductivity1_valley[1] + conductivity1_valley[2] + conductivity1_valley[3]
    conductivity[:,0] = e
    sigma_L1 = conductivity
    e = conductivity2_valley[1][:,0]
    conductivity = conductivity2_valley[1] + conductivity2_valley[2] + conductivity2_valley[3]
    conductivity[:,0] = e
    sigma_L2 = conductivity

    n = sigma_T1.shape[1]
    is_T1 = np.array([interpolate.interp1d(sigma_T1[:,0], sigma_T1[:,i])(x) for i in np.arange(1,n)])
    is_T2 = np.array([interpolate.interp1d(sigma_T2[:,0], sigma_T2[:,i])(x) for i in np.arange(1,n)])
    is_L1 = np.array([interpolate.interp1d(sigma_L1[:,0], sigma_L1[:,i])(x) for i in np.arange(1,n)])
    is_L2 = np.array([interpolate.interp1d(sigma_L2[:,0], sigma_L2[:,i])(x) for i in np.arange(1,n)])

    valley_resolved = []
    for i in np.arange(4):
        sigma = conductivity1_valley[i] + conductivity2_valley[i]
        if i == 0:
            sigma[:,0] = conductivity1_valley[0][:,0]
        else:
            sigma[:,0] = conductivity1_valley[1][:,0]
        tmp = np.concatenate([np.array([x]), np.array([interpolate.interp1d(sigma[:,0], sigma[:,i])(x) for i in np.arange(1,n)])])
        valley_resolved.append(tmp)
    tmp = np.concatenate([np.array([x]), np.sum(np.array([is_L1, is_L2]), axis=0)])
    valley_resolved.append(tmp)
    np.array(valley_resolved)

    conductivity_T = np.sum(np.array([is_T1, is_T2]), axis=0)
    conductivity_L = np.sum(np.array([is_L1, is_L2]), axis=0)
    conductivity_total = np.sum(np.array([conductivity_T, conductivity_L]), axis=0)
    conductivity_total = np.concatenate([np.array([x]), conductivity_total])
    conductivity = np.array([conductivity_T, conductivity_L])

    return valley_resolved, conductivity_total
# }}}

def plot_dispersion(point, full, effective, output): # dispersion {{{
    fig, axes = plt.subplots(1,3,figsize=(figure_size_unit*3,figure_size_unit))
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

    fig, ax = plt.subplots(1,1,figsize=(figure_size_unit, figure_size_unit))
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

    fig, axes = plt.subplots(1,3,figsize=(figure_size_unit*3, figure_size_unit))
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

def plot_conductivity_valley_resolved(titles, conductivity, conductivity_valley, output): # conductivity {{{
    figure_size_unit = 6

    n = 100e0
    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
    x = e_range / n

    n = conductivity_valley[0].shape[1]
    valley_resolved = []
    s = np.zeros((n-1, 99))
    for i in np.arange(4):
        sigma = conductivity_valley[i]
        if i == 0:
            sigma[:,0] = conductivity_valley[0][:,0]
        else:
            sigma[:,0] = conductivity_valley[1][:,0]
        s_tmp = np.array([interpolate.interp1d(sigma[:,0], sigma[:,i])(x) for i in np.arange(1,n)])
        tmp = np.concatenate([np.array([x]), s_tmp])
        valley_resolved.append(tmp)
        if i > 0:
            s += s_tmp
    tmp = np.concatenate([np.array([x]), s])
    valley_resolved.append(tmp)
    conductivity = np.array(valley_resolved)
    window = [[-1e5, 9e6], [-2e6,2e6]]

    plot_list_dia = {1:"xx", 5:"yy", 9:"zz"}
    plot_list_off = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}
    plot_list = [plot_list_dia, plot_list_off]

    col_num = 5
    row_num = 2
    fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
    axes = axes.T

    for ax_list in axes:
        for ax in ax_list:
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.grid()
            ax.set_axisbelow(True)
            ax.set_xlim(-float(cutoff), float(cutoff))
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.tick_params(axis='x', labelbottom=False)  # ラベルだけ消す
    for ax in axes.T[-1]:
        ax.set_xlabel("$\mu~[\mathrm{meV}]$")
        ax.tick_params(axis='x', labelbottom=True)

    axes[0][0].set_ylabel(r"$\sigma_{ii}~\mathrm{[/\Omega \, m]}$")
    axes[0][1].set_ylabel(r"$\sigma_{ij}~\mathrm{[/\Omega \, m]}$")

    loc = [["lower right", "upper right"], ["upper right", "center right"]]
    for n in np.arange(0,2):
        for ax in axes.T[n]:
            ax.set_ylim(window[n])

        for j in np.arange(5):
            i = 0
            sigma = conductivity[j]
            for key, val in plot_list[n].items():
                force = val[0]
                flow  = val[1]
                axes.T[n][j].scatter(sigma[0], sigma[key], marker=markers[i], s=20, edgecolors=colors[i+1], facecolor='None', label="$\sigma_{"+flow+force+"}$")
                i += 1
#            axes.T[n][j].legend(labelspacing=0.2, handletextpad=0.2, loc=loc[n][j], markerscale=2)

    axes.T[0][0].legend(labelspacing=0.2, loc=loc[0][1], markerscale=2)
    axes.T[1][0].legend(labelspacing=0.2, loc=loc[1][0], bbox_to_anchor = (1, 0.95), markerscale=2)

    ls = ["$h$", "$e_1$", "$e_2$", "$e_3$", "$e_1+e_2+e_3$"]
    for n in np.arange(5):
        axes.T[0][n].text(0.07, 0.8, ls[n], bbox=bbox, transform=axes.T[0][n].transAxes, alpha=0.8)
        axes.T[1][n].text(0.07, 0.8, ls[n], bbox=bbox, transform=axes.T[1][n].transAxes, alpha=0.8)

    fig.tight_layout()
    plt.savefig(png+output+"_valley.png", bbox_inches='tight', dpi=300)
    plt.savefig(pdf+output+"_valley.pdf")
    plt.close()
# }}}
#
#def plot_conductivity(conductivity, conductivity_valley, output): # individual + total electric conductivity {{{
#    n = 100e0
#    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
#    x = e_range / n
#
#    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
#    plot_list_conductivity_offdiagonal = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}
#
#    conductivity_i = []
#    for d in conductivity_valley:
#        tmp = np.array([interpolate.interp1d(d[:,0], d[:,i])(x) for i in np.arange(1,10)])
#        tmp = np.concatenate([np.array([x]), tmp])
#        conductivity_i.append(tmp)
#    conductivity_i = np.array(conductivity_i)
#
#    sigma_T = conductivity_valley[0]
#
#    e = conductivity_valley[1][:,0]
#    conductivity = conductivity_valley[1] + conductivity_valley[2] + conductivity_valley[3]
#    conductivity[:,0] = e
#    sigma_L = conductivity
#
#    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(x) for i in np.arange(1,10)])
#    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(x) for i in np.arange(1,10)])
#
#    conductivity = np.array([is_T, is_L])
#    conductivity_total = np.sum(conductivity, axis=0)
#    is_T = np.concatenate([np.array([x]), is_T])
#    is_L = np.concatenate([np.array([x]), is_L])
#    conductivity_total = np.concatenate([np.array([x]), conductivity_total])
#    conductivity = np.array([is_T, is_L, conductivity_total])
#
#    print("mu =", conductivity[2][0][49])
#    print("- - - - - - - - - - - - - - -")
#    print("T xx =", conductivity[0][1][49])
#    print("T yy =", conductivity[0][5][49])
#    print("T zz =", conductivity[0][9][49])
#    print("- - - - - - - - - - - - - - -")
#    print("L xx =", conductivity[1][1][49])
#    print("L yy =", conductivity[1][5][49])
#    print("L zz =", conductivity[1][9][49])
#    print("- - - - - - - - - - - - - - -")
#    print("sum xx =", conductivity[2][1][49])
#    print("sum yy =", conductivity[2][5][49])
#    print("sum zz =", conductivity[2][9][49])
#
#    window_diagonal = [0, 0]
#    max_offdiagonal = 0
#    for key in plot_list_conductivity_diagonal.keys():
#        window_diagonal[0] = min(window_diagonal[0], is_T[key].min(), is_L[key].min())
#        window_diagonal[1] = max(window_diagonal[1], is_T[key].max(), is_L[key].max())
#
#    for key in plot_list_conductivity_offdiagonal.keys():
#        max_offdiagonal = max(max_offdiagonal, np.abs(is_T[key].min()), np.abs(is_L[key].min()), np.abs(is_T[key].max()), np.abs(is_L[key].max()))
#
#    window_diagonal = np.array([-window_diagonal[1]*0.02, window_diagonal[1]*1.1])
#    window_offdiagonal = np.array([-max_offdiagonal, max_offdiagonal])*5
#
#    window_L_offdiagonal = [0, 0]
#    for valley in np.arange(1,4):
#        i = valley
#        for key in plot_list_conductivity_offdiagonal.keys():
#            window_L_offdiagonal[0] = min(window_L_offdiagonal[0], conductivity_i[i][key].min())
#            window_L_offdiagonal[1] = max(window_L_offdiagonal[1], conductivity_i[i][key].max())
#
#    window_L_offdiagonal = np.array(window_L_offdiagonal) * 1.1
#
#    col_num = 4
#    row_num = 3
#    fig, axes = plt.subplots(col_num,row_num,figsize=(figure_size_unit*row_num,figure_size_unit*col_num))
#    axes = axes.T
#
#    for ax_col in axes:
#        for ax in ax_col:
#            ax.set_xlim([-float(cutoff), float(cutoff)])
#            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
#            ax.grid()
#            ax.set_xticks(xticks)
#            ax.set_xticklabels(xticklabels)
#    for ax in axes.T[-1]:
#        ax.set_xlabel("$\mu~[\mathrm{eV}]$")
#
#    axes[0][0].set_ylabel("$\mathrm{longitudinal\ conductivity~[/\Omega \, m]}$")
#    axes[0][1].set_ylabel("$\mathrm{transverse\ conductivity  ~[/\Omega \, m]}$")
#    axes[0][2].set_ylabel("$\mathrm{longitudinal\ conductivity~[/\Omega \, m]}$")
#    axes[0][3].set_ylabel("$\mathrm{transverse\ conductivity  ~[/\Omega \, m]}$")
#
#    labels = ["$h$", "$e_1 + e_2 + e_3$", "$h + e_1 + e_2 + e_3$"]
#    for n in np.arange(0,3):
#        axes[n][0].text(0.07, 0.87, labels[n], bbox=bbox, transform=axes[n][0].transAxes)
#        axes[n][1].text(0.07, 0.87, labels[n], bbox=bbox, transform=axes[n][1].transAxes)
#        i = 1
#        for key, val in plot_list_conductivity_diagonal.items():
#            axes[n][0].set_ylim(window_diagonal)
#            axes[n][0].scatter(x, conductivity[n][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
#            i = i + 1
#
#        for key, val in plot_list_conductivity_offdiagonal.items():
#            axes[n][1].set_ylim(window_offdiagonal)
#            axes[n][1].scatter(x, conductivity[n][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
#            i = i + 1
#
#    labels = ["$e_1$", "$e_2$", "$e_3$"]
#    for valley in np.arange(1,4):
#        n = valley - 1
#        axes[n][2].text(0.07, 0.87, labels[n], bbox=bbox, transform=axes[n][2].transAxes)
#        axes[n][3].text(0.07, 0.87, labels[n], bbox=bbox, transform=axes[n][3].transAxes)
#
#        i = 1
#        axes[n][2].set_ylim(window_diagonal)
#        for key, val in plot_list_conductivity_diagonal.items():
#            axes[n][2].scatter(conductivity_i[valley][0], conductivity_i[valley][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
#            i = i + 1
#
#        axes[n][3].set_ylim(window_L_offdiagonal)
#        for key, val in plot_list_conductivity_offdiagonal.items():
#            axes[n][3].scatter(conductivity_i[valley][0], conductivity_i[valley][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
#            i = i + 1
#
##    ax_pos = axes[2][0].get_position()
##    fig.text(ax_pos.x1, ax_pos.y1, str_damping, ha='right', va='bottom')
##    fig.text(1, 1, str_damping, ha='right', va='bottom',transform=axes[2][0].transAxes)
#
#
#    for i in np.arange(col_num):
#        axes[0][i].legend(handletextpad=0.3, markerscale=2)
#
##    plt.show()
#    fig.align_ylabels()
#    fig.tight_layout()
#    plt.savefig(pdf+output+".pdf", bbox_inches='tight')
#    plt.close()
## }}}
#
def plot_conductivity(conductivity, conductivity_valley, output): # individual + total electric conductivity {{{
    n = 100e0
    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
    x = e_range / n

    plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
    plot_list_conductivity_offdiagonal = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}

    conductivity_i = []
    for d in conductivity_valley:
        tmp = np.array([interpolate.interp1d(d[:,0], d[:,i])(x) for i in np.arange(1,10)])
        tmp = np.concatenate([np.array([x]), tmp])
        conductivity_i.append(tmp)
    conductivity_i = np.array(conductivity_i)

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

    print("mu =", conductivity[2][0][49])
    print("- - - - - - - - - - - - - - -")
    print("T xx =", conductivity[0][1][49])
    print("T yy =", conductivity[0][5][49])
    print("T zz =", conductivity[0][9][49])
    print("- - - - - - - - - - - - - - -")
    print("L xx =", conductivity[1][1][49])
    print("L yy =", conductivity[1][5][49])
    print("L zz =", conductivity[1][9][49])
    print("- - - - - - - - - - - - - - -")
    print("sum xx =", conductivity[2][1][49])
    print("sum yy =", conductivity[2][5][49])
    print("sum zz =", conductivity[2][9][49])

    window_diagonal = [0, 0]
    max_offdiagonal = 0
    for key in plot_list_conductivity_diagonal.keys():
        window_diagonal[0] = min(window_diagonal[0], is_T[key].min(), is_L[key].min())
        window_diagonal[1] = max(window_diagonal[1], is_T[key].max(), is_L[key].max())

    for key in plot_list_conductivity_offdiagonal.keys():
        max_offdiagonal = max(max_offdiagonal, np.abs(is_T[key].min()), np.abs(is_L[key].min()), np.abs(is_T[key].max()), np.abs(is_L[key].max()))

    window_diagonal = np.array([-window_diagonal[1]*0.02, window_diagonal[1]*1.1])
    window_offdiagonal = np.array([-max_offdiagonal, max_offdiagonal])*5

    window_L_offdiagonal = [0, 0]
    for valley in np.arange(1,4):
        i = valley
        for key in plot_list_conductivity_offdiagonal.keys():
            window_L_offdiagonal[0] = min(window_L_offdiagonal[0], conductivity_i[i][key].min())
            window_L_offdiagonal[1] = max(window_L_offdiagonal[1], conductivity_i[i][key].max())

    window_L_offdiagonal = np.array(window_L_offdiagonal) * 1.1

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
    for ax in axes.T[-1]:
        ax.set_xlabel("$\mu~[\mathrm{meV}]$")

    axes[0][0].set_ylabel("$\mathrm{conductivity~[/\Omega \, m]}$")
    axes[0][1].set_ylabel("$\mathrm{conductivity~[/\Omega \, m]}$")

    labels = ["$h$", "$e_1 + e_2 + e_3$", "$h + e_1 + e_2 + e_3$"]
    for n in np.arange(0,3):
        axes[n][0].text(0.07, 0.84, labels[n], bbox=bbox, transform=axes[n][0].transAxes)
        axes[n][1].text(0.07, 0.84, labels[n], bbox=bbox, transform=axes[n][1].transAxes)
        i = 1
        for key, val in plot_list_conductivity_diagonal.items():
            axes[n][0].set_ylim(window_diagonal)
            axes[n][0].scatter(x, conductivity[n][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
            i = i + 1

        for key, val in plot_list_conductivity_offdiagonal.items():
            axes[n][1].set_ylim(window_offdiagonal)
            axes[n][1].scatter(x, conductivity[n][key], marker=markers[i-1], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+val+"}$")
            i = i + 1

    for i in np.arange(col_num):
        axes[0][i].legend(labelspacing=0.2, handletextpad=0.2, markerscale=2)

#    plt.show()
    fig.align_ylabels()
    fig.tight_layout()
    plt.savefig(pdf+output+".pdf", bbox_inches='tight')
    plt.close()
# }}}

def plot_spin_conductivity_valleys(sign, titles, conductivity, conductivity_valley, output): # spin conductivity {{{
    conductivity_i = []
    conductivity_i_sum = []

    fig, axes = plt.subplots(1,3,figsize=(figure_size_unit*3, figure_size_unit))
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

#    print("mu =", conductivity[2][0][49])
#    print("- - - - - - - - - - - - - - -")
#    print("T xyz =", conductivity[0][ 6][49])
#    print("T yzx =", conductivity[0][16][49])
#    print("T zxy =", conductivity[0][20][49])
#    print("- - - - - - - - - - - - - - -")
#    print("L xyz =", conductivity[1][ 6][49])
#    print("L yzx =", conductivity[1][16][49])
#    print("L zxy =", conductivity[1][20][49])
    print("- - - - - - - - - - - - - - -")
    print(sign, " xyz =", conductivity[2][ 6][49])
    print(sign, " yzx =", conductivity[2][16][49])
    print(sign, " zxy =", conductivity[2][20][49])

#    window = []
#    for j, plot_list in enumerate(plot_lists):
#        tmp = [0, 0]
#        for n in np.arange(0,3):
#            for key in plot_list.keys():
#                tmp[0] = min(tmp[0], conductivity[n][key].min())
#                tmp[1] = max(tmp[1], conductivity[n][key].max())
#        tmp = np.array(tmp) * 1.1e0
#        window.append(tmp)
#
#    window[3] = window[2]
#    window[4] = window[2]

    maximun = 0e0
    for i in np.arange(conductivity.shape[1]):
        for n in np.arange(0,3):
            maximun = max([maximun, np.abs(conductivity[n][i].max()), np.abs(conductivity[n][i].min())])
    window = [-maximun*1.1e0, maximun*1.1e0]

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

    labels = ["$h$", "$e_1 + e_2 + e_3$", "$h + e_1 + e_2 + e_3$"]
    for n in np.arange(0,3):
        k = 0
        for j, plot_list in enumerate(plot_lists):
            i = 0
            axes[n][k].text(0.07, 0.87, labels[n], bbox=bbox, fontsize="large", transform=axes[n][k].transAxes)
            for key, val in plot_list.items():
#                axes[n][j].set_ylim(window[j])
                axes[n][j].set_ylim(window)
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

    fig, axes = plt.subplots(1,3,figsize=(figure_size_unit*3, figure_size_unit))
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

def plot_compare(magnetic, angular, output): # {{{
    figure_size_unit = 6

    plot_list1  = { 6:"xyz", 16:"yzx", 20:"zxy"}
    plot_list2  = { 1:"xxx", 13:"yyx",  5:"xyy", 11:"yxy"}
    plot_lists = [plot_list1, plot_list2]

    sigma_m, magnetic_conductivity = interpolate_conductivity(magnetic[0], magnetic[1])
    sigma_a, angular_conductivity  = interpolate_conductivity(angular[0], angular[1])
    conductivity = [magnetic_conductivity, angular_conductivity]

    col_num = 2
    row_num = 2
    fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num*1.2,figure_size_unit*row_num))
    axes = axes.T

    for ax_list in axes:
        for ax in ax_list:
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.grid()
            ax.set_axisbelow(True)
            ax.set_xlim(-float(cutoff), float(cutoff))
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
    for ax in axes.T[-1]:
        ax.set_xlabel("$\mu~[\mathrm{meV}]$")

    axes[0][0].set_ylabel(r"$(e/\mu_{\mathrm{B}}) \sigma_{\mathrm{m}}~\mathrm{[/\Omega \, m]}$")
    axes[0][1].set_ylabel(r"$(2e/\hbar) \sigma_{\mathrm{a}}~\mathrm{[/\Omega \, m]}$")
#    axes[0][0].set_ylabel("$\mathrm{spin\ magnetic\ conductivity~[/\Omega \, m]}$")
#    axes[0][1].set_ylabel("$\mathrm{spin\ angular\ conductivity~[/\Omega \, m]}$")

    window = [0, 0]
    for i in np.arange(2):
        maximun = 0e0
        for plot_list in plot_lists:
            for key, val in plot_list.items():
                maximun = max([maximun, np.abs(conductivity[i][key].max()), np.abs(conductivity[i][key].min())])
        window[i] = [-maximun*1.4e0, maximun*1.4e0]

    window[0] = [-4.5e6, 0.1e6]
    window[1] = [-2.0e4, 2.2e4]

    signs = ["\mathrm{m}", "\mathrm{a}"]
    loc = [["lower right", "upper right"], ["upper right", "center right"]]
    for n in np.arange(0,2):
        for ax in axes.T[n]:
            ax.set_ylim(window[n])

        sign = signs[n]
        sigma = conductivity[n]
        for j, plot_list in enumerate(plot_lists):
            i = 0
            for key, val in plot_list.items():
                force = val[0]
                flow  = val[1]
                spin  = val[2]
                axes.T[n][j].scatter(sigma[0], sigma[key], marker=markers[i], s=20, edgecolors=colors[i+1], facecolor='None', label="$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
                i += 1
            axes.T[n][j].legend(labelspacing=0.2, handletextpad=0.2, loc=loc[n][j], markerscale=2)

    axes.T[1][0].legend(labelspacing=0.2, loc=loc[1][0], bbox_to_anchor = (1, 0.95), markerscale=2)

    left, bottom, width, height = [0.13, 0.1, 0.4, 0.3]
    ax = axes.T[0][1].inset_axes([left, bottom, width, height])
    ticks0 = [-0.06, -0.04, -0.02]
    window0x = [ticks0[0]-0.002, ticks0[-1]+0.002]
    ax.grid()
    ax.set_xlim(window0x)
    ax.set_xticks(ticks0)
    ax.set_xticklabels(["$-60$", "$-40$", "$-20$"])
    ax.tick_params(axis="both", which='major', labelsize=18)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
    ax.yaxis.offsetText.set_fontsize(18)
    i = 0
    for key, val in plot_list2.items():
        ax.scatter(conductivity[0][0], conductivity[0][key], s=20, marker=markers[i], edgecolors=colors[i+1], facecolor='None')
        i += 1
    axes.T[0][1].indicate_inset_zoom(ax, edgecolor="black")

    rect1 = patches.Rectangle((0.02, 0.02), width=0.56, height=0.4, edgecolor=None, facecolor='w', alpha=0.8, transform=axes.T[0][1].transAxes)
    axes.T[0][1].add_patch(rect1)

#    fig.text(1, 1, str_damping, ha='right', va='bottom',transform=axes[-1][0].transAxes)

#    plt.show()
    fig.tight_layout()
    plt.savefig(png+output+".png", bbox_inches='tight', dpi=300)
    plt.savefig(pdf+output+".pdf")
    plt.close()

    conductivity = [sigma_m, sigma_a]
    window[0] = [-4.4e6, 0.1e6]
    window[1] = [-1.2e4, 2.3e4]

    col_num = 5
    row_num = 2
    fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
    axes = axes.T

    for ax_list in axes:
        for ax in ax_list:
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.grid()
            ax.set_axisbelow(True)
            ax.set_xlim(-float(cutoff), float(cutoff))
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
    for ax in axes.T[-1]:
        ax.set_xlabel("$\mu~[\mathrm{meV}]$")

    axes[0][0].set_ylabel(r"$(e/\mu_{\mathrm{B}}) \sigma_{\mathrm{m}}~\mathrm{[/\Omega \, m]}$")
    axes[0][1].set_ylabel(r"$(2e/\hbar) \sigma_{\mathrm{a}}~\mathrm{[/\Omega \, m]}$")

    signs = ["\mathrm{m}", "\mathrm{a}"]
    loc = [["lower right", "upper right"], ["upper right", "center right"]]
    for n in np.arange(0,2):
        for ax in axes.T[n]:
            ax.set_ylim(window[n])

        sign = signs[n]
        sigma = conductivity[n]
        for j in np.arange(5):
            i = 0
            for key, val in plot_list1.items():
                force = val[0]
                flow  = val[1]
                spin  = val[2]
                axes.T[n][j].scatter(sigma[j][0], sigma[j][key], marker=markers[i], s=20, edgecolors=colors[i+1], facecolor='None', label="$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
                i += 1
#            axes.T[n][j].legend(labelspacing=0.2, handletextpad=0.2, loc=loc[n][j], markerscale=2)

    axes.T[0][0].legend(labelspacing=0.2, loc=loc[0][0], markerscale=2)
    axes.T[1][0].legend(labelspacing=0.2, loc=loc[1][0], bbox_to_anchor = (1, 0.95), markerscale=2)

    labels = ["$h$", "$e_1$", "$e_2$", "$e_3$", "$e_1+e_2+e_3$"]
    for n in np.arange(5):
        axes.T[0][n].text(0.07, 0.8, labels[n], bbox=bbox, transform=axes.T[0][n].transAxes, alpha=0.8)
        axes.T[1][n].text(0.07, 0.8, labels[n], bbox=bbox, transform=axes.T[1][n].transAxes, alpha=0.8)

    fig.tight_layout()
    plt.savefig(png+output+"_valley.png", bbox_inches='tight', dpi=300)
    plt.savefig(pdf+output+"_valley.pdf")
    plt.close()
# }}}

#
#def plot_compare(magnetic, angular, output): # {{{
#    figure_size_unit = 6
#    plot_list1  = { 6:"xyz", 16:"yzx", 20:"zxy"}
#    plot_list2  = { 1:"xxx", 13:"yyx",  5:"xyy", 11:"yxy"}
##    plot_lists = [plot_list1, plot_list2]
#    plot_lists = [plot_list1]
#
#    sigma_m, magnetic_conductivity = interpolate_conductivity(magnetic[0], magnetic[1])
#    sigma_a, angular_conductivity  = interpolate_conductivity(angular[0], angular[1])
#    conductivity = [magnetic_conductivity, angular_conductivity]
#
#    col_num = 2
#    row_num = 1
#    fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
#
#    for ax in axes:
#        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
#        ax.grid()
#        ax.set_axisbelow(True)
#        ax.set_xlim(-float(cutoff), float(cutoff))
#        ax.set_xticks(xticks)
#        ax.set_xticklabels(xticklabels)
#        ax.set_xlabel("$\mu~[\mathrm{meV}]$")
#
#    axes[0].set_ylabel(r"$(e/\mu_{\mathrm{B}}) \sigma_{\mathrm{m}}~\mathrm{[/\Omega \, m]}$")
#    axes[1].set_ylabel(r"$(2e/\hbar) \sigma_{\mathrm{a}}~\mathrm{[/\Omega \, m]}$")
#
#    window = [0, 0]
#    for i in np.arange(2):
#        maximun = 0e0
#        for plot_list in plot_lists:
#            for key, val in plot_list.items():
#                maximun = max([maximun, np.abs(conductivity[i][key].max()), np.abs(conductivity[i][key].min())])
#        window[i] = [-maximun*1.1e0, maximun*1.1e0]
#    window[0][1] = 0
#    window[1][0] = -1e4
#
#    signs = ["\mathrm{m}", "\mathrm{a}"]
#    loc = [["lower right", "upper right"], ["upper right", "center right"]]
#    for n in np.arange(0,2):
#        axes[n].set_ylim(window[n])
#
#        sign = signs[n]
#        sigma = conductivity[n]
#        for j, plot_list in enumerate(plot_lists):
#            i = 0
#            for key, val in plot_list.items():
#                force = val[0]
#                flow  = val[1]
#                spin  = val[2]
#                axes[n].scatter(sigma[0], sigma[key], marker=markers[i], s=20, edgecolors=colors[i+1], facecolor='None', label="$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
#                i += 1
#            axes[n].legend(labelspacing=0.2, handletextpad=0.2, loc=loc[n][j])
#
#    axes[1].legend(labelspacing=0.2, loc=loc[1][0], bbox_to_anchor = (1, 0.95))
#
##    fig.text(1, 1, str_damping, ha='right', va='bottom',transform=axes[-1][0].transAxes)
#
##    plt.show()
#    fig.tight_layout()
#    plt.savefig(png+output+".png", bbox_inches='tight', dpi=300)
#    plt.savefig(pdf+output+".pdf")
#    plt.close()
## }}}
#
def main():
    parameter = "T"+str(bandsT)+"bands_L"+str(bandsL)+"bands"+"_"+param+"eV"
#    label_damping_T = "_gamma0.5e-4"
#    label_damping_T = "_gamma1e-4"
#    label_damping_T = "_gamma0.8e-5"
    label_damping_T = "_gamma0.7e-4"
    label_damping_L = "_gamma1.6e-4"

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
    for i in np.arange(3):
        valley = key_L[i]
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

#    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'dos_eps'+label+'.csv'
    d = pd.read_csv(readfile,header=None).values
    dos.append(d)
    e_max.append(d[:,0].max())
    e_min.append(d[:,0].min())

    for i in np.arange(3):
        valley = key_L[i]
#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'dos_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=None).values
        e_max.append(d4[:,0].max())
        e_min.append(d4[:,0].min())

#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
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

#    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity.append(d)
    conductivity_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        conductivity.append(d4)

#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values
        conductivity.append(d6)

        conductivity_valley.append(np.append(d4, d6, axis=0))

    plot_conductivity_valleys(titles,conductivity,conductivity_valley,"conductivity"+parameter)
    plot_conductivity_valley_resolved(titles,conductivity,conductivity_valley,"conductivity"+parameter)
    plot_conductivity(conductivity,conductivity_valley,"conductivity_tensor"+parameter)
# }}}

# spin magnetic conductivity1 {{{
    conductivity1 = []
    conductivity1_valley = []

#    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
    print(readfile)
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity1.append(d)
    conductivity1_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-magnetic-conductivity1_eps'+label+'.csv'
        print(readfile)
        d4 = pd.read_csv(readfile,header=0).values
        conductivity1.append(d4)

#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
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

#    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
    d = pd.read_csv(readfile,header=0).values
    conductivity2.append(d)
    conductivity2_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        conductivity2.append(d4)

#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-magnetic-conductivity2_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values

        val = np.copy(d4[-1, :])
        val[0] = 0e0
        d6 = d6 + val

        conductivity2.append(d6)
        conductivity2_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("m", titles,conductivity2,conductivity2_valley,"spin_magnetic_conductivity2_"+parameter)
# }}}

    magnetic = [conductivity1_valley, conductivity2_valley]

# spin magnetic conductivity 1 + 2 {{{
    plot_spin_conductivity("m", conductivity1_valley, conductivity2_valley, "spin_magnetic_conductivity_tensor_"+parameter)
    plot_total_spin_conductivity("m", titles, conductivity1_valley, conductivity2_valley, "total_spin_magnetic_conductivity_"+parameter)
# }}}

# spin angular conductivity1 {{{
    conductivity1 = []
    conductivity1_valley = []

#    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity1.append(d)
    conductivity1_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        conductivity1.append(d4)

#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-angular-conductivity1_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values
        conductivity1.append(d6)

        conductivity1_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("a", titles,conductivity1,conductivity1_valley,"spin_angular_conductivity1_"+parameter)
# }}}

# spin angular conductivity2 {{{
    conductivity2 = []
    conductivity2_valley = []

#    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+'/'
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
    d = pd.read_csv(readfile,header=0).values
    conductivity2.append(d)
    conductivity2_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values
        conductivity2.append(d4)

#        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+'/'
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-angular-conductivity2_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values

        val = np.copy(d4[-1, :])
        val[0] = 0e0
        d6 = d6 + val

        conductivity2.append(d6)
        conductivity2_valley.append(np.append(d4, d6, axis=0))

    plot_spin_conductivity_valleys("a", titles,conductivity2,conductivity2_valley,"spin_angular_conductivity2_"+parameter)
# }}}

    angular = [conductivity1_valley, conductivity2_valley]

# spin angular conductivity 1 + 2 {{{
    plot_spin_conductivity("a", conductivity1_valley, conductivity2_valley, "spin_angular_conductivity_tensor_"+parameter)
    plot_total_spin_conductivity("a", titles, conductivity1_valley, conductivity2_valley, "total_spin_angular_conductivity_"+parameter)
# }}}

    plot_compare(magnetic, angular, "compare")

if __name__ == "__main__":
    main()
