import os
import numpy as np
import pandas as pd

from scipy.interpolate import griddata
from scipy import interpolate

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

from decimal import Decimal, Context

markers = ["o", "D", "v", "^", "<", ">", "s", "p", "*"]
colors =['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo', 'k', 'b', 'g', 'r']

matplotlib.rcParams["font.family"] = 'Times New Roman'
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams["font.size"] = 16
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

def decimal_normalize(f, r='0.001'):
    """数値fの誤差をrの位で四捨五入する"""
    f = Decimal(str(f)).quantize(Decimal(r))
    """数値fの小数点以下を正規化し、文字列で返す"""
    def _remove_exponent(d):
        return d.quantize(Decimal(1)) if d == d.to_integral() else d.normalize()
    a = Decimal.normalize(Decimal(str(f)))
    b = _remove_exponent(a)
    return str(b)


figure_size_unit = 8
bbox = {
    "facecolor" : "white",
    "edgecolor" : "black",
    "boxstyle" : "round, pad=0.6",
    "linewidth" : 2
}


base = "../"
dat_dir = base+"dat/"
pdf_base_dir = base+"fig/pdf/"

bandsT = 12
bandsL = 4
lower_band_L = 0
upper_band_L = 2

cutoff = 0.1e0
xticks = np.linspace(-cutoff, cutoff, 8+1)
xticklabels = ["$"+decimal_normalize(tick)+"$" for tick in xticks]

axises = {"x":0, "y":1, "z":2}

plot_list_conductivity_diagonal = {1:"xx", 5:"yy", 9:"zz"}
plot_list_conductivity_offdiagonal = {2:"xy", 3:"xz", 4:"yx", 6:"yz", 7:"zx", 8:"zy"}

plot_list_magnetic1 = {6:"xyz", 16:"yzx", 20:"zxy"}
plot_list_magnetic2 = {12:"yxz", 8:"xzy", 22:"zyx"}
plot_list_magnetic3 = {1:"xxx", 13:"yyx", 5:"xyy", 11:"yxy"}
plot_list_magnetic4a = {2:"xxy", 3:"xxz", 4:"xyx", 7:"xzx", 9:"xzz", 10:"yxx", 14:"yyy", 15:"yyz", 17:"yzy"}
plot_list_magnetic4b = {18:"yzz", 19:"zxx", 21:"zxz", 23:"zyy", 24:"zyz", 25:"zzx", 26:"zzy", 27:"zzz"}
plot_list_spin_magnetic_conductivity = [ plot_list_magnetic1, plot_list_magnetic2, plot_list_magnetic3, plot_list_magnetic4a, plot_list_magnetic4b ]

plot_list_sigma = [{5:"yy", 9:"zz"}, {12:"yxz", 8:"xzy"}]

def dos_T(params): # dos T {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    data = dat_dir+'T_'+str(bandsT)+'bands/band_index4/'+label_cutoff+"/"
    readfile = data+'dos_eps'+label_eps+'.csv'
    x = pd.read_csv(readfile,header=None)[0].values
    d = pd.read_csv(readfile,header=None)[1].values

    if params['output']:
        fig, ax = plt.subplots(1,1,figsize=(figure_size_unit,figure_size_unit))
        ax.set_xlim([-cutoff, cutoff])
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        ax.set_ylabel("$\mathrm{DOS}~[/\mathrm{eV m^3}]$")
        ax.grid()
        ax.plot(x, d, label="T", c="red")
        ax.legend()
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
#        plt.show()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

        fig.tight_layout()
        plt.savefig(pdf_dir+"dos_T"+str(bandsT)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return x, d
# }}}
def dos_L(params): # dos L {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])

    dos_lower = []
    dos_upper = []
    for valley in np.arange(1,4):
        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+label_cutoff+'/'
        readfile = data+'dos_eps'+label_eps+'.csv'
        x_lower = pd.read_csv(readfile,header=None)[0].values
        d_lower = pd.read_csv(readfile,header=None)[1].values
        dos_lower.append(d_lower)

        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/cutoff'+str(params['cutoff'])+'eV/'
        readfile = data+'dos_eps'+label_eps+'.csv'
        x_upper = pd.read_csv(readfile,header=None)[0].values
        d_upper = pd.read_csv(readfile,header=None)[1].values
        dos_upper.append(d_upper)

    if params['output']:
        fig, ax = plt.subplots(1,1,figsize=(figure_size_unit,figure_size_unit))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.set_xlim([-cutoff, cutoff])
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        ax.set_ylabel("$\mathrm{DOS}~[/\mathrm{eV m^3}]$")
        ax.grid()

        for valley in np.arange(1,4):
            i = valley - 1
            ax.plot(x_lower, dos_lower[i], label="L"+str(valley)+"-"+str(lower_band_L))
            ax.plot(x_upper, dos_upper[i], label="L"+str(valley)+"-"+str(upper_band_L))

        ax.legend()
#        plt.show()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

        fig.tight_layout()
        plt.savefig(pdf_dir+"dos_L"+str(bandsL)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    dos_lower = dos_lower[0] + dos_lower[1] + dos_lower[2]
    dos_upper = dos_upper[0] + dos_upper[1] + dos_upper[2]

    x = np.concatenate([x_lower, x_upper])
    dos = np.concatenate([dos_lower, dos_upper])

    return x, dos
# }}}
def dos(params): # individual + total dos {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])

    d_T = dos_T(params)
    d_L = dos_L(params)

    n = params['slice']
    e = np.linspace(-params['cutoff'], params['cutoff'], n+1)
    id_T = interpolate.interp1d(d_T[0], d_T[1])
    id_L = interpolate.interp1d(d_L[0], d_L[1])

    dos = np.array([id_T(e), id_L(e)])
    dos_total = np.sum(dos, axis=0)
    dos_total = np.stack([e, dos_total])

    if params['output']:
        fig, ax = plt.subplots(1,1,figsize=(figure_size_unit,figure_size_unit))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        ax.set_ylabel("$\mathrm{DOS}~[/\mathrm{eV \, m^3}]$")
        ax.grid()
        ax.set_xlim([-cutoff, cutoff])
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)

        ax.scatter(e, dos[0], label="T", marker=markers[0], s=70, edgecolors=colors[4], facecolor='None')
        ax.scatter(e, dos[1], label="L", marker=markers[1], s=70, edgecolors=colors[5], facecolor='None')
        ax.scatter(e, dos_total[1], label="Total", marker=markers[2], s=70, edgecolors=colors[6], facecolor='None')

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

        ax.legend()
        fig.tight_layout()
        plt.savefig(pdf_dir+"dos_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return dos_total
# }}}

def electric_conductivity_T(params):# electric conductivity T {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    data = dat_dir+'T_'+str(bandsT)+'bands/band_index4/'+label_cutoff+"/"
    readfile = data+'mu-dependence/conductivity_eps'+label_eps+'.csv'
    df = pd.read_csv(readfile,header=0)
    d = df.values

    titles = df.columns.values
    for key, val in plot_list_conductivity_diagonal.items():
        if titles[key].strip() != val:
            print(titles[key], val)
    for key, val in plot_list_conductivity_offdiagonal.items():
        if titles[key].strip() != val:
            print(titles[key], val)

    if params['output']:
        window = [-0.02e0*d[:,1:].max(), d[:,1:].max()*1.1e0]

        fig, axes = plt.subplots(1,2,figsize=(figure_size_unit*2,figure_size_unit))
        axes = axes.flatten()
        axes[0].set_ylabel("$\mathrm{conductivity~[/\Omega \, m]}$")
        for ax in axes:
            ax.set_xlim([-cutoff, cutoff])
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
            ax.grid()
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_ylim(window)

        axes[0].set_title("diagonal")
        i = 0
        for key, val in plot_list_conductivity_diagonal.items():
            axes[0].scatter(d[:,0], d[:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
            i = i + 1

        axes[1].set_title("off-diagonal")
        for key, val in plot_list_conductivity_offdiagonal.items():
            axes[1].scatter(d[:,0], d[:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
            i = i + 1

        ax_pos = axes[1].get_position()
        fig.text(ax_pos.x1 - 0.02, ax_pos.y1 + 0.1, str_damping)

#        axes[0].text(-cutoff*0.9, window[1]*0.93, "T", bbox=bbox, fontsize="large")

        for ax in axes:
            ax.legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"conductivity_T"+str(bandsT)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        #plt.show()
        plt.close()

    return d
# }}}
def electric_conductivity_L(params):# electric conductivity L {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    conductivity = []
    for valley in np.arange(1,4):
        col = valley - 1

        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+label_cutoff+'/'
        readfile = data+'mu-dependence/conductivity_eps'+label_eps+'.csv'
        df = pd.read_csv(readfile,header=0)
        d_lower = df.values

        titles = df.columns.values
        for key, val in plot_list_conductivity_diagonal.items():
            if titles[key].strip() != val:
                print(titles[key], val)
        for key, val in plot_list_conductivity_offdiagonal.items():
            if titles[key].strip() != val:
                print(titles[key], val)

        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/cutoff'+str(params['cutoff'])+'eV/'
        readfile = data+'mu-dependence/conductivity_eps'+label_eps+'.csv'
        df = pd.read_csv(readfile,header=0)
        d_upper = df.values

        titles = df.columns.values
        for key, val in plot_list_conductivity_diagonal.items():
            if titles[key].strip() != val:
                print(titles[key], val)
        for key, val in plot_list_conductivity_offdiagonal.items():
            if titles[key].strip() != val:
                print(titles[key], val)

        sigma = np.concatenate([d_lower, d_upper])
#        sigma[np.abs(sigma[:,0]) > cutoff] = np.nan
#        sigma = sigma[~np.isnan(sigma).any(axis=1)]
        conductivity.append(sigma)

    if params['output']:
        window_diagonal = [0, 0]
        window_offdiagonal = [0, 0]

        for valley in np.arange(1,4):
            i = valley - 1
            for key in plot_list_conductivity_diagonal.keys():
                window_diagonal[0] = min(window_diagonal[0], conductivity[i][:,key].min())
                window_diagonal[1] = max(window_diagonal[1], conductivity[i][:,key].max())

            for key in plot_list_conductivity_offdiagonal.keys():
                window_offdiagonal[0] = min(window_offdiagonal[0], conductivity[i][:,key].min())
                window_offdiagonal[1] = max(window_offdiagonal[1], conductivity[i][:,key].max())

        window_diagonal = np.array([-window_diagonal[1]*0.02, window_diagonal[1]*1.1])
        window_offdiagonal = np.array(window_offdiagonal)*1.1

        col_num = 3
        row_num = 2
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.T

        for ax_col in axes:
            for ax in ax_col:
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[1]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")

        axes[0][0].set_ylabel("longitudinal conductivity$\mathrm{~[/\Omega \, m]}$")
        axes[0][1].set_ylabel("transverse conductivity  $\mathrm{~[/\Omega \, m]}$")

        for valley in np.arange(1,4):
            col = valley - 1
            i = 0
            for key, val in plot_list_conductivity_diagonal.items():
#                axes[col][0].set_yscale("log")
                axes[col][0].set_ylim(window_diagonal)
                axes[col][0].scatter(conductivity[col][:,0], conductivity[col][:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1

            for key, val in plot_list_conductivity_offdiagonal.items():
                axes[col][1].set_ylim(window_offdiagonal)
                axes[col][1].scatter(conductivity[col][:,0], conductivity[col][:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1

            axes[col][0].text(0.05, 0.9, "L"+str(valley), bbox=bbox, fontsize="large", transform=axes[col][0].transAxes)

        ax_pos = axes[2][0].get_position()
        fig.text(ax_pos.x1 - 0.01, ax_pos.y1 + 0.11, str_damping)

        axes[0][0].legend()
        axes[0][1].legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"conductivity_L"+str(bandsL)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    e = conductivity[0][:,0]
    conductivity = conductivity[0] + conductivity[1] + conductivity[2]
    conductivity[:,0] = e

    return conductivity
# }}}
def electric_conductivity(params): # individual + total electric conductivity {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    sigma_T = electric_conductivity_T(params)
    sigma_L = electric_conductivity_L(params)

    n = params['slice']
    e = np.linspace(-params['cutoff'], params['cutoff'], n+1)
    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(e) for i in np.arange(1,10)])
    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(e) for i in np.arange(1,10)])

    conductivity = np.array([is_T, is_L])
    conductivity_total = np.sum(conductivity, axis=0)
    is_T = np.concatenate([np.array([e]), is_T])
    is_L = np.concatenate([np.array([e]), is_L])
    conductivity_total = np.concatenate([np.array([e]), conductivity_total])
    conductivity = np.array([is_T, is_L, conductivity_total])

    if params['output']:
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
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[1]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")

        axes[0][0].set_ylabel("longitudinal conductivity$\mathrm{~[/\Omega \, m]}$")
        axes[0][1].set_ylabel("transverse conductivity  $\mathrm{~[/\Omega \, m]}$")

        labels = ["T", "L", "Total"]
        for n in np.arange(0,3):
            axes[n][0].text(-cutoff*0.9, window_diagonal[1]*0.93, labels[n], bbox=bbox, fontsize="large")
            i = 0
            for key, val in plot_list_conductivity_diagonal.items():
                axes[n][0].set_ylim(window_diagonal)
                axes[n][0].scatter(e, conductivity[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1

            for key, val in plot_list_conductivity_offdiagonal.items():
                axes[n][1].set_ylim(window_offdiagonal)
                axes[n][1].scatter(e, conductivity[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1

        axes[0][0].legend()
        axes[0][1].legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"conductivity_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return conductivity
# }}}

def spin_magnetic_conductivity_T(params, N):# spin magnetic conductivity T {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    data = dat_dir+'T_'+str(bandsT)+'bands/band_index4/'+label_cutoff+"/"
    readfile = data+'mu-dependence/spin-magnetic-conductivity'+str(N)+'_eps'+label_eps+'.csv'
    df = pd.read_csv(readfile,header=0)
    d = df.values

    titles = df.columns.values
    for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
        for key, val in plot_list.items():
            if titles[key].strip() != val:
                print(titles[key].strip(), val)

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            tmp = [0, 0]
            for key in plot_list.keys():
                tmp[0] = min(tmp[0], d[:,key].min())
                tmp[1] = max(tmp[1], d[:,key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        col_num = 1
        row_num = 5
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.flatten()
        for ax in axes:
            ax.set_xlim([-cutoff, cutoff])
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
            ax.grid()
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            i = 0
            for key, val in plot_list.items():
                axes[j].set_ylim(window[j])
                axes[j].scatter(d[:,0], d[:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1

        ax_pos = axes[1].get_position()
        fig.text(ax_pos.x1 - 0.02, ax_pos.y1 + 0.1, str_damping)

#        axes[0].text(-cutoff*0.9, window[1]*0.93, "T", bbox=bbox, fontsize="large")

        for ax in axes:
            ax.legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_magnetic_conductivity"+str(N)+"_T"+str(bandsT)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        #plt.show()
        plt.close()

    return d
# }}}
def spin_magnetic_conductivity_L(params, N):# spin magnetic conductivity L {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    conductivity = []
    for valley in np.arange(1,4):
        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+label_cutoff+'/'
        readfile = data+'mu-dependence/spin-magnetic-conductivity'+str(N)+'_eps'+label_eps+'.csv'
        df = pd.read_csv(readfile,header=0)
        d_lower = df.values

        titles = df.columns.values
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            for key, val in plot_list.items():
                if titles[key].strip() != val:
                    print(titles[key].strip(), val)

        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/cutoff'+str(params['cutoff'])+'eV/'
        readfile = data+'mu-dependence/spin-magnetic-conductivity'+str(N)+'_eps'+label_eps+'.csv'
        df = pd.read_csv(readfile,header=0)
        d_upper = df.values
        if N == 2:
            val = np.copy(d_lower[-1, :])
            val[0] = 0e0
            d_upper = d_upper + val

        titles = df.columns.values
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            for key, val in plot_list.items():
                if titles[key].strip() != val:
                    print(titles[key].strip(), val)

        sigma = np.concatenate([d_lower, d_upper])
        conductivity.append(sigma)

    if params['output']:
        col_num = 3
        row_num = 5
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.T

        for ax_list in axes:
            for ax in ax_list:
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[3]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            tmp = [0, 0]
            for valley in np.arange(1,4):
                for key in plot_list.keys():
                    tmp[0] = min(tmp[0], conductivity[valley-1][:,key].min())
                    tmp[1] = max(tmp[1], conductivity[valley-1][:,key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        for valley in np.arange(1,4):
            k = valley - 1
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[k][j].set_ylim(window[j])
                    axes[k][j].scatter(conductivity[k][:,0], conductivity[k][:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[k][0].text(0.05, 0.9, "L"+str(valley), bbox=bbox, fontsize="large", transform=axes[k][0].transAxes)

        ax_pos = axes[2][0].get_position()
        fig.text(0.8, 1.01, str_damping, transform=axes[2][0].transAxes)

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_magnetic_conductivity"+str(N)+"_L"+str(bandsL)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    e = conductivity[0][:,0]
    conductivity = conductivity[0] + conductivity[1] + conductivity[2]
    conductivity[:,0] = e

    return conductivity
# }}}
def spin_magnetic_conductivity(params, N): # individual spin magnetic conductivity {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    sigma_T = spin_magnetic_conductivity_T(params, N)
    sigma_L = spin_magnetic_conductivity_L(params, N)

    n = 100
    e = np.linspace(-params['cutoff'], params['cutoff'], n+1)
    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(e) for i in np.arange(1,28)])
    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(e) for i in np.arange(1,28)])

    conductivity = np.array([is_T, is_L])
    conductivity_total = np.sum(conductivity, axis=0)
    is_T = np.concatenate([np.array([e]), is_T])
    is_L = np.concatenate([np.array([e]), is_L])
    conductivity_total = np.concatenate([np.array([e]), conductivity_total])
    conductivity = np.array([is_T, is_L, conductivity_total])

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
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
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[3]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        labels = ["T", "L", "Total"]
        for n in np.arange(0,3):
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[n][j].set_ylim(window[j])
                    axes[n][j].scatter(e, conductivity[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[n][0].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][0].transAxes)

        axes[0][0].legend()
        axes[0][1].legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_magnetic_conductivity"+str(N)+"_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return conductivity

# }}}
def total_spin_magnetic_conductivity(params): # total spin magnetic conductivity {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    sigma1 = spin_magnetic_conductivity(params, 1)
    sigma2 = spin_magnetic_conductivity(params, 2)

    n = params['slice']
    e = np.linspace(-params['cutoff'], params['cutoff'], n+1)

    sigma = np.sum([sigma1, sigma2], axis=0)
    sigma[:,0] = e
    conductivity = sigma

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            tmp = [0, 0]
            for key in plot_list.keys():
                tmp[0] = min(tmp[0], sigma[2][key].min())
                tmp[1] = max(tmp[1], sigma[2][key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        window[2] = window[2]*5e0
        window[3] = window[3]*10e0
        window[4] = window[4]*10e0

        col_num = 3
        row_num = 5
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.T

        for ax_list in axes:
            for ax in ax_list:
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[3]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        labels = ["T", "L", "Total"]
        for n in np.arange(0,3):
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[n][j].set_ylim(window[j])
                    axes[n][j].scatter(e, sigma[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[n][0].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][0].transAxes)

        axes[0][0].legend()
        axes[0][1].legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_magnetic_conductivity_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return conductivity
# }}}

def spin_angular_conductivity_T(params, N):# spin magnetic conductivity T {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    data = dat_dir+'T_'+str(bandsT)+'bands/band_index4/'+label_cutoff+"/"
    readfile = data+'mu-dependence/spin-angular-conductivity'+str(N)+'_eps'+label_eps+'.csv'
    df = pd.read_csv(readfile,header=0)
    d = df.values

    titles = df.columns.values
    for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
        for key, val in plot_list.items():
            if titles[key].strip() != val:
                print(titles[key].strip(), val)

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            tmp = [0, 0]
            for key in plot_list.keys():
                tmp[0] = min(tmp[0], d[:,key].min())
                tmp[1] = max(tmp[1], d[:,key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        col_num = 1
        row_num = 5
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.flatten()
        for ax in axes:
            ax.set_xlim([-cutoff, cutoff])
            ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
            ax.grid()
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels)
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            i = 0
            for key, val in plot_list.items():
                axes[j].set_ylim(window[j])
                axes[j].scatter(d[:,0], d[:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1

        ax_pos = axes[1].get_position()
        fig.text(ax_pos.x1 - 0.02, ax_pos.y1 + 0.1, str_damping)

#        axes[0].text(-cutoff*0.9, window[1]*0.93, "T", bbox=bbox, fontsize="large")

        for ax in axes:
            ax.legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_angular_conductivity"+str(N)+"_T"+str(bandsT)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        #plt.show()
        plt.close()

    return d
# }}}
def spin_angular_conductivity_L(params, N):# spin magnetic conductivity L {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    conductivity = []
    for valley in np.arange(1,4):
        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+label_cutoff+'/'
        readfile = data+'mu-dependence/spin-angular-conductivity'+str(N)+'_eps'+label_eps+'.csv'
        df = pd.read_csv(readfile,header=0)
        d_lower = df.values

        titles = df.columns.values
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            for key, val in plot_list.items():
                if titles[key].strip() != val:
                    print(titles[key].strip(), val)

        data = dat_dir+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/cutoff'+str(params['cutoff'])+'eV/'
        readfile = data+'mu-dependence/spin-angular-conductivity'+str(N)+'_eps'+label_eps+'.csv'
        df = pd.read_csv(readfile,header=0)
        d_upper = df.values
        if N == 2:
            val = np.copy(d_lower[-1, :])
            val[0] = 0e0
            d_upper = d_upper + val

        titles = df.columns.values
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            for key, val in plot_list.items():
                if titles[key].strip() != val:
                    print(titles[key].strip(), val)

        sigma = np.concatenate([d_lower, d_upper])
        conductivity.append(sigma)

    if params['output']:
        col_num = 3
        row_num = 5
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.T

        for ax_list in axes:
            for ax in ax_list:
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[3]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            tmp = [0, 0]
            for valley in np.arange(1,4):
                for key in plot_list.keys():
                    tmp[0] = min(tmp[0], conductivity[valley-1][:,key].min())
                    tmp[1] = max(tmp[1], conductivity[valley-1][:,key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        for valley in np.arange(1,4):
            k = valley - 1
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[k][j].set_ylim(window[j])
                    axes[k][j].scatter(conductivity[k][:,0], conductivity[k][:,key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[k][0].text(0.05, 0.9, "L"+str(valley), bbox=bbox, fontsize="large", transform=axes[k][0].transAxes)

        ax_pos = axes[2][0].get_position()
        fig.text(0.8, 1.01, str_damping, transform=axes[2][0].transAxes)

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_angular_conductivity"+str(N)+"_L"+str(bandsL)+"bands_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    e = conductivity[0][:,0]
    conductivity = conductivity[0] + conductivity[1] + conductivity[2]
    conductivity[:,0] = e

    return conductivity
# }}}
def spin_angular_conductivity(params, N): # individual spin magnetic conductivity {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    sigma_T = spin_magnetic_conductivity_T(params, N)
    sigma_L = spin_magnetic_conductivity_L(params, N)

    n = 100
    e = np.linspace(-params['cutoff'], params['cutoff'], n+1)
    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(e) for i in np.arange(1,28)])
    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(e) for i in np.arange(1,28)])

    conductivity = np.array([is_T, is_L])
    conductivity_total = np.sum(conductivity, axis=0)
    is_T = np.concatenate([np.array([e]), is_T])
    is_L = np.concatenate([np.array([e]), is_L])
    conductivity_total = np.concatenate([np.array([e]), conductivity_total])
    conductivity = np.array([is_T, is_L, conductivity_total])

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
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
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[3]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        labels = ["T", "L", "Total"]
        for n in np.arange(0,3):
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[n][j].set_ylim(window[j])
                    axes[n][j].scatter(e, conductivity[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[n][0].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][0].transAxes)

        axes[0][0].legend()
        axes[0][1].legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_angular_conductivity"+str(N)+"_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return conductivity

# }}}
def total_spin_angular_conductivity(params): # total spin magnetic conductivity {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    sigma1 = spin_magnetic_conductivity(params, 1)
    sigma2 = spin_magnetic_conductivity(params, 2)

    n = params['slice']
    e = np.linspace(-params['cutoff'], params['cutoff'], n+1)

    sigma = np.sum([sigma1, sigma2], axis=0)
    sigma[:,0] = e
    conductivity = sigma

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
            tmp = [0, 0]
            for key in plot_list.keys():
                tmp[0] = min(tmp[0], sigma[2][key].min())
                tmp[1] = max(tmp[1], sigma[2][key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        window[2] = window[2]*5e0
        window[3] = window[3]*10e0
        window[4] = window[4]*10e0

        col_num = 3
        row_num = 5
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.T

        for ax_list in axes:
            for ax in ax_list:
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[3]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        labels = ["T", "L", "Total"]
        for n in np.arange(0,3):
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[n][j].set_ylim(window[j])
                    axes[n][j].scatter(e, sigma[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[n][0].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][0].transAxes)

        axes[0][0].legend()
        axes[0][1].legend()

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"spin_angular_conductivity_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()

    return conductivity
# }}}

def plot_total(params): # {{{
    label_cutoff = 'cutoff'+str(params['cutoff'])+'eV'
    label_eps = "{:.6f}".format(params['eps'])
    str_damping = "$\gamma = "+decimal_normalize(params['eps']*1e3)+"~\mathrm{meV}$"

    sigma_e = electric_conductivity(params)
    sigma_m = total_spin_magnetic_conductivity(params)

    sigma = [sigma_e, sigma_m]

    if params['output']:
        window = [0, 0]
        for j, plot_list in enumerate(plot_list_sigma):
            for key in plot_list.keys():
                window[0] = min(window[0], sigma[j][2][key].min())
                window[1] = max(window[1], sigma[j][2][key].max())
        window = np.array(window) * 1.1e0

        col_num = 1
        row_num = 1
        fig, ax = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))

        ax.set_xlim([-cutoff, cutoff])
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
        ax.minorticks_on()
        ax.grid()
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels)
        ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        ax.set_ylabel("conductivity$\mathrm{~[/\Omega \, m]}$")
        ax.set_ylim(window)

        i = 0
        for j, plot_list in enumerate(plot_list_sigma):
            for key, val in plot_list.items():
                ax.scatter(sigma[j][2][0], sigma[j][2][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                i = i + 1
#        axes[n][0].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][0].transAxes)

        pdf_dir = pdf_base_dir + label_cutoff + "/"
        os.makedirs(pdf_dir, exist_ok=True)

        ax.legend()
#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"total_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()
# }}}

def cutoff_dependences(params): # {{{
    d = []
    sigma_e = []
    sigma_m = []
    for c in params['cutoffs']:
        params['cutoff'] = c

        d.append(dos(params))
        sigma_e.append(electric_conductivity(params))
        sigma_m.append(total_spin_magnetic_conductivity(params))

    sigma = [sigma_e, sigma_m]

    exit()

    if params['output']:
        window = []
        for j, plot_list in enumerate(plot_list_sigma):
            tmp = [0, 0]
            for key in plot_list.keys():
                tmp[0] = min(tmp[0], sigma[j][key].min())
                tmp[1] = max(tmp[1], sigma[j][key].max())
            tmp = np.array(tmp) * 1.1e0
            window.append(tmp)

        col_num = 3
        row_num = 1
        fig, axes = plt.subplots(row_num,col_num,figsize=(figure_size_unit*col_num,figure_size_unit*row_num))
        axes = axes.T

        for ax_list in axes:
            for ax in ax_list:
                ax.set_xlim([-cutoff, cutoff])
                ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
                ax.grid()
                ax.set_xticks(xticks)
                ax.set_xticklabels(xticklabels)
        for ax in axes.T[0]:
            ax.set_xlabel("$\mu~[\mathrm{eV}]$")
        for ax in axes[0]:
            ax.set_ylabel("spin conductivity$\mathrm{~[/\Omega \, m]}$")

        labels = ["T", "L", "Total"]
        for n in np.arange(0,3):
            for j, plot_list in enumerate(plot_list_spin_magnetic_conductivity):
                i = 0
                for key, val in plot_list.items():
                    axes[n][j].set_ylim(window[j])
                    axes[n][j].scatter(e, sigma[n][key], marker=markers[i], s=70, edgecolors=colors[i], facecolor='None', label=val)
                    i = i + 1
                axes[0][j].legend()
            axes[n][0].text(0.05, 0.9, labels[n], bbox=bbox, fontsize="large", transform=axes[n][0].transAxes)

        axes[0][0].legend()
        axes[0][1].legend()

#        plt.show()
        fig.tight_layout()
        plt.savefig(pdf_dir+"cutoff_dependences_gamma"+label_eps+".pdf", bbox_inches = 'tight', dpi=300)
        plt.close()
# }}}


if __name__ == "__main__":
    params = {}
    params['output'] = True
    params['eps'] = 5e-4
    params['slice'] = 100
    params['cutoffs'] = [0.08, 0.09, 0.1]
    params['cutoffs'] = [0.08]
#    cutoff_dependences(params)

    params['cutoff'] = 0.08
    spin_angular_conductivity_T(params, 1)
    spin_angular_conductivity_T(params, 2)
