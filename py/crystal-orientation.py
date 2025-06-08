import numpy as np
import pandas as pd
import configparser as cnf

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from decimal import Decimal, Context

from scipy import interpolate

angstrom = 1e-10 # [m]
hbar     = 6.582119569e-16 # [eV s]
mass     = 9.10938356e-31 # [kg]
charge   = 1.60217662e-19 # [C]
muB      = 9.2740100783e-24 # [J/T]
muBeV    = 5.7883818060e-5 # [eV/T]

markers = ["o", "v", "D", ",", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

home = "../"
data0 = "../dat/"
figure_size_unit = 5

config = cnf.ConfigParser()
config.read(home+'config.ini', encoding='utf-8')

physics = config['physics']
a0 = float(physics.get('a'))  #  4.5332e0 # angstrom
c0 = float(physics.get('c'))  # 11.7967e0 # angstrom
g0 = float(physics.get('g0')) # 1.3861e0 # angstrom^-1
bands = int(physics.get('bands'))
bandsT = int(physics.get('bandsT'))
bandsL = int(physics.get('bandsL'))
lowest_T = int(physics.get('lowest_band_T'))
lowest_L = int(physics.get('lowest_band_L'))
lower_band_L = 0
upper_band_L = 2
#bandsL = 12
#lower_band_L = 4
#upper_band_L = 6
valleys = {"T":0, "L1":1, "L2":2, "L3":3}

cutoffs = ["0.08", "0.10"]
cutoff = cutoffs[0]
param = "cutoff"+cutoff
epsilon_T = 0.7e-4
epsilon_L = 1.6e-4
label_damping_T = "_gamma0.7e-4"
label_damping_L = "_gamma1.6e-4"
parameter = "T"+str(bandsT)+"bands_L"+str(bandsL)+"bands"+"_"+param+"eV"

key_L = [2 , 3, 1]

a = np.array([[ -a0/2e0, -a0/(2e0*np.sqrt(3e0)), c0/3e0 ],
              [  a0/2e0, -a0/(2e0*np.sqrt(3e0)), c0/3e0 ],
              [     0.0,  a0/(    np.sqrt(3e0)), c0/3e0 ]])

b = np.array([[- g0, -    g0/np.sqrt(3e0), (a0/c0)*g0 ],
              [  g0, -    g0/np.sqrt(3e0), (a0/c0)*g0 ],
              [ 0e0,  2e0*g0/np.sqrt(3e0), (a0/c0)*g0 ]])

def check():
    for i in np.arange(0,3):
        for j in np.arange(0,3):
            print("a_{}*b_{} = {}".format(i+1,j+1,np.dot(a[i], b[j])/(2e0*np.pi)))

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
bbox = {
    "facecolor" : "white",
    "edgecolor" : "black",
    "boxstyle" : "round, pad=0.5",
    "linewidth" : 1
}

def plane_normal(hkl):
    """(hkl) → 単位法線ベクトル (直交座標)"""
    h, k, l = hkl
    G = h*b[0] + k*b[1] + l*b[2]
    return G / np.linalg.norm(G)

def get_index(axis):
    from itertools import product

    axes = ['x', 'y', 'z']          # 座標ラベル
    idx = {}
    labels = []                     # ['xxx', 'xxy', ...] を自前で作る場合
    idx2ija = {}                    # {index: (i, j, α)}
    lab2ija = {}                    # {'xxx':(0,0,0), ...}

    for n, (i, j, a) in enumerate(product(axes, repeat=3)):
        label = f'{i}{j}{a}'        # 例: 'xxy'
        labels.append(label)

        # 数値インデックス (0,1,2) に変換して保存
        idx_tuple = (axes.index(i), axes.index(j), axes.index(a))
        idx2ija[n]   = idx_tuple
        lab2ija[label] = idx_tuple
        idx[label] = n

#    # -------- 使い方例 --------
#    print(labels)           # → ['xxx', 'xxy', 'xxz', 'xyx', ... , 'zzz']
#    print(idx2ija[13])      # 13番目の成分は (1,1,1)＝'yyy'
#    print(lab2ija['zxy'])   # 'zxy' は (2,0,1)
#    print(idx['zxy'])       # 19

    filtered = [idx[lab] for lab in labels if lab[2] == axis]

    return filtered


def load_electric_conductivity():
    conductivity_valley = []

    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values
    d = df.values
    conductivity_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
        d4 = pd.read_csv(readfile,header=0).values

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/conductivity_eps'+label+'.csv'
        d6 = pd.read_csv(readfile,header=0).values

        conductivity_valley.append(np.append(d4, d6, axis=0))

    x, conductivity, conductivity_at_Fermi_level = interpolate_conductivity(conductivity_valley)

    return conductivity, conductivity_at_Fermi_level

def load_spin_conductivity(operator, index):
    conductivity_valley = []
    data = data0+'T_'+str(bandsT)+'bands/band_index4/'+param+label_damping_T+'/'
    label = "{:.6f}".format(epsilon_T)
    readfile = data+'mu-dependence/spin-'+operator+'-conductivity'+index+'_eps'+label+'.csv'
#    print(readfile)
    df = pd.read_csv(readfile,header=0)
    titles = df.columns.values[1:]

#    idxs = [get_index("x"), get_index("y"), get_index("z")]
#    for j, idx in enumerate(idxs):
#        print(np.array([titles[i] for i in idx]).reshape(3,3))

    d = df.values
    conductivity_valley.append(d)

    for i in np.arange(3):
        valley = key_L[i]
        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(lower_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-'+operator+'-conductivity'+index+'_eps'+label+'.csv'
#        print(readfile)
        d4 = pd.read_csv(readfile,header=0).values

        data = data0+'L'+str(valley)+'_'+str(bandsL)+'bands/band_index'+str(upper_band_L)+'/'+param+label_damping_L+'/'
        label = "{:.6f}".format(epsilon_L)
        readfile = data+'mu-dependence/spin-'+operator+'-conductivity'+index+'_eps'+label+'.csv'
#        print(readfile)
        d6 = pd.read_csv(readfile,header=0).values

        if index == "2":
            val = np.copy(d4[-1, :])
            val[0] = 0e0
            d6 = d6 + val

        conductivity_valley.append(np.append(d4, d6, axis=0))

    return titles, conductivity_valley

def load_magnetic_conductivity():
    titles, sigma1_valley = load_spin_conductivity("magnetic", "1")
    _,      sigma2_valley = load_spin_conductivity("magnetic", "2")
    conductivity = [sigma1_valley, sigma2_valley]

#    plot_raw_spin_conductivity_valleys("m", titles, sigma1_valley, "spin_magnetic_conductivity1_"+parameter)
#    plot_raw_spin_conductivity_valleys("m", titles, sigma2_valley, "spin_magnetic_conductivity2_"+parameter)

    conductivity, conductivity_at_Fermi_level = interpolate_spin_conductivity(*conductivity)

#    plot_spin_conductivity("m", conductivity, "spin_magnetic_conductivity_"+parameter)

    return conductivity, conductivity_at_Fermi_level

def load_angular_conductivity():
    titles, sigma1_valley = load_spin_conductivity("angular", "1")
    _,      sigma2_valley = load_spin_conductivity("angular", "2")
    conductivity = [sigma1_valley, sigma2_valley]

#    plot_raw_spin_conductivity_valleys("a", titles, sigma1_valley, "spin_angular_conductivity1_"+parameter)
#    plot_raw_spin_conductivity_valleys("a", titles, sigma2_valley, "spin_angular_conductivity2_"+parameter)

    conductivity, conductivity_at_Fermi_level = interpolate_spin_conductivity(*conductivity)

#    plot_spin_conductivity("a", conductivity, "spin_angular_conductivity_"+parameter)

    return conductivity, conductivity_at_Fermi_level

def plot_raw_spin_conductivity_valleys(sign, titles, conductivity, output):
    xticks = np.linspace(-float(cutoff), float(cutoff), 4+1)
    xticklabels = ["$"+decimal_normalize(tick*1e3)+"$" for tick in xticks]

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
#        ax.set_ylim(window)

    plot_lists = [6, 16, 20]
    for s, i in enumerate(plot_lists):
        for key, val in valleys.items():
            axes[s].scatter(conductivity[val][:,0], conductivity[val][:,i], s=4, label=key)
    axes[-1].legend()
    plt.show()

#    plt.savefig(png+output+".png", bbox_inches='tight', dpi=300)
#    plt.savefig(pdf+output+".pdf")
    #plt.show()
    plt.close()

def plot_spin_conductivity(sign, conductivity, output):
    plot_list1  = { 6:"xyz", 16:"yzx", 20:"zxy"}
    plot_list2  = {12:"yxz",  8:"xzy", 22:"zyx"}
    plot_list3  = { 1:"xxx", 13:"yyx",  5:"xyy", 11:"yxy"}
    plot_list4a = { 2:"xxy",  3:"xxz",  4:"xyx",  7:"xzx",  9:"xzz", 10:"yxx", 14:"yyy", 15:"yyz", 17:"yzy"}
    plot_list4b = {18:"yzz", 19:"zxx", 21:"zxz", 23:"zyy", 24:"zyz", 25:"zzx", 26:"zzy", 27:"zzz"}
    plot_lists = [ plot_list1, plot_list2, plot_list3, plot_list4a, plot_list4b ]

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
                axes[n][j].scatter(conductivity[n][0], conductivity[n][key], marker=markers[i], s=20, edgecolors=colors[i], facecolor='None', label="$\sigma_{"+sign+", "+flow+force+"}"+"^{"+spin+"}$")
                i += 1
            axes[0][j].legend()
            k += 1

    plt.show()
    fig.tight_layout()
#    fig.savefig(pdf+output+".pdf", bbox_inches='tight')
    plt.close()

def interpolate_conductivity(conductivity_valley):
    n = 100
    e_range = np.linspace(-float(cutoff)*n, float(cutoff)*n, int(n)-1)
    x = e_range / n

    sigma_T = conductivity_valley[0]

    e = conductivity_valley[1][:,0]
    conductivity = conductivity_valley[1] + conductivity_valley[2] + conductivity_valley[3]
    conductivity[:,0] = e
    sigma_L = conductivity

    n = sigma_T.shape[1]
    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(x) for i in np.arange(1,n)])
    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(x) for i in np.arange(1,n)])

    conductivity = np.array([is_T, is_L])

    y = 0e0
    is_T = np.array([interpolate.interp1d(sigma_T[:,0], sigma_T[:,i])(y) for i in np.arange(1,n)])
    is_L = np.array([interpolate.interp1d(sigma_L[:,0], sigma_L[:,i])(y) for i in np.arange(1,n)])
    conductivity_at_Fermi_level = np.sum(np.array([is_T, is_L]), axis=0)

    return x, conductivity, conductivity_at_Fermi_level

def interpolate_spin_conductivity(conductivity1_valley, conductivity2_valley):
    x, is1, f1 = interpolate_conductivity(conductivity1_valley)
    x, is2, f2 = interpolate_conductivity(conductivity2_valley)

    conductivity = np.sum(np.array([is1, is2]), axis=0)
    conductivity_total = np.sum(conductivity, axis=0)

    pad = np.zeros((conductivity.shape[0], 1, conductivity.shape[2]))   # (2,1,99)
    pad[0] = x[None, None, :]
    pad[1] = x[None, None, :]
    conductivity = np.concatenate([pad, conductivity], axis=1)

    conductivity_total = np.concatenate([np.array([x]), conductivity_total])

    conductivity_at_Fermi_level = np.sum(np.array([f1, f2]), axis=0)

    return conductivity_total, conductivity_at_Fermi_level

def orientation(hkl, sigma_e, sigma_m, sigma_a):
    axis = ['x', 'y', 'z']
    idxs = [get_index("x"), get_index("y"), get_index("z")]

    m = []
    a = []
    print("- - - - - - - - - - - - - - - - - - - -")
    n = plane_normal(hkl)
    print("n({}{}{}) = {}".format(*hkl, n))

    sigma_e_perp = float(n @ sigma_e @ n)
#    sigma_e_perp = 2.4e5
    print("σ_e({}{}{}) = {:.4f} S/m".format(*hkl,sigma_e_perp))

    for j, idx in enumerate(idxs):
        sigma = np.array([sigma_m[i] for i in idx]).reshape(3,3)
        sigma_m_perp = float(n @ sigma @ n)
        print("σ_m^{}({}{}{}) = {:.4f} S/m, theta = {:.4f}".format(axis[j], *hkl, sigma_m_perp, sigma_m_perp/sigma_e_perp))
        m.append(sigma_m_perp)

    for j, idx in enumerate(idxs):
        sigma = np.array([sigma_a[i] for i in idx]).reshape(3,3)
        sigma_a_perp = float(n @ sigma @ n)
        print("σ_a^{}({}{}{}) = {:.4f} S/m, theta = {:.6f}".format(axis[j], *hkl, sigma_a_perp, sigma_a_perp/sigma_e_perp))
        a.append(sigma_a_perp)

#    m = np.array(m)
#    a = np.array(a)
#    print(np.dot(m,n)/sigma_e_perp)
#    print(np.dot(a,n)/sigma_e_perp)

# ------------------------------------------------------------
# 1. (h k l) → 回転行列 R を作るユーティリティ
# ------------------------------------------------------------
def rotation_matrix_from_hkl(h, k, l, prefer_x_axis=True):
    """
    面 (h k l) を z' (=面外) に持つ試料座標系を作る。
      - x' : 面内の任意の単位ベクトル
      - y' : x' × z'
      - z' : 面法線 [h k l] / |[h k l]|
    戻り値 R の列が新基底 (結晶座標で表現)。
    """
    n = plane_normal(hkl)
    if np.allclose(n, 0):
        raise ValueError("(h,k,l)=(0,0,0) は面を定義できません")
    z_p = n / np.linalg.norm(n)

    # ref を手動指定
    if prefer_x_axis:
        ref = np.array([1.0, 0.0, 0.0])
        if np.allclose(np.abs(z_p), [1,0,0], atol=1e-7):   # ref が z′ と平行になる例外
            ref = np.array([0.0, 1.0, 0.0])
    else:
        ref = np.array([0.0, 1.0, 0.0])

    # 面内の x' を作り，正規化
    x_p = ref - np.dot(ref, z_p) * z_p
    x_p /= np.linalg.norm(x_p)

    # 右手系で y' を決定
    y_p = np.cross(z_p, x_p)

    # 列が新基底
    return np.vstack([x_p, y_p, z_p]).T          # shape (3,3)

# ------------------------------------------------------------
# 2. テンソル回転（rank-2, rank-3 両対応）
# ------------------------------------------------------------
def transform_rank2(tensor, R):
    return R @ tensor @ R.T

def transform_rank3(tensor, R):
    return np.einsum('ia,jb,kc,abc->ijk', R, R, R, tensor)

# ------------------------------------------------------------
# 3. 面内完全多結晶平均で法線方向のスピン Hall 伝導度
# ------------------------------------------------------------
def spin_hall_perp_average(n, sigma_rank3_sample, n_phi=720):
    z = 2                                  # z' インデックス
    phi = np.linspace(0, 2e0*np.pi, n_phi, endpoint=False)

    acc = []
    for ang in phi:
        E = np.array([np.cos(ang), np.sin(ang), 0e0])
        zxE = - np.array([-np.sin(ang), np.cos(ang), 0e0])
        js_z = [0e0, 0e0, 0e0]
        for i in np.arange(3):
            js_z[i] = np.dot(sigma_rank3_sample[:, z, i], E)
        jspin = np.dot(js_z, zxE)
        acc.append(jspin)

    fig, ax = plt.subplots(1,1,figsize=(figure_size_unit,figure_size_unit))
    ax.plot(phi,acc)
    plt.show()

    acc = np.sum(np.array(acc), axis=0)

    return acc / n_phi

# ------------------------------------------------------------
# 4. まとめ関数：結晶テンソル → ⟨σSH⟩(hkl)
# ------------------------------------------------------------
def sigma_spin_hall_perp(sign, h, k, l, n, sigma_crystal, *, n_phi=720):
    R = rotation_matrix_from_hkl(h, k, l)
    sigma_sample = transform_rank3(sigma_crystal, R)
    sigma_perp = spin_hall_perp_average(n, sigma_sample, n_phi=n_phi)

    print("σ_{}({}{}{}) = {:.4e} S/m".format(sign, *hkl, sigma_perp))

if __name__ == "__main__":
#    check()

    _, sigma_e = load_electric_conductivity()
    _, sigma_m = load_magnetic_conductivity()
    _, sigma_a = load_angular_conductivity()

    sigma_e = sigma_e.reshape(3,3)

    m = []
    a = []
    idxs = [get_index("x"), get_index("y"), get_index("z")]
    for j, idx in enumerate(idxs):
        sigma = np.array([sigma_m[i] for i in idx]).reshape(3,3)
        m.append(sigma)
        sigma = np.array([sigma_a[i] for i in idx]).reshape(3,3)
        a.append(sigma)

    sigma_m = np.array(m)
    sigma_a = np.array(a)

    n = np.array([1e0,0e0,0e0])

    hkl = (1,1,0)
##    orientation(hkl, sigma_e, sigma_m, sigma_a)
    sigma_spin_hall_perp("m", *hkl, n, sigma_m)
    sigma_spin_hall_perp("a", *hkl, n, sigma_a)

    hkl = (1,1,1)
#    orientation(hkl, sigma_e, sigma_m, sigma_a)
    sigma_spin_hall_perp("m", *hkl, n, sigma_m)
    sigma_spin_hall_perp("a", *hkl, n, sigma_a)

#
#    hkl = (1,1,1)
#    orientation(hkl, sigma_e, sigma_m, sigma_a)

    hkl = (-1,-1,2)
#    orientation(hkl, sigma_e, sigma_m, sigma_a)

