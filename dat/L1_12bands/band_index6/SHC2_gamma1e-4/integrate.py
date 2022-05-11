import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from scipy.optimize import curve_fit
import glob
from matplotlib.ticker import ScalarFormatter

home = "./"
data = "../dat/"
png = "../fig/png/"
svg = "../fig/svg/"

markers = ["o", ",", "D", "v", "^", "<", ">", "s", "p", "1", "2"]
colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']

SHC_files1 = sorted(glob.glob("./spin_conductivity2_mu-*"),reverse=True)
SHC_files2 = sorted(glob.glob("./spin_conductivity2_mu0*"))
SHC_files = SHC_files1 + SHC_files2
mu_full  = []
SHC_full = []

def func(x, b):
    return x / (x**2 + 1e-4**2)**2 * b

sigma_daikei = []
sigma_simpson = []
for file in SHC_files:
    mu = float(file.replace("./spin_conductivity2_mu", "").replace(".csv", ""))
    mu_full.append(mu)
    SHC = pd.read_csv(file,header=0).values
    N = len(SHC[:,0])
    print(mu, N)
    SHC_full.append(SHC)

    param, cov = curve_fit(func, SHC[:,0], SHC[:,6], bounds=((-1e1), (1e1)))

    fig, ax = plt.subplots(1,2,figsize=(15,7))
    ax = ax.flatten()
    ax[0].xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax[0].ticklabel_format(style="sci",  axis="x", scilimits=(0,0))
    ax[0].set_xlim(-1e-3,1e-3)
    n = 100
    x = np.arange(-n,n) / float(n)  * 1e-3
#    ax[0].set_ylim(-20,20)
    ax[0].scatter(SHC[:,0], SHC[:,6], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None')
    f = func(x, *param)
    ax[0].plot(x, f, label='fit: b=%7.5e' % tuple(param))
#    f = func(SHC[:,0], *param)
#    ax[0].plot(SHC[:,0], f, label='fit: b=%7.5e' % tuple(param))

    ax[1].scatter(SHC[:,0], SHC[:,6], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None')
#    f = func(SHC[:,0], *param)
#    ax[1].plot(SHC[:,0], f, label='fit: a=%7.5f, b=%7.5f' % tuple(param))
    f = func(x, *param)
    ax[1].plot(x, f, label='fit: b=%7.5e' % tuple(param))
    plt.legend()
    plt.show()
    plt.close()


    # daikei {{{
    i = 0
    de = (SHC[i+1,0] - SHC[i,0])*5e-1
    conductivity = SHC[i,1:] * de

    N = len(SHC[:,0])
    for i in np.arange(1,N-1):
        de = (SHC[i+1,0] - SHC[i-1,0])*5e-1
        conductivity = conductivity + SHC[i,1:] * de

    i = N-1
    de = (SHC[i,0] - SHC[i-1,0])*5e-1
    conductivity = conductivity + SHC[i,1:] * de

    conductivity = np.append([mu], conductivity)
    sigma_daikei.append(conductivity)
    # }}}

    # simpson {{{
    conductivity = np.empty(27)
    for i in np.arange(1,N-1,3):
        x0 = SHC[i-1,0]
        x1 = SHC[i  ,0]
        x2 = SHC[i+1,0]
        x20 = x2 - x0
        x10 = x1 - x0
        x21 = x2 - x1
        s1 = x20 * (- x2 - 2e0*x0 + 3e0*x1)/x10 * SHC[i-1,1:] / 6e0
        s2 = x20**3 / (x10 * x21) * SHC[i,1:] / 6e0
        s3 = x20 * (2e0*x2 + x0 - 3e0*x1) / x21 * SHC[i+1,1:] / 6e0
        conductivity = conductivity + s1 + s2 + s3

    conductivity = np.append([mu], conductivity)
    sigma_simpson.append(conductivity)
    # }}}

sigma_daikei = np.array(sigma_daikei)
sigma_simpson = np.array(sigma_simpson)

fig, ax = plt.subplots(1,1,figsize=(7,7))
#ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
#ax.ticklabel_format(style="sci",  axis="y", scilimits=(0,0))
ax.scatter(sigma_daikei[:,0], sigma_daikei[:,6], s=3, marker=markers[0], edgecolors=colors[3], facecolor='None', label="daikei")
ax.plot(sigma_daikei[:,0], sigma_daikei[:,6], color=colors[3])
ax.scatter(sigma_simpson[:,0], sigma_simpson[:,6], s=3, marker=markers[1], edgecolors=colors[2], facecolor='None', label="simpson")
ax.plot(sigma_simpson[:,0], sigma_simpson[:,6], color=colors[2])
ax.legend()
plt.show()
plt.close()
