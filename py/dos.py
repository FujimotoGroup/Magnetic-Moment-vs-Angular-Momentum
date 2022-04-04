import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import griddata
import pandas as pd
import configparser as cnf
from scipy.optimize import curve_fit
import scipy as spy

home = "../"
data = "../dat/"
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

numeric = config['numeric']

colors =['k', 'b', 'g', 'r', 'tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'navy', 'indigo']


delta = (EL2-EL0)/2.
dos_file = data+'DOS_isotropic_plus.csv'
df = pd.read_csv(dos_file,header=None).values
e   = df[:,0]
dos = df[:,1]
n = 1000.0
def func(x, p):
    return abs(x)*np.lib.scimath.sqrt(x*x - delta*delta)*p/(2.0*np.pi*np.pi)

ei = np.linspace(e.min()*n, e.max()*n, 1000) / n
parameter_initial = np.array([1.0])
popt, pcov = curve_fit(func, e, dos, p0=parameter_initial)
print("c =", popt[0], "c =", pcov[0])
fit_c = popt[0]
#fit = abs(ei)*.sqrt(ei*ei - delta*delta)*0.2766/(2.0*np.pi*np.pi)


fig = plt.figure(figsize=(12,5))
ax1 = fig.add_subplot(1,1,1)
ax1.plot(ei, fit)
ax1.scatter(e, dos)
#ax1.plot(e_range, fit, 'r-', label= 'fit: r=%1.1f, - %2.2f * x^%2.2f' % tuple([ np.exp(fit_c)]))
plt.show()
