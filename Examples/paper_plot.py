from RVTraceEstimator import RVTraceEstimator

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from PyAstronomy import pyasl
import tayph.system_parameters as sp
import astropy.units as u
import radvel
from PyAstronomy import modelSuite as ms
import seaborn as sns

import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

plt.rc('font', size=12)          # controls default text sizes
#plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('legend', fontsize=14)    # legend fontsize
#plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["ytick.major.width"] = 2
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.major.size"] = 5
plt.rcParams["ytick.major.size"] = 5
plt.rcParams["xtick.minor.size"] = 3.5
plt.rcParams["ytick.minor.size"] = 3.5
#print(plt.rcParams.keys())
plt.rcParams["ytick.minor.width"] = 1
plt.rcParams["xtick.minor.width"] = 1
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['legend.facecolor'] = 'white'
plt.rcParams['legend.edgecolor'] = 'white'
plt.rcParams['legend.framealpha'] = 1


# read RV both times for star and planet without saving the plot.
rv1 = RVTraceEstimator('.', 'paranal')
rv1.plot_doppler_trace('star', save=False)

rv2 = RVTraceEstimator('.', 'paranal')
rv2.plot_doppler_trace('planet', save=False)


# create figure with all traces
fig, ax = plt.subplots(2, 1, sharey=True, sharex=True, figsize=(6,8), constrained_layout=True)

colors=sns.color_palette('colorblind',4 )

ax[0].plot(rv1.RVs_planet, rv1.true_phases, c=colors[0], lw=4, label='planet')
ax[0].plot(rv1.RVs_RM[rv1.transit!=1], rv1.true_phases[rv1.transit!=1], c=colors[2], lw=3, ls='dashed', label='RM')
ax[0].plot(rv1.RVs_star, rv1.true_phases, c=colors[1], lw=3, label='star')
ax[0].plot(rv1.RVs_tel, rv1.true_phases, c=colors[3], lw=4, ls='dotted', label='tellurics')
ax[0].axhline(rv1.true_phases[rv1.transit != 1][0], c='k', ls='dashed')

ax[0].minorticks_on()
ax[0].xaxis.set_tick_params(which='both', top=True, labelsize=15)
ax[0].yaxis.set_tick_params(which='both', right=True, labelsize=15)


ax[1].plot(rv2.RVs_planet, rv2.true_phases, c=colors[0], lw=4)
ax[1].plot(rv2.RVs_RM[rv2.transit!=1], rv2.true_phases[rv2.transit!=1], c=colors[2], lw=3, ls='dashed')
ax[1].plot(rv2.RVs_star, rv2.true_phases, c=colors[1], lw=3)
ax[1].plot(rv2.RVs_tel, rv2.true_phases, c=colors[3], lw=4, ls='dotted')

ax[1].axhline(rv2.true_phases[rv2.transit != 1][0], c='k', ls='dashed')

ax[1].minorticks_on()
ax[1].xaxis.set_tick_params(which='both', top=True, labelsize=15)
ax[1].yaxis.set_tick_params(which='both', right=True, labelsize=15)

ax[1].xaxis.set_major_locator(MultipleLocator(50))
ax[1].set_xlim(np.min(rv1.RVs_planet)-10, np.max(rv2.RVs_star)+10)

fig.supxlabel('Radial velocity [km/s]', fontsize=18)
fig.supylabel('Orbital phase', fontsize=18)

ax[0].legend(loc='lower right')

txt = ax[0].text(x=0.96, y=0.95, s='Stellar rest frame', va='top', ha='right', fontweight=600, transform=ax[0].transAxes, fontsize=13, color='k')
txt = ax[1].text(x=0.96, y=0.95, s='Planetary rest frame', va='top', ha='right', fontweight=600, transform=ax[1].transAxes, fontsize=13, color='k')


plt.savefig('TraceEstimate.pdf')


