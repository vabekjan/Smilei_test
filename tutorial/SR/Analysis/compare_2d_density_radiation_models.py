# ____________________________________________________________________________
#
# Show the particle and photon local density
#
# _____________________________________________________________________________

# _____________________________________________________________________________
# Importations

import happi
from matplotlib.pyplot import *
import numpy as np
from matplotlib.colors import LogNorm

# ______________________________________________________________________________
# RCparams

rcParams['font.size'] = 20
rcParams['figure.facecolor'] = 'w'
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['axes.labelsize'] = 20

rcParams['xtick.major.size'] = 10
rcParams['ytick.major.size'] = 10

rcParams['xtick.minor.size'] = 5
rcParams['ytick.minor.size'] = 5

rcParams['axes.linewidth'] = 1.5

rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5

rcParams['interactive'] = True

# ______________________________________________________________________________
# Parameters

# Time step for the diagnotics
timestep = 6500

# Read from command line
if len(sys.argv) > 1:
    timestep = int(sys.argv[1])

# ______________________________________________________________________________
# Figure

fig0 = figure(figsize=(16, 9))
gs = GridSpec(13,5)
ax0 = subplot(gs[0:4,:])
ax1 = subplot(gs[4:8,:])
ax2 = subplot(gs[8:12,:])

# ______________________________________________________________________________
# Electron density for the corrected Landau-Lifshitz model

# open Smilei results
S = happi.Open("../Radiation_corrected_Landau_Lifshitz/", verbose=False)

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep).get()

f = np.array(PartDiag["data"][0]).T
x = np.array(PartDiag["x"])
y = np.array(PartDiag["y"])

im0 = ax0.pcolormesh(x,y,f,
                    cmap="jet",
                    shading='nearest')

# ______________________________________________________________________________
# Electron density for the corrected Landau-Lifshitz model

# open Smilei results
S = happi.Open("../Radiation_Niel/", verbose=False)

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep).get()

f = np.array(PartDiag["data"][0]).T
x = np.array(PartDiag["x"])
y = np.array(PartDiag["y"])

im1 = ax1.pcolormesh(x,y,f,
                    cmap="jet",
                    shading='nearest')

# ______________________________________________________________________________
# Electron density for the corrected Landau-Lifshitz model

# open Smilei results
S = happi.Open("../Radiation_Monte_Carlo/", verbose=False)

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep).get()

f = np.array(PartDiag["data"][0]).T
x = np.array(PartDiag["x"])
y = np.array(PartDiag["y"])

im2 = ax2.pcolormesh(x,y,f,
                    cmap="jet",
                    shading='nearest')

# ______________________________________________________________________________
# Figure properties

im0.set_norm(LogNorm())
cb0 = colorbar(im0,format='%g',ax=ax0)
#ax0.set_xlabel(r'$\omega_r x / c$')
ax0.set_ylabel(r'$\omega_r y / c$')
ax0.set_ylim([y[0],y[-1]])
ax0.set_xlim([x[0],x[-1]])
cb0.set_label(r'$n_- / n_c$')
t = ax0.set_title(r'Corrected Landau-Lifshitz at iteration {}'.format(timestep))
t.set_y(1.02)

im1.set_norm(LogNorm())
cb1 = colorbar(im1,format='%g',ax=ax1)
#ax1.set_xlabel(r'$\omega_r x / c$')
ax1.set_ylabel(r'$\omega_r y / c$')
ax1.set_ylim([y[0],y[-1]])
ax1.set_xlim([x[0],x[-1]])
cb1.set_label(r'$n_- / n_c$')
t = ax1.set_title(r'Niel et al. at iteration {}'.format(timestep))
t.set_y(1.02)

im2.set_norm(LogNorm())
cb2 = colorbar(im2,format='%g',ax=ax2)
ax2.set_xlabel(r'$\omega_r x / c$')
ax2.set_ylabel(r'$\omega_r y / c$')
ax2.set_ylim([y[0],y[-1]])
ax2.set_xlim([x[0],x[-1]])
cb2.set_label(r'$n_{-} / n_c$')
t = ax2.set_title(r'Monte-Carlo at iteration {}'.format(timestep))
t.set_y(1.02)

fig0.tight_layout()

show()
