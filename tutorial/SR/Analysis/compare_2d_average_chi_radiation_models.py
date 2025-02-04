# ____________________________________________________________________________
#
# Show the particle and photon local average energy
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

rcParams['font.size'] = 15
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
timestep = 6000

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
# Electron quantum parameter for the corrected Landau-Lifshitz model

# open Smilei results
S = happi.Open("../Radiation_corrected_Landau_Lifshitz/", verbose=False)

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep).get()

f_density = np.array(PartDiag["data"][0]).T
x_density = np.array(PartDiag["x"])
y_density = np.array(PartDiag["y"])

PartDiag = S.ParticleBinning(diagNumber=2,timesteps=timestep).get()

f_chi = np.array(PartDiag["data"][0]).T
x_chi = np.array(PartDiag["x"])
y_chi = np.array(PartDiag["y"])

f_chi[f_density>0] = f_chi[f_density>0] / f_density[f_density>0]

print("Maximal quantum parameter: {}".format(f_chi.max()))

im0 = ax0.pcolormesh(x_chi,y_chi,f_chi,
                    cmap="jet",
                    shading='nearest')

# ______________________________________________________________________________
# Electron quantum parameter for the Niel model

# open Smilei results
S = happi.Open("../Radiation_Niel/", verbose=False)

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep).get()

f_density = np.array(PartDiag["data"][0]).T
x_density = np.array(PartDiag["x"])
y_density = np.array(PartDiag["y"])

PartDiag = S.ParticleBinning(diagNumber=2,timesteps=timestep).get()

f_chi = np.array(PartDiag["data"][0]).T
x_chi = np.array(PartDiag["x"])
y_chi = np.array(PartDiag["y"])

f_chi[f_density>0] = f_chi[f_density>0] / f_density[f_density>0]

print("Maximal quantum parameter: {}".format(f_chi.max()))

im1 = ax1.pcolormesh(x_chi,y_chi,f_chi,
                    cmap="jet",
                    shading='nearest')

# ______________________________________________________________________________
# Electron quantum parameter for the Monte-Carlo model

# open Smilei results
S = happi.Open("../Radiation_Monte_Carlo/", verbose=False)

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep).get()

f_density = np.array(PartDiag["data"][0]).T
x_density = np.array(PartDiag["x"])
y_density = np.array(PartDiag["y"])

PartDiag = S.ParticleBinning(diagNumber=2,timesteps=timestep).get()

f_chi = np.array(PartDiag["data"][0]).T
x_chi = np.array(PartDiag["x"])
y_chi = np.array(PartDiag["y"])

f_chi[f_density>0] = f_chi[f_density>0] / f_density[f_density>0]

print("Maximal quantum parameter: {}".format(f_chi.max()))

im2 = ax2.pcolormesh(x_chi,y_chi,f_chi,
                    cmap="jet",
                    shading='nearest')

# ______________________________________________________________________________
# Figure properties

im0.set_norm(LogNorm())
im0.set_clim([1.,1e-3])
cb0 = colorbar(im0,format='%g',ax=ax0)
ax0.set_ylabel(r'$\omega_r y / c$')
ax0.set_ylim([y_chi[0],y_chi[-1]])
ax0.set_xlim([x_chi[0],x_chi[-1]])
cb0.set_label(r'$\chi_-$')
t = ax0.set_title(r'Corrected Landau-Lifshitz at iteration {}'.format(timestep))
t.set_y(1.02)

im1.set_norm(LogNorm())
im1.set_clim([1.,1e-3])
cb1 = colorbar(im1,format='%g',ax=ax1)
ax1.set_ylabel(r'$\omega_r y / c$')
ax1.set_ylim([y_chi[0],y_chi[-1]])
ax1.set_xlim([x_chi[0],x_chi[-1]])
cb1.set_label(r'$\chi_-$')
t = ax1.set_title(r'Niel et al. at iteration {}'.format(timestep))
t.set_y(1.02)

im2.set_norm(LogNorm())
im2.set_clim([1.,1e-3])
cb2 = colorbar(im2,format='%g',ax=ax2)
ax2.set_xlabel(r'$\omega_r x / c$')
ax2.set_ylabel(r'$\omega_r y / c$')
ax2.set_ylim([y_chi[0],y_chi[-1]])
ax2.set_xlim([x_chi[0],x_chi[-1]])
cb2.set_label(r'$\chi_-$')
t = ax2.set_title(r'Monte-Carlo at iteration {}'.format(timestep))
t.set_y(1.02)

fig0.tight_layout()

show()
