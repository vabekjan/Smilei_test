# ____________________________________________________________________________
#
# Show the particle and photon local average energy
#
# _____________________________________________________________________________

import happi
from matplotlib.pyplot import *
import numpy as np
from matplotlib.colors import LogNorm

# ______________________________________________________________________________
# Parameters

# Path to the simulation directory
path = "./"

# Time step for the diagnotics
timestep = 5500

# Read from command line
if len(sys.argv) > 1:
    timestep = int(sys.argv[1])

# ______________________________________________________________________________
# Figure

fig0 = figure()
gs = GridSpec(13,5)
ax0 = subplot(gs[:,:])

# ______________________________________________________________________________
# Smilei general information

# open Smilei results
S = happi.Open(path, verbose=False)

# Parameters from the simulation
dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]
simulation_time = S.namelist.Main.simulation_time

# ______________________________________________________________________________
# Electron Diagnostics

PartDiag = S.ParticleBinning(diagNumber=0,timesteps=timestep)

f_density = PartDiag.getData()[0].T
x_density = PartDiag.getAxis("x")
y_density = PartDiag.getAxis("y")

PartDiag = S.ParticleBinning(diagNumber=2,timesteps=timestep)

f_chi = PartDiag.getData()[0].T
x_chi = PartDiag.getAxis("x")
y_chi = PartDiag.getAxis("y")

f_chi[f_density>0] = f_chi[f_density>0] / f_density[f_density>0]
print("Maximal electron quantum parameter: {}".format(f_chi.max()))
f_chi[f_density<=0] = np.nan


im0 = ax0.imshow(
f_chi, cmap="jet", extent=[x_chi[0],x_chi[-1],y_chi[0],y_chi[-1]], aspect="auto"
)

# ______________________________________________________________________________
# Figure properties

im0.set_norm(LogNorm())
im0.set_clim([1e-3,1.])
cb0 = colorbar(im0,format='%g',ax=ax0)
ax0.set_xlabel(r'$\omega_r x / c$')
ax0.set_ylabel(r'$\omega_r y / c$')
ax0.set_ylim([y_chi[0],y_chi[-1]])
cb0.set_label(r'$\chi_-$')
t = ax0.set_title(r'Iteration {}'.format(timestep))
t.set_y(1.02)

fig0.tight_layout()

show()
