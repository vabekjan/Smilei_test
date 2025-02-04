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
gs = GridSpec(9,5)
ax0 = subplot(gs[0:4,:])
ax1 = subplot(gs[4:8,:])

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
# Ey electric field

FieldDiag = S.Field(0, "Ey", timesteps=timestep)

f = FieldDiag.getData()[0].T
x = FieldDiag.getAxis("x")
y = FieldDiag.getAxis("y")

im0 = ax0.imshow(
 f, cmap="jet", extent=[x[0],x[-1],y[0],y[-1]], aspect="auto"
)

# ______________________________________________________________________________
# Bz electric field

FieldDiag = S.Field(0, "Bz", timesteps=timestep)

f = FieldDiag.getData()[0].T
x = FieldDiag.getAxis("x")
y = FieldDiag.getAxis("y")

im1 = ax1.imshow(
 f, cmap="jet", extent=[x[0],x[-1],y[0],y[-1]], aspect="auto"
)

# ______________________________________________________________________________
# Figure properties

cb0 = colorbar(im0,format='%g',ax=ax0)
ax0.set_xlabel(r'$\omega_r x / c$')
ax0.set_ylabel(r'$\omega_r y / c$')
ax0.set_ylim([y[0],y[-1]])
ax0.set_xlim([x[0],x[-1]])
cb0.set_label(r'$e E_y / m \omega c$')

cb1 = colorbar(im1,format='%g',ax=ax1)
ax1.set_xlabel(r'$\omega_r x / c$')
ax1.set_ylabel(r'$\omega_r y / c$')
ax1.set_ylim([y[0],y[-1]])
ax1.set_xlim([x[0],x[-1]])
cb1.set_label(r'$e B_z / m \omega$')

#fig0.tight_layout()

show()
