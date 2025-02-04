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
ax0 = subplot(gs[0:4,:])
ax1 = subplot(gs[4:8,:])
ax2 = subplot(gs[8:12,:])

# ______________________________________________________________________________
# Smilei general information

# open Smilei results
S = happi.Open(path, verbose=False)

# Parameters from the simulation
dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]
simulation_time = S.namelist.Main.simulation_time

def my_plot(i1,i2,ax,name):
    PartDiag = S.ParticleBinning(diagNumber=i1,timesteps=timestep)
    
    f_density = PartDiag.getData()[0].T
    x_density = PartDiag.getAxis("x")
    y_density = PartDiag.getAxis("y")
    
    PartDiag = S.ParticleBinning(diagNumber=i2,timesteps=timestep)
    
    f_chi = PartDiag.getData()[0].T
    x_chi = PartDiag.getAxis("x")
    y_chi = PartDiag.getAxis("y")
    
    f_chi[f_density>0] = f_chi[f_density>0] / f_density[f_density>0]
    f_chi[f_density<=0] = np.nan
    
    print("Maximal {} quantum parameter: {}".format(name,f_chi.max()))
    
    return ax.imshow(
        f_chi, cmap="jet", extent=[x_chi[0],x_chi[-1],y_chi[0],y_chi[-1]], aspect="auto", interpolation="nearest"
    )


# ______________________________________________________________________________
# Electron Diagnostics

im0 = my_plot(0,6,ax0,"electron")

# ______________________________________________________________________________
# Positron Diagnostics

im1 = my_plot(1,7,ax1,"positron")

# ______________________________________________________________________________
# Photon Diagnostics

im2 = my_plot(2,8,ax2,"photon")

# ______________________________________________________________________________
# Figure properties

im0.set_norm(LogNorm())
#im0.set_clim([1.,1e-3])
cb0 = colorbar(im0,format='%g',ax=ax0)
#ax0.set_xlabel(r'$\omega_r x / c$')
ax0.set_ylabel(r'$\omega_r y / c$')
cb0.set_label(r'$\chi_-$')
t = ax0.set_title(r'Iteration {}'.format(timestep))
t.set_y(1.02)

im1.set_norm(LogNorm())
#im1.set_clim([1.,1e-3])
cb1 = colorbar(im1,format='%g',ax=ax1)
#ax1.set_xlabel(r'$\omega_r x / c$')
ax1.set_ylabel(r'$\omega_r y / c$')
cb1.set_label(r'$\chi_+$')

im2.set_norm(LogNorm())
#im2.set_clim([1.,1e-3])
cb2 = colorbar(im2,format='%g',ax=ax2)
ax2.set_xlabel(r'$\omega_r x / c$')
ax2.set_ylabel(r'$\omega_r y / c$')
cb2.set_label(r'$\chi_{\gamma}$')

fig0.tight_layout()

show()
