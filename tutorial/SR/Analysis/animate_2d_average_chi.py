# ____________________________________________________________________________
#
# Show the particle and photon local average energy
#
# _____________________________________________________________________________

import happi
from matplotlib.pyplot import *
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.animation as animation

# ______________________________________________________________________________
# Parameters

# Path to the simulation directory
path = "./"


# open Smilei results
S = happi.Open(path, verbose=False)

# Parameters from the simulation
dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]
simulation_time = S.namelist.Main.simulation_time

PartDiag = S.ParticleBinning(diagNumber=0)

timesteps = PartDiag.getTimesteps()
nt = len(timesteps)

print("Number of iterations: {}".format(nt))

# ______________________________________________________________________________
# Animation function

def animate(i):

    ax.cla()


    # ______________________________________________________________________________
    # Electron Diagnostics

    PartDiag = S.ParticleBinning(diagNumber=0,timesteps=i*500)

    f_density = PartDiag.getData()[0].T
    x_density = PartDiag.getAxis("x")
    y_density = PartDiag.getAxis("y")

    PartDiag = S.ParticleBinning(diagNumber=2,timesteps=i*500)

    f_chi = PartDiag.getData()[0].T
    x_chi = PartDiag.getAxis("x")
    y_chi = PartDiag.getAxis("y")

    f_chi[f_density>0] = f_chi[f_density>0] / f_density[f_density>0]
    print("Iteration {} - Maximal electron quantum parameter: {}".format(timesteps[i],f_chi.max()))
    f_chi[f_density<=0] = np.nan

    im = ax.imshow(
     f_chi, cmap="jet", extent=[x_chi[0],x_chi[-1],y_chi[0],y_chi[-1]], aspect="auto"
    )

    im.set_norm(LogNorm())
    im.set_clim([1e-3,1.])
    ax.set_xlabel(r'$\omega_r x / c$')
    ax.set_ylabel(r'$\omega_r y / c$')
    ax.set_ylim([y_chi[0],y_chi[-1]])
    if i==0:
        cb0 = colorbar(im,format='%g')
        cb0.set_label(r'$\chi_-$')
    t = ax.set_title(r'Iteration {}'.format(i*500))
    t.set_y(1.02)

    return im

# ______________________________________________________________________________
# Figure

fig = figure()
gs = GridSpec(13,5)
ax = subplot(gs[:,:])

animate(0)

fig.tight_layout()

ani = animation.FuncAnimation(fig, animate, np.arange(1, nt),
                              interval=500, blit=False)

show()
