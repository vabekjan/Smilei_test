# ____________________________________________________________________________
#
# Show the time evolution of the number of macro-particles for each species 
# during the simulation.
#
# _____________________________________________________________________________

# _____________________________________________________________________________
# Importations

import happi
from matplotlib.pyplot import *
import numpy as np

# ______________________________________________________________________________
# Parameters

# Path to the simulation directory
path = "./"

# ______________________________________________________________________________
# Smilei general information

# open Smilei results
S = happi.Open(path, verbose=False)

# Parameters
dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]
simulation_time = S.namelist.Main.simulation_time

print(' Space steps:%f %f'%(dx,dy))


# ______________________________________________________________________________
# Figure

fig0 = figure(figsize=(10, 6))
gs = GridSpec(2, 2)
ax0 = subplot(gs[:,:])

# ______________________________________________________________________________
# Scalar discontinuous

times = np.array(S.Scalar("Ntot_electron").get()["times"])

N_electron = np.array(S.Scalar("Ntot_electron").get()["data"])

N_positron = np.array(S.Scalar("Ntot_positron").get()["data"])

N_photon = np.array(S.Scalar("Ntot_photon").get()["data"])


ax0.plot(times*dt,N_electron,color='b',label="Electron")
ax0.plot(times*dt,N_positron,color='r',label="Positron")
ax0.plot(times*dt,N_photon,color='g',label="Photon")

# ______________________________________________________________________________
# Figure properties

ax0.set_xlabel(r'$\omega_r t$')
ax0.set_ylabel('Number of particles')
ax0.set_yscale('log')
ax0.set_xlim([200,315])

# Legend properties
ax0.legend(loc='lower left', fontsize=15)
fig0.tight_layout()

show()
