# ____________________________________________________________________________
#
# Show the particle energy spectrum at a given timestep
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

# Iteration number
timestep = 8000

# Read from command line
if len(sys.argv) > 1:
    timestep = int(sys.argv[1])

# ______________________________________________________________________________
# Smilei general information

# open Smilei results
S = happi.Open(path, verbose=False)

# Parameters
dt = S.namelist.Main.timestep
dx = S.namelist.Main.cell_length[0]
dy = S.namelist.Main.cell_length[1]
simulation_time = S.namelist.Main.simulation_time



# ______________________________________________________________________________
# Figure

fig0 = figure()
gs = GridSpec(2, 3)
ax0 = subplot(gs[:,:])

# ______________________________________________________________________________
# Particle diagnostics

PartDiag = S.ParticleBinning(diagNumber=9,timesteps = timestep)
gamma = np.array(PartDiag.getAxis("gamma"))*0.511
n = np.array(PartDiag.getData()[0])

print(np.shape(n),np.shape(gamma))

ax0.plot(gamma,n,color='b',label='Electron',lw=2)

PartDiag = S.ParticleBinning(diagNumber=10,timesteps = timestep)
gamma = np.array(PartDiag.getAxis("gamma"))*0.511
n = np.array(PartDiag.getData()[0])
ax0.plot(gamma,n,color='r',label='Positron',lw=2)

PartDiag = S.ParticleBinning(diagNumber=11,timesteps = timestep)
gamma = np.array(PartDiag.getAxis("gamma"))*0.511
n = np.array(PartDiag.getData()[0])
ax0.plot(gamma,n,color='g',label='Photon',lw=2)

# ______________________________________________________________________________
# Figure properties

t = ax0.set_title('Energy spectrum, iteration {}'.format(timestep))
ax0.set_xlabel(r'$\varepsilon\ (MeV)$')
ax0.set_ylabel(r'$Weight$')
ax0.set_xscale('log')
ax0.set_yscale('log')

ax0.legend(loc='best', fontsize=15)

fig0.tight_layout()

show()
