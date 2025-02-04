# ____________________________________________________________________________
#
# Show the energy balance (time evolution of the kinetic energy)
#
# _____________________________________________________________________________

# _____________________________________________________________________________
# Importations

import happi
from matplotlib.pyplot import *
import numpy as np

# ______________________________________________________________________________
# RCparams

rcParams['figure.facecolor'] = 'w'
rcParams['font.size'] = 20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['axes.labelsize'] = 25

rcParams['xtick.major.size'] = 10
rcParams['ytick.major.size'] = 10

rcParams['xtick.minor.size'] = 5
rcParams['ytick.minor.size'] = 5

rcParams['axes.linewidth'] = 1.5

rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2

rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5

rcParams['interactive'] = False

# ______________________________________________________________________________
# Figure

fig0 = figure(figsize=(18, 6))
gs = GridSpec(2, 3)
ax0 = subplot(gs[:,:])

# ______________________________________________________________________________
# Corrected Landau-Lifshitz

# open Smilei results
S = happi.Open("../Radiation_corrected_Landau_Lifshitz/", verbose=False)

dt = S.namelist.Main.timestep

times = np.array(S.Scalar("Ukin_electron").get()["times"])

# Kinetic electrone energy
ukin_electron = np.array(S.Scalar("Ukin_electron").get()["data"])

# Radiated energy wihtout the macro-photons
urad = np.array(S.Scalar("Urad").get()["data"])

# Total kinetic energy
utot = urad + ukin_electron

ax0.plot(times*dt,ukin_electron/utot,color='b',lw=2,ls='-',label="Electron kinetic energy")
ax0.plot(times*dt,urad/utot,color='purple',lw=2,ls='-',label="Radiation")
ax0.plot(times*dt,utot/utot[0],color='red',lw=2,ls='-',label="Total kinetic energy")

# ______________________________________________________________________________
# Landau-Lifshitz

# open Smilei results
S = happi.Open("../Radiation_Niel/", verbose=False)

dt = S.namelist.Main.timestep

times = np.array(S.Scalar("Ukin_electron").get()["times"])

ukin_electron = np.array(S.Scalar("Ukin_electron").get()["data"])

# Radiated energy wihtout the macro-photons
urad = np.array(S.Scalar("Urad").get()["data"])

# Total kinetic energy
utot = urad + ukin_electron

ax0.plot(times*dt,ukin_electron/utot,color='b',lw=2,ls='--')
ax0.plot(times*dt,urad/utot,color='purple',lw=2,ls='--')
ax0.plot(times*dt,utot/utot[0],color='red',lw=2,ls='--')

# ______________________________________________________________________________
# Monte-Carlo

# open Smilei results
S = happi.Open("../Radiation_Monte_Carlo/", verbose=False)

dt = S.namelist.Main.timestep

times = np.array(S.Scalar("Ukin_electron").get()["times"])

ukin_electron = np.array(S.Scalar("Ukin_electron").get()["data"])

# Radiated energy wihtout the macro-photons
urad = np.array(S.Scalar("Urad").get()["data"])

# Total kinetic energy
utot = urad + ukin_electron

ax0.plot(times*dt,ukin_electron/utot,color='b',lw=2,ls=':')
ax0.plot(times*dt,urad/utot,color='purple',lw=2,ls=':')
ax0.plot(times*dt,utot/utot[0],color='red',lw=2,ls=':')

# ______________________________________________________________________________
# Figure properties

ax0.plot([0,0],[0,0],ls='-',color='k',label='Corrected Landau-Lifshitz')
ax0.plot([0,0],[0,0],ls='--',color='k',label='Niel')
ax0.plot([0,0],[0,0],ls=':',color='k',label='Monte-Carlo')
t = ax0.set_title('Energy balance')
ax0.set_xlabel(r'$\omega_r t$')
ax0.set_ylabel(r'$\varepsilon / \varepsilon_{tot}$')
ax0.set_xlim([200,300])
ax0.set_ylim([0,1.1])
#ax0.set_yscale('log')

ax0.legend(loc='best', fontsize=15)

fig0.tight_layout()

show()
