# ----------------------------------------------------------------------------------------
# 					SIMULATION PARAMETERS FOR THE PIC-CODE SMILEI
# ----------------------------------------------------------------------------------------

from math import pi, sqrt
import numpy as np
from numpy import s_

lambda0 = 1.0e-6           # laser wavelength in SI
c0 = 3e8                 # speed of light in SI
l0 = 2.*np.pi             # laser wavelength [in code units]
t0 = l0                  # optical cycle
Lsim = [10.*l0, 10.*l0, 10.*l0]  # length of the simulation
Tsim = 18.*l0             # duration of the simulation
resx = 64.               # nb of cells in one laser wavelength
resy = resx
resz = resx
rest = resx*np.sqrt(3.)/0.95 # nb of timesteps in one optical cycle 

xf = Lsim[0]/2		  # focus distance from xmin
w0 = 0.5*l0              # waist
a0 = 4.7 * 1 * 935 * 1.22 / 1.41 * 2.48 * 1.25
FWHMtinI = 0.5*t0
delay = 18*FWHMtinI
eps = l0/(np.pi * w0)
phi01 = 0.0
phi02 = 0.0

# Linearni polarizace
'''
def By_profile(y,z,t):
    return 0

def Bz_profile(y,z,t):
    return a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * np.cos(t + phi0) * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)
'''

# Radialni polarizace
def By_profile(y,z,t):
    factor = (z-Lsim[2]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (15)
    return a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi01)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi01)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)

def Bz_profile(y,z,t):
    factor = (y-Lsim[1]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (15)
    return - a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi01)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi01)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)

def By_profile_counter(y,z,t):
    factor = (z-Lsim[2]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (13)
    return a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi02)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi02)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)

def Bz_profile_counter(y,z,t):
    factor = (y-Lsim[1]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (15)
    return - a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi02)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi02)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)


Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resy,l0/resz],
    grid_length  = Lsim,
    
    number_of_patches = [ 8, 8, 8 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
     
    EM_boundary_conditions = [
        ['silver-muller'],
        ['silver-muller'],
        ['silver-muller'],
    ],
    
    random_seed = smilei_mpi_rank,

    reference_angular_frequency_SI = 2.*np.pi*c0/lambda0
)

LaserOffset(
    box_side               = "xmin",
    space_time_profile     = [ By_profile, Bz_profile ],
    offset                 = xf,
#    extra_envelope          = ,
#    keep_n_strongest_modes = 200,
    angle = 0,
    number_of_processes    = 64,
#    file                   = './LaserOffset0.h5',
    fft_time_window        = 1*Tsim
)

LaserOffset(
    box_side               = "xmax",
    space_time_profile     = [ By_profile_counter, Bz_profile_counter ],
    offset                 = xf,
#    extra_envelope          = ,
#    keep_n_strongest_modes = 200,
    angle = 0,
    number_of_processes    = 64,
#    file                   = './LaserOffset1.h5',
    fft_time_window        = 1.0*Tsim
)


DiagScalar(
    every = [1.*rest,50.*rest,1],
)

DiagFields(
    every = [4.7*rest,9.55*rest,1*rest], 
    fields = ['Ex','Ey','Ez','Bx','By','Bz'],
#    subgrid = s_[96:224,96:224,96:224]
)

