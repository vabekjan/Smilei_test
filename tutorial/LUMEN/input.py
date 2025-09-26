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
Lsim = [5.*l0, 10.*l0, 10.*l0]  # length of the simulation
Tsim = 12.*l0 #20            # duration of the simulation
resx = 64.               # nb of cells in one laser wavelength
resy = resx
resz = resx
rest = resx*np.sqrt(3.)/0.95 # nb of timesteps in one optical cycle 

xf = Lsim[0]/2		  # focus distance from xmin
w0 = 0.5*l0 #1.59              # waist
a0 = 4.7 * 1 * 3000
FWHMtinI = 0.5*t0
delay = 18*FWHMtinI
eps = l0/(np.pi * w0)
phi01 = 0
phi02 = 0 # np.pi # 0 = max b center; np.pi = max e center

n0 =   0.000001 #0.1
electron_gev = 1
gamma = electron_gev*1e3/0.511          # Electron beam gamma factor
v = math.sqrt(1 - 1./gamma**2)          # electron beam initial velocity

#def time_envelope(t):
#    return np.exp(-2.*np.log(2.)*(t-xf-delay)**2/FWHMtinI**2)
#
#def carrier(t):
#    return np.cos(t + phi)
#
#def space_envelope(y,z):
#    return np.exp(-((y-Lsim[1]/2.0)**2.0+(z-Lsim[2]/2.0)**2.0)/w0**2.0) #gauss

#def my_filter(particles):
#    return particles.id>0

def n0_electron(x,y,z):
#        if ((9.8*l0 < x < 10*l0) and (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2) < 3.15*l0)):
#        if (np.sqrt((x-Lsim[0]/2)**2+(y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2) < 3.15*l0):
        if (np.sqrt((x-Lsim[0]/2)**2+(y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2) < 0.02*l0):
#        if ((4.9*l0 < x < 5.1*l0) and (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2) < 0.1*l0)):
#        if ((np.abs(x-Lsim[0]/2) < 2*l0) and (np.abs(y-Lsim[1]/2) < 2*l0) and (np.abs(z-Lsim[2]/2) < 2*l0)):
                return n0
        else:
                return 0.

# LP
'''
def By_profile(y,z,t):
    return 0

def Bz_profile(y,z,t):
    return a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * np.cos(t + phi0) * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)
'''


def By_profile(y,z,t):
    factor = (z-Lsim[2]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (15)
    return -a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi01)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi01)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)

def Bz_profile(y,z,t):
    factor = (y-Lsim[1]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (15)
    return a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi01)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi01)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi01)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi01)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)

def By_profile_counter(y,z,t):
    factor = (z-Lsim[2]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (13)
    return -a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi02)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi02)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)

def Bz_profile_counter(y,z,t):
    factor = (y-Lsim[1]/2)/np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)
    # 5. rad (15)
    return a0 * np.exp(-((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0**2) * (eps * (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02) + eps**3 * ((np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/2 +(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/2 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/4) + eps**5 * (3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0) * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**3 * np.cos(t+phi02)/8 + 3*(np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**5 * np.cos(t+phi02)/16 - (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**7 * np.cos(t+phi02)/4 + (np.sqrt((y-Lsim[1]/2)**2+(z-Lsim[2]/2)**2)/w0)**9 * np.cos(t+phi02)/32)) * factor  * np.exp(-2.*np.log(2.)*(t-delay)**2/FWHMtinI**2)
    
    
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
    #extra_envelope          = ,
    #keep_n_strongest_modes = 200, # def
    angle = 0,
    number_of_processes    = 64,
#    file                   = './LaserOffset0.h5',
#    fft_time_window        = Tsim/2 # aby se impuls dostal za polovinu simulacni oblasti. Po te se zacen impuls opakovat. Funguje urcite pro 1T.
#    fft_time_window        = 1.0*Tsim # aby se impuls dostal za polovinu simulacni oblasti. Po te se zacen impuls opakovat. Funguje urcite pro 1T.
    fft_time_window        = 1*Tsim # aby se impuls dostal za polovinu simulacni oblasti. Po te se zacen impuls opakovat. Funguje urcite pro 1T.
)

LaserOffset(
    box_side               = "xmax",
    space_time_profile     = [ By_profile_counter, Bz_profile_counter ],
    offset                 = xf,
    #extra_envelope          = ,
    #keep_n_strongest_modes = 200, # def
    angle = 0,
    number_of_processes    = 64,
#    file                   = './LaserOffset1.h5',
#    fft_time_window        = Tsim/2 # aby se impuls dostal za polovinu simulacni oblasti. Po te se zacen impuls opakovat. Funguje urcite pro 1T.
    fft_time_window        = 1.0*Tsim # aby se impuls dostal za polovinu simulacni oblasti. Po te se zacen impuls opakovat. Funguje urcite pro 1T.
)

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 1, #10  #128
    c_part_max = 1.,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.],
    pusher = "vay",
    radiation_model = "mc",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,
#    time_frozen = 4.1*t0, #5.46*t0,
    boundary_conditions = [
        ["remove","remove"],
        ["remove","remove"],
        ["remove","remove"],
    ]
    #track_every = 2,
    #track_flush_every = 100,
    #is_test = False
)

Species(
    name = "bwelectron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0, #10  #128
    c_part_max = 1.,
    mass = 1.0,
    charge = -1.0,
    charge_density = 0.0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.],
    pusher = "vay",
    radiation_model = "mc",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,
#    time_frozen = 4.1*t0, #5.46*t0,
    boundary_conditions = [
        ["remove","remove"],
        ["remove","remove"],
        ["remove","remove"],
    ]
    #track_every = 2,
    #track_flush_every = 100,
    #is_test = False
)

Species(
    name = "positron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 1.0,
    mass = 1.0,
    charge = 1.0,
    charge_density = 0.0,
    mean_velocity = [0.0, 0.0, 0.0],
    temperature = [0.],
    pusher = "vay",
    radiation_model = "mc",
    radiation_photon_species = "photon",
    radiation_photon_sampling = 1,
    radiation_photon_gamma_threshold = 2,
    boundary_conditions = [
        ["remove","remove"],
        ["remove","remove"],
        ["remove","remove"],
    ]
    #track_every = 2,
    #track_flush_every = 100,
    #is_test = False
)

Species(
    name = "photon",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 0,
    c_part_max = 20.0,
    mass = 0,
    charge = 0.,
    number_density = 0.0,
    mean_velocity = [0.0 ,0.0, 0.0],
    temperature = [0.],
    pusher = "norm",
    multiphoton_Breit_Wheeler = ["bwelectron","positron"],
    multiphoton_Breit_Wheeler_sampling = [1,1],
    boundary_conditions = [
        ["remove","remove"],
        ["remove","remove"],
        ["remove","remove"],
    ]
    #track_every = 2,
    #track_flush_every = 100,
    #is_test = False
)

RadiationReaction(
#    minimum_chi_continuous = 1e-3,
#    minimum_chi_discontinuous = 1e-2,
#    table_path = "/media/jirkama1/Data/smilei/smilei-new/nics_256x256_1e-4_2e3",
#    table_path = "/gpfs/scratch/work/user/mjirka/mjirka/smilei/smilei-v4.6/nics_256x256_1e-4_2e3",
)

MultiphotonBreitWheeler(
)

DiagScalar(
    every = [1.*rest,25.*rest,1],
)
'''
DiagFields(
#    every = [7.*rest,9.35*rest,4*rest], # pred srazkou
    every = [8.75*rest,9.35*rest,0.125*rest], # celkovy pohled
    fields = ['Ex','Ey','Ez','Bx','By','Bz'],
#    fields = ['By'],
#    subgrid = s_[96:224,96:224,96:224]
)
'''
'''
DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [5*l0, 5*l0, 5*l0],
    vector = [1., 0., 0.],
    direction = "backward",
    deposited_quantity = "weight_chi",
    species = ["electron"],
    axes = [["y", 0, 10.*l0, 640],
            ["z", 0, 10.*l0, 640]],
    every = [1*rest,25*rest,1*rest]
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [5*l0, 5*l0, 5*l0],
    vector = [1., 0., 0.],
    direction = "backward",
    deposited_quantity = lambda p: p.chi,
    species = ["electron"],
    axes = [["y", 0, 10.*l0, 640],
            ["z", 0, 10.*l0, 640]],
    every = [1*rest,25*rest,1*rest]
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [0.5*l0, 5*l0, 5*l0],
    vector = [1., 0., 0.],
    direction = "backward",
    deposited_quantity = "weight",
    species = ["electron"],
    axes = [["y", 0, 10.*l0, 640],
            ["z", 0, 10.*l0, 640]],
    every = [1*rest,25*rest,1*rest]
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [0.5*l0, 5*l0, 5*l0],
    vector = [1., 0., 0.],
    direction = "backward",
    deposited_quantity = "weight_ekin",
    species = ["electron"],
    axes = [["y", 0, 10.*l0, 640],
            ["z", 0, 10.*l0, 640]],
    every = [1*rest,25*rest,1*rest]
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [5*l0, 5*l0, 5*l0],
    vector = [1., 0., 0.],
    direction = "backward",
    deposited_quantity = lambda p: p.px,
    species = ["electron"],
    axes = [["y", 0, 10.*l0, 640],
            ["z", 0, 10.*l0, 640]],
    every = [1*rest,25*rest,1*rest]
)

DiagScreen(
    #name = "my screen",
    shape = "plane",
    point = [0.5*l0, 5*l0, 5*l0],
    vector = [1., 0., 0.],
    direction = "backward",
    deposited_quantity = lambda p: p.px,
    species = ["electron"],
    axes = [["y", 0, 10.*l0, 640],
            ["z", 0, 10.*l0, 640]],
    every = [1*rest,25*rest,1*rest]
)
'''
'''
DiagTrackParticles(
    species = "electron",
#    every = [5*rest,10*rest,0.1*rest],
    every = [0*rest,12*rest,0.1*rest],
#    flush_every = 100,
#    filter = my_filter,
    attributes = ["x","y","z","px","py","pz","w","chi","Ex","Ey","Ez","Bx","By","Bz"],
##    attributes = ["x","y","z","px","py","pz","w"]
)

DiagTrackParticles(
    species = "positron",
#    every = [5*rest,10*rest,0.1*rest],
    every = [0*rest,12*rest,0.1*rest],
##    flush_every = 100,
##    filter = my_filter,
##    attributes = ["x","px","chi"],
    attributes = ["x","y","z","px","py","pz","w","chi","Ex","Ey","Ez","Bx","By","Bz"],
)

DiagTrackParticles(
    species = "bwelectron",
#    every = [5*rest,10*rest,0.1*rest],
    every = [0*rest,12*rest,0.1*rest],
##    flush_every = 100,
##    filter = my_filter,
##    attributes = ["x","px","chi"],
    attributes = ["x","y","z","px","py","pz","w","chi","Ex","Ey","Ez","Bx","By","Bz"],
)
'''
'''
DiagParticleBinning(
    deposited_quantity = "weight",
    every = [1*rest,25*rest,1*rest],
    time_average = 1,
    species = ["electron"],
    axes = [ ["ekin",    0.05,    2000.,   4000, ] ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = [1*rest,25*rest,1*rest],
    time_average = 1,
    species = ["photon"],
    axes = [ ["ekin",    0.05,    2000.,   4000, ] ]
)

DiagParticleBinning(
    deposited_quantity = "weight",
    every = [1*rest,25*rest,1*rest],
    time_average = 1,
    species = ["positron"],
    axes = [ ["ekin",    0.05,    2000.,   4000, ] ]
)
'''
