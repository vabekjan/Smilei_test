# _____________________________________________________________________________
#
# Collision electron beam with a counter-propagating plane wave
#
# _____________________________________________________________________________

from math import pi, sqrt

# _____________________________________________________________________________
# Main parameters

c = 299792458
lambdar = 1e-6
wr = 2*pi*c/lambdar

l0 = 2.0*pi              # laser wavelength
t0 = l0                       # optical cicle
Lx = 30*l0                    # Longitudinal length
Ly = 4*l0                     # Transverse length

n0 = 1e-5                     # particle density

Tsim = 60.*t0                 # duration of the simulation
resx = 64.                    # nb of cells in one laser wavelength
resy = 64.

dx = l0/resx                            # space step
dy = l0/resy
dt  = 0.95 * 1./(1./dx + 1./dy)         # timestep (0.95 x CFL)

start = 0                               # Laser start
fwhm = 10*t0                            # Gaussian time fwhm
duration = 50*t0                        # Laser duration
center = duration*0.5                   # Laser profile center
order = 4                               # Laser order

gamma = 1000./0.511                     # Electron beam gamma factor
v = sqrt(1 - 1./gamma**2)          # electron beam initial velocity

# _____________________________________________________________________________
# Density profile for inital location of the particles and photons

def n0_electron(x,y):
        if ((0.97*Lx < x < 0.99*Lx)and(0.48*Ly < y < 0.52*Ly)):
                return n0
        else:
                return 0.

# ______________________________________________________________________________
# Namelist Main

Main(
    geometry = "2Dcartesian",

    interpolation_order = 4 ,

    cell_length = [dx,dy],
    grid_length  = [Lx,Ly],

    number_of_patches = [4,4],

    timestep = dt,
    simulation_time = Tsim,

    EM_boundary_conditions = [
        ["silver-muller", "silver-muller"],
        ["periodic", "periodic"]
    ],

    reference_angular_frequency_SI = wr,

    random_seed = 0

)

# ______________________________________________________________________________
# Namelist Laser

LaserGaussian2D(
    box_side         = "xmin",
    a0              = 100.,
    omega           = 1.,
    focus           = [0.5*Lx, 0.5*Ly],
    waist           = 1E9,
    incidence_angle = 0.,
    polarization_phi = 0.,
    ellipticity     = 0,
    time_envelope  = tgaussian(start=start,duration=duration,
                               fwhm=fwhm,
                               center=center,
                               order=order)
)

# ______________________________________________________________________________
# Namelist Species

Species(
    name = "electron",
    position_initialization = "random",
    momentum_initialization = "cold",
    particles_per_cell = 32,
    c_part_max = 1.,
    mass = 1.0,
    charge = -1.0,
    charge_density = n0_electron,
    mean_velocity = [-v, 0.0, 0.0],
    temperature = [0.],
    pusher = "vay",
   radiation_model = "Monte-Carlo",
#    radiation_photon_species = "photon",
#    radiation_photon_sampling = 1,
#    radiation_photon_gamma_threshold = 2,
    time_frozen = 29*t0,
    boundary_conditions = [["remove","remove"],
                            ["periodic","periodic"]]
    #track_every = 2,
    #track_flush_every = 100,
    #is_test = False
)



# ______________________________________________________________________________
# Namelist configuration of the radiation reaction models

RadiationReaction(
   minimum_chi_continuous = 1e-3,
   minimum_chi_discontinuous = 1e-3,
#    table_path = "/gpfshome/mds/staff/mlobet/smilei/databases/"
)

# ______________________________________________________________________________
# Namelist scalar diagnostics

DiagScalar(
    every = 10,
    vars=['Uelm','Ukin','Utot','Uexp','Ubal','Urad',
          'Ukin_electron',
          'Ntot_electron']
)

# ______________________________________________________________________________
# Namelist Fields

DiagFields(
    every = 500,
    fields = ['Ex','Ey','Ez','By','Bz']
)

# ______________________________________________________________________________
# Namelist particle binning diagnostics

# 1. 2D grid of the weight

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 1000]
    ]
)


# 2. 2D grid of the weight x normalized kinetic energy (gamma factor - 1)

DiagParticleBinning(
    deposited_quantity = "weight_ekin",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 1000]
    ]
)


# 3. Quantum parameter x weight
# Quantum parameter particle binning can be switched off if no radiation losses

DiagParticleBinning(
    deposited_quantity = "weight_chi",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["x", 0., Lx, 1000],
        ["y", 0., Ly, 1000]
    ]
)


# 4. Energy distributions

DiagParticleBinning(
    deposited_quantity = "weight",
    every = 500,
    time_average = 1,
    species = ["electron"],
    axes = [
        ["gamma", 10., 3000, 200,'logscale']
    ]
)
