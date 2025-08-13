from math import pi, cos, sin

l0 = 2.0*pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [7.*l0,10.*l0,10.*l0]  # length of the simulation
Tsim = 8.*t0                 # duration of the simulation
resx = 16.                    # nb of cells in one laser wavelength
rest = 30.                    # nb of timesteps in one optical cycle 

angle1 = 0.2
angle2 = 0.1

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 4 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 8,8,8 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [ ['silver-muller'] ],
    EM_boundary_conditions_k = [
        [cos(angle1)*cos(angle2), sin(angle2), -sin(angle1)*cos(angle2)],
        [-cos(angle1)*cos(angle2), -sin(angle2), sin(angle1)*cos(angle2)],
        [0., 1., 0.],[0., -1., 0.],
        [0., 0., 1.],[0., 0., -1.],],
)

LaserGaussian3D(
    a0              = 1.,
    omega           = 1.,
    focus           = [0.9*Lsim[0], 0.6*Lsim[1], 0.3*Lsim[2]],
    waist           = l0,
    incidence_angle = [angle1, angle2],
#    time_envelope   = tgaussian()
)

Species(
	name = "eon",
	position_initialization = "regular",
	momentum_initialization = "cold",
	particles_per_cell = 0.0001,
	mass = 1.0,
	charge = -1.0,
	number_density = 0.0001,
	mean_velocity = [0.,0.,0.],
	boundary_conditions = [
		["remove", "remove"],
	],
	is_test = True
)

Vectorization(
   mode="off",
)

globalEvery = int(rest)

DiagScalar(
    every=globalEvery
)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez']
)
from numpy import s_
DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez'],
    subgrid = s_[4:100:3, 5:400:10, 6:300:80]
)
DiagFields(
    every = 2*globalEvery,
    fields = ['Ex','Ey','Ez'],
    time_average = globalEvery
)

DiagProbe(
    every = 10,
    origin = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [30],
    origin = [0.1*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]],
    corners = [[0.9*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]]],
    fields = []
)

DiagProbe(
    every = 100,
    number = [10, 10],
    origin = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    corners = [
        [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
    ],
    fields = []
)

DiagProbe(
    every = 100,
    number = [4, 4, 4],
    origin = [0.1*Lsim[0], 0.*Lsim[1], 0.5*Lsim[2]],
    corners = [
        [0.9*Lsim[0], 0. *Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0.9*Lsim[1], 0.5*Lsim[2]],
        [0.1*Lsim[0], 0. *Lsim[1], 0.9*Lsim[2]],
    ],
    fields = []
)

DiagTrackParticles(
    species = "eon",
    every = 10,
	attributes = ["x", "y", "z", "px", "py", "pz", "Ex", "Ey", "Ez", "Bx", "By", "Bz"]
)
