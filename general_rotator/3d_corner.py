"""
This example attempts to use general rotations to go beyind the tst3d_20_y_em_propagation.py exmaple from
Smilei benchmarks. It shows how to define a general B(x,y,z,t) field (also non-paraxial cases considered)
and then to transform it. The transformation uses general translation and rotation to project given
components onto the boundary planes of the Cartesian box.
"""
from math import pi, cos, sin, tan, asin, atan
import numpy as np


l0 = 2.0*np.pi              # laser wavelength
t0 = l0                       # optical cicle
Lsim = [10.*l0,8.*l0,5.*l0]  # length of the simulation
Tsim = 18.*t0                 # duration of the simulation
resx = 12.                    # nb of cells in one laser wavelength
rest = 22.                    # nb of timesteps in one optical cycle 


# ang = [-pi/7., pi/6.] # angle as defined by the smilei tutorial
ang = [0.,0.] # angle as defined by the smilei tutorial
focus = [0.5*Lsim[0], 0.5*Lsim[1], 0.5*Lsim[2]] # the position of the focus in the Smilei coordinates = offset in the 

### General rotation matrices
# Rotation matrices
def RotM(axis,angle):
    """
    Rotation matrices in standard convention.
    Args:
        axis : specify rotation axis string in 'x', 'y', 'z'.
        angle (float): Rotation angle.

    Returns:
        3x3 np-array: rotation matrix.
    """
    if (axis=='x'):
        return np.asarray([[1, 0            ,  0            ],
                           [0, np.cos(angle), -np.sin(angle)],
                           [0, np.sin(angle),  np.cos(angle)]])
    elif (axis=='y'):
        return np.asarray([[ np.cos(angle), 0,  np.sin(angle) ],
                           [ 0            , 1,  0             ],
                           [-np.sin(angle), 0,  np.cos(angle) ]])
    elif (axis=='z'):
        return np.asarray([[ np.cos(angle), -np.sin(angle),  0 ],
                           [ np.sin(angle),  np.cos(angle),  0 ],
                           [ 0            ,  0            ,  1 ]]) 

### Transformation function
def transform_vector_field(R,offset,vector,selection=(0,1,2)):
    """
    Tranform a vector field V=`vector` according to the general rotation V'(x')=R*V(RT*x'+offset),
    coordinate notation follows the documents.
    R - rotation matrix
    offset - offset vector
    vector - vector field represented by a function with 3 components taking (x,t) arguments, x is 3-component spatial, t is time
    selection - select output projections, labeled by components 0, 1, 2, allows permutations and dimension reduction

    output - vector field represented by a function inlcuding selection
    """
    # def vector_transformed(x_):
    #     return np.asarray([(R@vector((np.transpose(R)@x_) + offset))[idx] for idx in selection])
    def vector_field_transformed(x_,t_):
        return np.asarray([(np.transpose(R)@vector(R@(x_ - offset),t_))[idx] for idx in selection])
    return vector_field_transformed

### General Gaussian beam represented by B-field
waist = l0
a0    = 1.
omega = 1.

fwhm = t0*0.6
def time_envelope(t):
    sigma = (0.5*fwhm)**2/np.log(2.0)
    return np.exp( -(t)**2 / sigma )

# derived quatities
k0 = 1.
zR = omega * waist**2/2.

def B_Gauss_lin(x, t):

    r2 = x[1]**2 + x[2]**2
    w = l0 * np.sqrt(1.0 + (x[0]/zR)**2)
    invR = x[0] / (x[0]**2 + zR**2)
    gouy = np.arctan(x[0] / zR)
    
    phase = (omega*t - k0*x - 0.5 * r2 * invR  + gouy)
    envelope = ((waist/w)*np.exp(-r2/w**2)*time_envelope(t - x[0])
    )

    return a0 * envelope * np.cos(phase)
    
def Bfield(x,t):
    return np.asarray([
        0,
        B_Gauss_lin(x, t),
        0        
        ])

xfix = 0.
def Bfield_xplane(xfix):
    def B1(y_,z_,t_):
        return transform_vector_field(R('z',ang[1])@R('y',ang[0]),Bfield,selection=(1,))(xfix,y_,z_,t_)
    def B2(y_,z_,t_):
        return transform_vector_field(R('z',ang[1])@R('y',ang[0]),Bfield,selection=(2,))(xfix,y_,z_,t_)
    return [B1, B2]
    
    

Main(
    geometry = "3Dcartesian",
    
    interpolation_order = 2 ,
    
    cell_length = [l0/resx,l0/resx,l0/resx],
    grid_length  = Lsim,
    
    number_of_patches = [ 4,4,4 ],
    
    timestep = t0/rest,
    simulation_time = Tsim,
    
    EM_boundary_conditions = [ ['PML'] ],
)

Laser(
    box_side = "xmin",
    space_time_profile = Bfield_xplane(0.) #[ By_profile, Bz_profile ],
)

# time_envelope = tgaussian(fwhm=t0*6, center=t0*9)

# LaserGaussian3D(
#     box_side        = "xmin",
#     a0              = 1.,
#     omega           = 1.,
#     focus           = focus,
#     waist           = l0,
#     incidence_angle = ang,
#     time_envelope   = time_envelope,
# )

# LaserGaussian3D(
#     box_side        = "ymin",
#     a0              = 1.,
#     omega           = 1.,
#     focus           = focus,
#     waist           = l0,
#     incidence_angle = [ang[0], pi/2.-ang[1]],
#     time_envelope   = time_envelope
# )

# LaserGaussian3D(
#     box_side        = "zmin",
#     a0              = 1.,
#     omega           = 1.,
#     focus           = focus,
#     waist           = l0,
#     incidence_angle = [-asin(cos(ang[0])*sin(ang[1])), -atan(cos(ang[1])/tan(ang[0]))],
#     time_envelope   = time_envelope,
#     polarization_phi = -pi/2,
# )


globalEvery = int(rest)

DiagScalar(
    every=globalEvery
)

DiagFields(
    every = globalEvery,
    fields = ['Ex','Ey','Ez','Bx','By','Bz']
)
