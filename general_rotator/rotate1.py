import numpy as np
import units
import mynumerics as mn


# 1 Define rotation matrix
a = 5.
# T = Rz*Rx*Rz

# 2 Construct the field


def GaussianBeamEfield(r,z,t,E0,w0,tFWHM,lambd, comoving = True):
  # Gaussian beam in the comoving frame with c
  if (not(comoving)): t = t - z/units.c_light
  omega0 = mn.ConvertPhoton(lambd, 'lambdaSI', 'omegaSI')
  k0 = 2.0*np.pi/lambd
  zR = np.pi*w0**2/lambd
  w=w0*np.sqrt(1.0+(z/zR)**2)
  phase_Gouy = np.arctan(z/zR)
  phase_curv = mn.GaussianBeamCurvaturePhase(r,z,k0,zR)
  return E0*(w0/w)*np.exp(-(r/w)**2)*np.exp(-(2.0*np.log(2.0)*t/tFWHM)**2)*np.cos(omega0*t -phase_curv + phase_Gouy)

# 3 Fourier transform it (k-space) and (omega-space)

# k-space
# Bkykz = np.fft.fft2

# 4 Apply amplitude corrections

# 5 Fourier it back

# 6 Project onto the plane, keeping only the 2 planar components