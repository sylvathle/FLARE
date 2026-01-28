import numpy as np

mev2joules = 1.60218e-13
c = 299792458 #m/s
mass_proton = 1.67262192e-27*c**2/mev2joules

def rigidity(kE,m,q):
  return float(np.sqrt(kE**2+2*m*kE)/q * 1e-3)

def rigidity2Energy(r,m,q):
  return float(-m+np.sqrt(m**2+(q*r/1e-3)**2))

def beta(kE,m):
  p = np.sqrt(kE*(kE+2*m))
  return float(p / np.sqrt(p**2+m**2))

  # prompt: Write a function that converts kinetic energy to momentum

def kinetic_energy_to_momentum(kinetic_energy, mass):
  momentum = mass * np.sqrt(2 * kinetic_energy / mass)
  return momentum



