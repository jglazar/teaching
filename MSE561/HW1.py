#Author: James Glazar
#Date: 8 July 2020
#Purpose: HW1 Solution -- Find Lennard-Jones energy as function of lattice
#parameter for triangular 2D lattice.

import matplotlib.pyplot as plt

EPSILON = 0.010323
SIGMA = 3.405
RCUT = 7.5
R0 = 7.0
A = -6.8102128*10**-3
B = -5.564087*10**-3

def eng(r):
  if r < R0:
    return 4*EPSILON*( (SIGMA/r)**12 - (SIGMA/r)**6 )
  elif r < RCUT:
    return A*(r-RCUT)**3 + B*(r-RCUT)**2
  else:
    return 0

x = [0.0001*i for i in range(35*10**3,40*10**3+1,1)]
y = []
min_a = 0
min_energy = 0

for a in x:
  energy = 6 * ( eng(a) + eng(a*3**0.5) + eng(a*2) )
  if energy < min_energy:
    min_energy = energy
    min_a = a
  y.append(energy)

print('Minimum Energy: ', min_energy)
print('Equilibrium Lattice Parameter: ', min_a)

plt.xlabel('Lattice Parameter [Ang]')
plt.ylabel('Energy [eV]')

plt.plot(x, y, color='teal', linewidth=3.0)
plt.show()
