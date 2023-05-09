#!/usr/bin/python3
#Author: James Glazar
#Date: 8 July 2020
#Purpose: HW1 Solution -- Find Lennard-Jones energy as function of lattice
#parameter for triangular 2D lattice.

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# parameters from lecture notes
EPSILON = 0.010323  
SIGMA = 3.405
RCUT = 7.5
R0 = 7.0
A = -6.8102128*10**-3
B = -5.564087*10**-3

# Calculate pairwise Lennard-Jones energy using smooth correction
def LJ(r):  # r = distance from central atom
  if r < R0:
    return 4*EPSILON*( (SIGMA/r)**12 - (SIGMA/r)**6 )
  elif r < RCUT:
    return A*(r-RCUT)**3 + B*(r-RCUT)**2
  else:
    return 0

def Et(l, n):  # l = list of a,E1/2/3  n = neighbors to include
  if n == 3:
    return l['E1'] + l['E2'] + l['E3']
  elif n == 2:
    return l['E1'] + l['E2']
  elif n == 1:
    return l['E1']
  else:
    return 0

x = np.arange(3.5, 4.0, 0.0001)
eng_df = pd.DataFrame(
  {'a':x, 
   'E1':[6.*LJ(i) for i in x], 
   'E2':[6.*LJ(i*np.sqrt(3)) for i in x],
   'E3':[6.*LJ(i*2) for i in x]
  }
)
eng_df['Et_3'] = eng_df.apply(Et, args=[3], axis=1)
eng_df['Et_2'] = eng_df.apply(Et, args=[2], axis=1)
eng_df['Et_1'] = eng_df.apply(Et, args=[1], axis=1)


print('---Including up to 3rd neighbors---')
print('Equilibrium Lattice Parameter and its Energy: ',
eng_df[ eng_df['Et_3'] == eng_df['Et_3'].min()].a.item(),
eng_df[ eng_df['Et_3'] == eng_df['Et_3'].min()].Et_3.item() 
 )

print('---Including up to 2nd neighbors---')
print('Equilibrium Lattice Parameter and its Energy: ',
eng_df[ eng_df['Et_2'] == eng_df['Et_2'].min()].a.item(),
eng_df[ eng_df['Et_2'] == eng_df['Et_2'].min()].Et_2.item()
 )

print('---Including up to 1st neighbors---')
print('Equilibrium Lattice Parameter and its Energy: ',
eng_df[ eng_df['Et_1'] == eng_df['Et_1'].min()].a.item(),
eng_df[ eng_df['Et_1'] == eng_df['Et_1'].min()].Et_1.item()
 )

plt.xlabel('Lattice Parameter [Ang]')
plt.ylabel('Energy [eV]')
plt.plot('a', 'Et_1', data=eng_df, color='paleturquoise', linewidth=3.0, label='up to 1st neighbors')
plt.plot('a', 'Et_2', data=eng_df, color='mediumturquoise', linewidth=3.0, label='up to 2nd neighbors')
plt.plot('a', 'Et_3', data=eng_df, color='teal', linewidth=3.0, label='up to 3rd neighbors')

plt.legend()
plt.show()
