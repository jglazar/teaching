#!/usr/bin/python3
#Author: James Glazar
#Date: 6 August 2020
#Purpose: HW4 Solution -- relax square 20 x 20 block of atoms
#interacting via Lennard-Jones potential using Monte Carlo

import matplotlib.pyplot as plt
from random import uniform
from math import exp

''' Parameters '''
EPSILON = 0.010323  
SIGMA = 3.405
INIT_R = SIGMA*2**(1.0/6.0)
RCUT = 7.5
R0 = 7.0
A = -6.8102128*10**-3
B = -5.564087*10**-3
RDF_MAX = 12.0
AREA = (19.*INIT_R)**2
RHO = 400. / AREA
PI = 3.1415926536
MAX_STEP = 50
DISP = 0.5
KB = 8.617*10**-5
TEMP = 5.

''' End of parameters '''

''' Useful functions '''
# Calculate pairwise Lennard-Jones energy using smooth correction
def LJ(r):  # r = distance from central atom
  if r < R0:
    return 4*EPSILON*( (SIGMA/r)**12 - (SIGMA/r)**6 )
  elif r < RCUT:
    return A*(r-RCUT)**3 + B*(r-RCUT)**2
  else:
    return 0

# Calculate derivative of above potential
def dLJ(r):  # r = distance from central atom
  if r < R0:
    return 4*EPSILON*( -12.0/r *(SIGMA/r)**12 + 6.0/r * (SIGMA/r)**6 )
  elif r < RCUT:
    return 3*A*(r-RCUT)**2 + 2*B*(r-RCUT)
  else:
    return 0

# Calculate distance between two coordinates
def dist(a, b):
  return ((a[0] - b[0])**2 + (a[1] - b[1])**2)**0.5

# Get lists of coords of nearest neighbors up to cutoff
def NN(c, cut):  # c = (x,y) coordinate, cut = cutoff distance
  return [i for i in coords if i!=c and dist(c,i) < cut]  

# plot the rdf the formal way
def rdf(l):  # c = list of coords
  del_r = RDF_MAX/100.
  x = [i*del_r for i in range(1,101)]
  g = [0 for _ in range(1,101)]
  for a in l:
    f = [0 for _ in range(1,101)]
    dists = [dist(a, i) for i in NN(a, RDF_MAX)]
    for ind, d in enumerate(x):
      count = sum( [1 for i in dists if d-del_r<i<=d] )
      f[ind] = count / (2*PI*d*del_r)
    g = [sum(i) for i in zip(f,g)]
  return [i/(400.*RHO) for i in g]

# accept or reject proposed change?
def accept(e):  # e = energy difference (final - initial)
  if e <= 0 or uniform(0,1) <= exp(-e/(KB*TEMP)):
    return True
  else:
    return False

'''End of useful functions '''

coords = [ (INIT_R*float(i), INIT_R*float(j)) for i in range(0,20) for j in range(0,20) ]
center = coords[10*20 + 10]

# plot the initial grid of atoms
plt.scatter([i[0] for i in coords], [i[1] for i in coords])
plt.xlabel('x coordinate [Ang]')
plt.ylabel('y coordinate [Ang]')
plt.show()

# plot the initial rdf the formal way
plt.plot([i*RDF_MAX/100. for i in range(1,101)], rdf(coords))
plt.xlabel('Distance [Ang]')
plt.ylabel('Radial distribution function [# atoms]')
plt.show()

acc_counter = 0
step = 1
while step <= MAX_STEP:
  for ind, c in enumerate(coords):
    eng_i = 0
    for n in NN(c, RCUT):
      eng_i += LJ(dist(c,n))      
    c_f = ( c[0] + uniform(-DISP,DISP), c[1] + uniform(-DISP,DISP) )
    eng_f = 0
    for n in [i for i in NN(c_f, RCUT) if i!=c]:
      eng_f += LJ(dist(c_f,n)) 
    if accept(eng_f - eng_i):
      acc_counter += 1
      coords[ind] = c_f
  step += 1

# print fraction of moves that were accepted
print(acc_counter / (400*MAX_STEP))

print(coords)

# plot the final grid of atoms
plt.scatter([i[0] for i in coords], [i[1] for i in coords])
plt.xlabel('x coordinate [Ang]')
plt.ylabel('y coordinate [Ang]')
plt.show()

# plot the final rdf the formal way
plt.plot([i*RDF_MAX/100. for i in range(1,101)], rdf(coords))
plt.xlabel('Distance [Ang]')
plt.ylabel('Radial distribution function [# atoms]')
plt.show()
