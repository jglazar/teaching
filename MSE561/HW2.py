#!/usr/bin/python3
#Author: James Glazar
#Date: 13 July 2020
#Purpose: HW2 Solution -- Set up square 20 x 20 block of atoms
#interacting via Lennard-Jones potential

# Note -- this includes approximate and fully formal methods
# to compute the RDF and stress tensor. They show similar behavior. 

import matplotlib.pyplot as plt

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
''' End of parameters '''

''' Useful functions '''
# Calculate pairwise Lennard-Jones energy using smooth correction
def LJ(r):  # r = distance from central atom
  if r < R0:
    return 4*EPSILON*( (SIGMA/r)**12 - (SIGMA/r)**6 )
  if r < RCUT:
    return A*(r-RCUT)**3 + B*(r-RCUT)**2
  return 0

#Calculate derivative of above potential
def dLJ(r):  # r = distance from central atom
  if r < R0:
    return 4*EPSILON*( -12.0/r *(SIGMA/r)**12 + 6.0/r * (SIGMA/r)**6 )
  if r < RCUT:
    return 3*A*(r-RCUT)**2 + 2*B*(r-RCUT)
  return 0

# Calculate distance between two coordinates
def dist(a, b):
  return ((a[0] - b[0])**2 + (a[1] - b[1])**2)**0.5

# Get lists of coords of nearest neighbors up to n'th shell
def NN(c, n):  # c = (x,y) coordinate, n = # nearest neighbors
  out = [ ]
  neighbors = [i for i in coords if i!=c]
  dists = [dist(c, i) for i in neighbors]
  for i in range(n):
    closest = [i for i in neighbors if dist(c,i)==min(dists)]
    out.append(closest)
    neighbors = [i for i in neighbors if i not in closest]
    dists = [dist(c,i) for i in neighbors]
  return out

# Get list of coords for nearest neighbors up to cutoff distance
def neigh(c, cut):  # c = (x,y) coordinate, cut = cutoff distance
  return [i for i in coords if i!=c and dist(c,i) < cut]

# plot the rdf the formal way
def rdf(l):  # c = list of coords
  del_r = RDF_MAX/100.
  x = [i*del_r for i in range(1,101)]
  g = [0 for _ in range(1,101)]
  for a in l:
    f = [0 for _ in range(1,101)]
    dists = [dist(a, i) for i in neigh(a, RDF_MAX)]
    for ind, d in enumerate(x):
      count = sum( [1 for i in dists if d-del_r<i<=d] )
      f[ind] = count / (2*PI*d*del_r)
    g = [sum(i) for i in zip(f,g)]
  return [i/(400.*RHO) for i in g]

# calculate stress tensor the formal way (average over all atoms)
def stress(l):  # l = list of coordinates
  stress_xx = 0
  stress_xy = 0
  stress_yx = 0
  stress_yy = 0
  for a in l:
    for n in neigh(a, RCUT):
      d = dist(a, n)
      x = a[0] - n[0]
      y = a[1] - n[1]
      stress_xx += dLJ(d) * x**2 / d
      stress_xy += dLJ(d) * x*y / d
      stress_yx += dLJ(d) * y*x / d
      stress_yy += dLJ(d) * y**2 / d
  stress_xx /= AREA
  stress_xy /= AREA
  stress_yx /= AREA
  stress_yy /= AREA
  return [stress_xx, stress_xy, stress_yx, stress_yy]


'''End of useful functions '''

coords = [ (INIT_R*float(i), INIT_R*float(j)) for i in range(0,20) for j in range(0,20) ]
center = coords[10*20 + 10]

# get neighbors of middle atom (10,10)
neighbors = NN( center, 3 )

# plot the initial grid of atoms
plt.scatter([i[0] for i in coords], [i[1] for i in coords])
plt.xlabel('x coordinate [Ang]')
plt.ylabel('y coordinate [Ang]')
plt.show()

# plot the radial distribution function for the central atom using a histogram
dists = [dist(center, i) for i in [j for j in coords if j!=center]]
plt.hist(dists, range=(0,RDF_MAX), rwidth=0.5, bins=50)
plt.xlabel('Distance [Ang]')
plt.ylabel('Radial distribution function [# atoms]')
plt.show()

# plot the initial rdf the formal way
plt.plot([i*RDF_MAX/100. for i in range(1,101)], rdf(coords))
plt.xlabel('Distance [Ang]')
plt.ylabel('Radial distribution function [# atoms]')
plt.show()

#center = coords[9]
#print(center)

# calculate the stress tensor for the central atom (incorrect scaling, but that's OK)
stress_xx = 0
stress_xy = 0
stress_yx = 0
stress_yy = 0

for a in [i for i in coords if i!=center]:
  d = dist(center, a)    
  stress_xx += dLJ(d) * (center[0] - a[0])**2 / d
  stress_xy += dLJ(d) * (center[0] - a[0])*(center[1] - a[1]) / d
  stress_yx += dLJ(d) * (center[1] - a[1])*(center[0] - a[0]) / d
  stress_yy += dLJ(d) * (center[1] - a[1])**2 / d

  if dLJ(d) != 0:
    print(d, dLJ(d))

print('stress xx, xy, yx, yy: ', stress_xx, stress_xy, stress_yx, stress_yy)

# calculate the stress tensor the formal way
s_all = stress(coords)
print('formal stress xx, xy, yx, yy: ', s_all[0], s_all[1], s_all[2], s_all[3])
