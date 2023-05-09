
#!/usr/bin/python3
#Author: James Glazar
#Date: 13 July 2020
#Purpose: HW3 Solution -- relax square 20 x 20 block of atoms
#interacting via Lennard-Jones potential

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
MAX_STEP = 20000
LAMBDA = 1.0

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

# Total force vector acting on given atom
def force(c):  # c = (x,y) coordinate  
  out = [0., 0.]
  for n in NN(c, RCUT):  # sum forces over neighbors within potential cutoff
    d = dist(c,n)
    F = -1.*dLJ(d)  
    a = float( c[0] - n[0] )
    b = float( c[1] - n[1] )
    x = F * a / d
    y = F * b / d
    out = [sum(i) for i in zip(out, [x,y])]
  return out

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

# calculate stress tensor the formal way (average over all atoms)
def stress(l):  # l = list of coordinates
  stress_xx = 0
  stress_xy = 0
  stress_yx = 0
  stress_yy = 0
  for a in l:
    for n in NN(a, RCUT):
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

# plot the initial grid of atoms
plt.scatter([i[0] for i in coords], [i[1] for i in coords])
plt.xlabel('x coordinate [Ang]')
plt.ylabel('y coordinate [Ang]')
plt.show()

# plot the initial radial distribution function for the central atom using a histogram
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

step = 1
tot_energy = []
max_force = []
s_xx = []
s_xy = []
s_yx = []
s_yy = []
while step <= MAX_STEP:
  if step % int(MAX_STEP/10) == 0:
    max_f = 0
    eng_t = 0
    for c in coords:
      fo = force(c)
      max_f = max(max_f, (fo[0]**2 + fo[1]**2)**0.5)
      for n in NN(c, RCUT):
        eng_t += LJ(dist(c,n))      
    eng_t /= 2.  # correct double-counting
    tot_energy.append(eng_t)
    max_force.append(max_f)
    s_all = stress(coords)
    s_xx.append(s_all[0])
    s_xy.append(s_all[1])
    s_yx.append(s_all[2])
    s_yy.append(s_all[3])

  for ind, c in enumerate(coords):
    f = force(c)
    coords[ind] = (c[0]+LAMBDA*f[0], c[1]+LAMBDA*f[1])

  step += 1

# plot the final grid of atoms
plt.scatter([i[0] for i in coords], [i[1] for i in coords])
plt.xlabel('x coordinate [Ang]')
plt.ylabel('y coordinate [Ang]')
plt.show()

# plot the final radial distribution function for the central atom using a histogram
dists = [dist(center, i) for i in [j for j in coords if j!=center]]
plt.hist(dists, range=(0,RDF_MAX), rwidth=0.5, bins=20)
plt.xlabel('Distance [Ang]')
plt.ylabel('Radial distribution function [# atoms]')
plt.show()

# plot the final rdf the formal way
plt.plot([i*RDF_MAX/100. for i in range(1,101)], rdf(coords))
plt.xlabel('Distance [Ang]')
plt.ylabel('Radial distribution function [# atoms]')
plt.show()

# plot the energy and force as a function of step
x_axis = [i for i in range(int(MAX_STEP/10), int(MAX_STEP+1), int(MAX_STEP/10))]
fig,ax = plt.subplots()
ax.plot(x_axis, tot_energy, color="red", marker="o")
ax.set_xlabel("Iteration")
ax.set_ylabel("Total Energy [eV]",color="red")
ax2=ax.twinx()
ax2.plot(x_axis, max_force, color="blue",marker="o")
ax2.set_ylabel("Max Force [eV / Ang]",color="blue")
plt.show()

# plot the stress tensor elements as a function of step
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(x_axis, s_xx)
axs[0, 0].set_title('XX')
axs[0, 1].plot(x_axis, s_yy, 'tab:orange')
axs[0, 1].set_title('YY')
axs[1, 0].plot(x_axis, s_xy, 'tab:green')
axs[1, 0].set_title('XY')
axs[1, 1].plot(x_axis, s_yx, 'tab:red')
axs[1, 1].set_title('YX')
for ax in axs.flat:
    ax.set(xlabel='Iteration', ylabel='Stress')
# Hide x labels and tick labels for top plots and y ticks for right plots.
for ax in axs.flat:
    ax.label_outer()
plt.show()
