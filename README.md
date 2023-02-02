# Astro 401 - Homework #1

**Questions:** Consider two bodies of mass m1 and m2 separated by a distance d traveling in circular orbits around their mutual center of mass. Take the angular momentum vector of the system to be pointing in the +z direction. In a frame that co-rotates with the orbital motion there are five locations (called the Lagrange points) where the effective acceleration felt by a test particle vanishes. The acceleration arises from the gravity due to m1 and m2 plus the rotation of the system as a whole. If the masses are at points x1 < x2 on the x axis (sod=x2−x1),thenbyconventionL3 <x1 <L1 <x2 <L2. TheL4 andL5 pointslie off the x axis and form equilateral triangles with the points x1 and x2. Conventionally, L4 is taken to be in the +y direction, L5 in the −y direction.

(1). (30 points) Write a program to explore this system by computing the effective gravity (a vector) and potential (a scalar) at specific grid points in the xy-plane. Use units such that the gravitational constant '''G ≡ 1'''. The program should take as input the mass of both bodies and their separation. You may “hardwire” the grid dimensions into your code if you wish. The output should be the potential and x and y components of the acceleration at each grid point.

''' python

########### Libraries ############
import numpy as np
import matplotlib.pyplot as plt
import math

########### Loading Files ############
d = np.loadtxt('lagrange_mass.dat', unpack=1)
r = np.loadtxt('lagrange_grid.dat', unpack=1)
s = np.loadtxt('lagrange_root.dat', unpack=1)
n = int((len(r[0]))**0.5)
x = r[0].reshape((n,n))
y = r[1].reshape((n,n))
p = r[2].reshape((n,n))
gx = r[3].reshape((n,n))
gy = r[4].reshape((n,n))

########### Plot Size ############
fig,ax = plt.subplots(figsize=(8,8))
ax.set_aspect('equal', 'box')

########### Contour Lines ############
cs = ax.contour(x, y, np.log10(-p), 40)

########### Labels ############
ax.set_xlabel('X', size='x-large')
ax.set_ylabel('Y', size='x-large')

########### Dots (Colors/Labels) ############
    #Lagrange Point 
ax.plot(s[0], s[1], marker='o', linestyle='none', markersize=7, markerfacecolor='r', markeredgecolor='k')
ax.plot(s[0][0], s[1][0], marker='o', linestyle='none', markersize=7, markerfacecolor='blue', markeredgecolor='k')
ax.plot(s[0][1], s[1][1], marker='o', linestyle='none', markersize=7, markerfacecolor='lightgreen', markeredgecolor='k')
ax.plot(s[0][2], s[1][2], marker='o', linestyle='none', markersize=7, markerfacecolor='orange', markeredgecolor='k')
ax.plot(s[0][3], s[1][3], marker='o', linestyle='none', markersize=7, markerfacecolor='magenta', markeredgecolor='k')
    
for i in range(len(s[0])):
    ax.text(s[0][i], s[1][i]+ 0.05, 'L%d'%(i+1), color='r', size='x-large', horizontalalignment='center', verticalalignment='bottom')
    
    # m1 = Sun
ax.text(d[2], 0.0, 'M1', color='k', size='large', horizontalalignment='center', verticalalignment='bottom')
ax.plot(d[2], 0.0, marker='o', linestyle='none', markersize=20, markerfacecolor='yellow', markeredgecolor='yellow')
    
    # m2 = Earth
ax.text(d[3], 0.0, 'M2', color='k', size='large', horizontalalignment='center', verticalalignment='bottom')
ax.plot(d[3], 0.0, marker='o', linestyle='none', markersize=10, markerfacecolor='teal', markeredgecolor='teal')

########### Calculations ############
g = (gx**2 + gy**2)**0.5
gx *= g**0.25/g
gy *= g**0.25/g

########### Gravitiy Arrows ############
ax.quiver(x, y, gx, gy, scale=150 ) # used 150 for case 1 and 190 for case 2

########### Save and Plot ############
plt.savefig('case1.jpeg')
plt.show()


'''
