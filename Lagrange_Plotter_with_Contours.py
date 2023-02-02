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
ax.set_title('Case 1: m1 = 3, m2 = 1, d = 1', size='x-large')
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
