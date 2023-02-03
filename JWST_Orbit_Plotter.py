
#on Jupyter, run plg.py
#on xterm, python plg.py

import numpy as np

import matplotlib.pyplot as plt
import math
import turtle


d = np.loadtxt('lagrange_mass.dat', unpack=1)
r = np.loadtxt('lagrange_grid.dat', unpack=1)
s = np.loadtxt('lagrange_root.dat', unpack=1)
n = int((len(r[0]))**0.5)
x = r[0].reshape((n,n))
y = r[1].reshape((n,n))
p = r[2].reshape((n,n))
gx = r[3].reshape((n,n))
gy = r[4].reshape((n,n))

fig,ax = plt.subplots(figsize=(10,10))
ax.set_aspect('equal', 'box')
cs = ax.contour(x, y, np.log10(-p), 1)
ax.set_xlabel('X', size='x-large')
ax.set_ylabel('Y', size='x-large')
ax.set_title('JWSTs Orbit', size = 'xx-large')
## dots
ax.plot(s[0], s[1], marker='o', linestyle='none', markersize=10,markerfacecolor='r', markeredgecolor='r')
ax.plot(s[0][0], s[1][0], marker='o', linestyle='none', markersize=10, markerfacecolor='blue', markeredgecolor='blue')
ax.plot(s[0][1], s[1][1], marker='o', linestyle='none', markersize=10, markerfacecolor='lightgreen', markeredgecolor='lightgreen')
ax.plot(s[0][2], s[1][2], marker='o', linestyle='none', markersize=10, markerfacecolor='orange', markeredgecolor='orange')
ax.plot(s[0][3], s[1][3], marker='o', linestyle='none', markersize=10, markerfacecolor='magenta', markeredgecolor='magenta')
    
## m1
ax.text(d[2], 0.0, 'M1', color='k', size='large',horizontalalignment='center')
ax.plot(d[2], 0.0, marker='o', linestyle='none', markersize=50, markerfacecolor='yellow',markeredgecolor='yellow')
## m2
ax.text(d[3], 0.0, 'M2', color='k', size='large',horizontalalignment='center')
ax.plot(d[3], 0.0, marker='o', linestyle='none', markersize=10, markerfacecolor='teal',markeredgecolor='teal')


for i in range(len(s[0])):
    ax.text(s[0][i], s[1][i]+ 0.05, 'L%d'%(i+1), color='r', size='x-large',
            horizontalalignment='right', verticalalignment='bottom')
    
## https://dev.to/baronsindo/plotting-a-circle-in-python-with-angles-m0c
# used link above to figure out how to draw a circle

pi=3.14159265

angles = np.linspace(0*pi, 2*pi,150)
print(angles)
xs = np.cos(angles)*1.15
ys = np.sin(angles)*1.15
xr = np.cos(angles)
yr = np.sin(angles)

plt.plot(xs,ys,color = 'lightgreen')
plt.plot(xr,yr,color = 'teal')

plt.gca().set_aspect('equal')


#circle1 = plt.Circle((d[2], 0.0), 1.5, color='r')
#ax.add_patch(circle1)


#plt.Circle((s[0],0), radius=1.01)



g = (gx**2 + gy**2)**0.5
gx *= g**0.25/g
gy *= g**0.25/g

#ax.quiver(x, y, gx, gy, scale=50)
plt.savefig('lg4.jpeg')
plt.show()

##planet = circle(x,y,8)
#Implementing ellipse equations to generate the values needed to plot an ellipse
#Using only the planet's min (m) and max (M) distances from the sun
#Equations return '2a' (the ellipses width) and '2b' (the ellipses height)
"""
def PlanetOrbit(Name, M, m):
    w, h = OrbitLength(M, m)
    Xoffset= ((M+m)/2)-m
    Name = Ellipse(xy=((Xoffset),0), width=w, height=h, angle=0, linewidth=1, fill=False)
    ax.add_artist(Name)
    
    
def OrbitLength(M, m):
    a=(M+m)/2
    c=a-m
    e=c/a
    b=a*(1-e**2)**0.5
    print(a)
    print(b)
    return 2*a, 2*b

PlanetOrbit('Earth', 69.8, 46.0)
#This function uses the returned 2a and 2b for the ellipse function's variables


#These are the arguments taken from hyperphysics.phy-astr.gsu.edu/hbase/solar/soldata2.html
#They are the planet names, max and min distances, and their longitudinal angle
#Also included is Halley's Comet, used to show different scale  and eccentricity






#Also generating the orbit offset (putting the sun at a focal point) using M and m
def PlanetOrbit(Name, M, m):
    w, h = OrbitLength(M, m)
    Xoffset= ((M+m)/2)-m
    Name = Ellipse(xy=((Xoffset),0), width=w, height=h, angle=0, linewidth=1, fill=False)
    ax.add_artist(Name)




"""