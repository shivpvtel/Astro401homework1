# Astro 410 - Homework #1

****Question:**** <br />
Consider two bodies of mass `m1` and `m2` separated by a distance `d` traveling in circular orbits around their mutual center of mass. Take the angular momentum vector of the system to be pointing in the `+z` direction. In a frame that co-rotates with the orbital motion there are five locations (called the Lagrange points) where the effective acceleration felt by a test particle vanishes. The acceleration arises from the gravity due to `m1` and `m2` plus the rotation of the system as a whole. If the masses are at points `x1 < x2` on the x axis (so `d = x2 −  x1`),then by convention `L3 < x1 < L1 < x2 < L2`. The `L4` and `L5` points lie off the x axis and form equilateral triangles with the points `x1` and `x2`. Conventionally, `L4` is taken to be in the `+y` direction, `L5` in the `−y` direction.

1. (30 points) Write a program to explore this system by computing the effective gravity (a vector) and potential (a scalar) at specific grid points in the xy-plane. Use units such that the gravitational constant `G ≡ 1`. The program should take as input the mass of both bodies and their separation. You may “hardwire” the grid dimensions into your code if you wish. The output should be the potential and x and y components of the acceleration at each grid point.

So I first started my code in C++ but I could not get the gnuplot to work. My code was also not outputting the correct lagrangian points and format for Gnuplot (or the python code) to understand. During tuesday class Dr.Gulis' said that I need switch to python because she does not know how to use gnuplot. However, taking only CMPSC 121/122, I had no choice but to learn python within a matter of hours to figure out how to use it for my code. I used Dr.Gulis' Python code `two-body-potential-and-lagrangianpoints.py` now renamed to `Lagrange_Points_Data_Calc.py` to get the correct outputs as well as the `plot-potential-and-forcevector.py` which has been renamed to `Lagrange_Plotter_with_Contours.py `. For the the JWST's Orbit i edited `Lagrange_Plotter_with_Contours.py` to better graph the orbit for the earth and JWST (L2).

I versions I used in Jupyter notebook is:
- IPython          : 7.31.1
- ipykernel        : 6.15.2
- ipywidgets       : 7.6.5
- jupyter_client   : 7.3.4
- jupyter_core     : 4.11.1
- jupyter_server   : 1.18.1
- jupyterlab       : 3.4.4
- nbclient         : 0.5.13
- nbconvert        : 6.4.4
- nbformat         : 5.5.0
- notebook         : 6.4.12
- qtconsole        : 5.3.2
- traitlets        : 5.1.1

Dr.Gulis' code was writen using Python but im not sure if it was Python3 (I dont know anything in python). My code will have comments that explain the equations.

``` python
"""
usage: python lps.py m1 m2 d
m1: mass of body 1
m2: mass of body 2
d:  distance between two bodies
"""
########### Libraries ############

import matplotlib.pyplot as plt
from numpy import *
import sys

########### Files ############

GRID_FILE = 'lagrange_grid.dat'
ROOT_FILE = 'lagrange_root.dat'
MASS_FILE = 'lagrange_mass.dat'

########### Grids ############

XMIN = -2.0
XMAX = 2.0
DX = 0.08

YMIN = -2.0
YMAX = 2.0
DY = 0.08

FEPS = 1e-6

########### User Input ############

m1 = float(input("What is the mass of m1: ")) # changed from the original code
m2 = float(input("What is the mass of m2: ")) # changed from the original code
d = float(input("Distance from m1 to m2: ")) # changed from the original code

mSun = 1 # in solar units
mEarth = 1/333000 # in solar units 
mTotal = mSun + mEarth

########### Calculations ############

m = m1 + m2  ## This is the total Mass which is dominated by the suns mass
w2 = m / (d ** 3) ## this is the angular momentum 
x1 = -m2 * d / m  ## center of mass for x1
x2 = m1 * d / m   ## center of mass for x2
xc = 0.5 * (x1 + x2) #3 midpoint between x1 and x2

########### Root finding algorithm ############

"""
def RootFindingAlgo(f, xOne, xTwo, tol): 
    if np.sign(f(xOne)) == np.sign(f(xTwo)):
        raise Exception( "root cannot be found")

    root = (xOne + xTwo)/2

    if np.sign(f(xOne)) == np.sign(f(root)): 
        return RootFindingAlgo(f, root, xTwo, tol)
    elif np.sign(f(xTwo)) == np.sign(f(root)):
        return RootFindingAlgo(f, xOne, root, tol)
    elif np.abs(f(root)) < tol:
        return root
""" ### This is a rootfinding algorithm that I tried to make but i did not implement it

def rtbis(f, xlo, xhi, eps):
    flo = f(xlo)    
    iter = 0
    while (xhi - xlo > 0.5 * eps * abs(xlo + xhi)):
        iter = iter + 1
        xm = 0.5 * (xlo + xhi)
        fm = f(xm)
        if abs(fm) < 1e-20:
            break
        if (flo < 0):
            if fm < 0:
                xlo = xm
            else:
                xhi = xm
        else:
            if fm < 0:
                xhi = xm
            else:
                xlo = xm
        if iter > 2000:

            print('maxiter reached in rtbis: %g %g %g %g %g %g' % (
            xlo, xhi, xm, xhi - xlo, 0.5 * eps * abs(xhi + xlo), fm))
            break
    return xm

########### Grav x function ############
def gravx(ir13, ir23, x, y):
    return -m1 * (x - x1) * ir13 - m2 * (x - x2) * ir23 + w2 * x
    
    ## This equation is the gravitational acceleration in the x direction

########### Grav y function ############
def gravy(ir13, ir23, x, y):
    return -m1 * y * ir13 - m2 * y * ir23 + w2 * y
  ## This equation is the gravitational acceleration in the x direction
  
  
def fx(x):
    r12 = (x - x1) ** 2 
    r22 = (x - x2) ** 2
    return gravx(1 / r12 ** 1.5, 1 / r22 ** 1.5, x, 0.0)

def fy(y):
    r12 = (xc - x1) ** 2 + y ** 2
    r22 = (xc - x2) ** 2 + y ** 2
    return gravy(1 / r12 ** 1.5, 1 / r22 ** 1.5, xc, y)

########### Writing Data to File ############

    """NEW CODE"""  ## this set of code was part of the new python file Dr.Gulis sent out
print('writting mass data to %s'%MASS_FILE)
with open(MASS_FILE, 'w') as f:
    f.write('%12.5E\n'%m1)
    f.write('%12.5E\n'%m2)
    f.write('%12.5E\n'%x1)
    f.write('%12.5E\n'%x2)
    """NEW CODE ENDS"""

print('Writing grid data to %s' % GRID_FILE)  
## This function iterates through the plot and its components and writes the output of the calculations to the .dat file.
xg = arange(XMIN, XMAX + 0.1 * DX, DX)
yg = arange(YMIN, YMAX + 0.1 * DY, DY)
with open(GRID_FILE, 'w') as f:
    for x in xg:
        for y in yg:
            y2 = y ** 2
            r2 = x * x + y * y
            r12 = (x - x1) ** 2 + y2
            r22 = (x - x2) ** 2 + y2
            r12 = max(r12, FEPS)
            r22 = max(r22, FEPS)
            isqrtr12 = 1 / sqrt(r12)
            isqrtr22 = 1 / sqrt(r22)
            ir13 = isqrtr12 / r12
            ir23 = isqrtr22 / r22
            gx = gravx(ir13, ir23, x, y)
            gy = gravy(ir13, ir23, x, y)
            p = - m1 * isqrtr12 - m2 * isqrtr22 - 0.5 * w2 * r2  ## This is the Gravitational accelaration equation
            f.write('%15.8E %15.8E %15.8E %15.8E %15.8E\n' % (x, y, p, gx, gy))


print('writing lagrange points to %s' % ROOT_FILE)
## This outputs the calculated lagrangian points using the root finding algorithm
with open(ROOT_FILE, 'w') as f:
    # L1: on x axis, between x1 and x2
        print('L1')
        lp = rtbis(fx, x1 + FEPS, x2 - FEPS, FEPS)
        f.write('%15.8E %15.8E\n' % (lp, 0.))
    # L2: on x-axis, to the right of x2
        print('L2')
        lp = rtbis(fx, x2 + FEPS, XMAX, FEPS)
        f.write('%15.8E %15.8E\n' % (lp, 0.))
    # L3: on x-axis, to the left of x1
        print('L3')
        lp = rtbis(fx, XMIN, x1 - FEPS, FEPS)
        f.write('%15.8E %15.8E\n' % (lp, 0.))
    # L4: on positive x=xc line
        print('L4')
        lp = rtbis(fy, FEPS, YMAX, FEPS)
        f.write('%15.8E %15.8E\n' % (xc, lp))
    # L5: on negative x=xc line
        print('L5')
        lp = rtbis(fy, YMIN, -FEPS, FEPS)
        f.write('%15.8E %15.8E\n' % (xc, lp))

```

2. (20 points) Use your favorite plotting program to plot vectors (for the effective acceleration) and contours (for the effective potential) for the cases where m1 = 3, m2 = 1, d = 1 and m1 = 100 , m2 = 1, d = 1.

- For this question im going to be using Matplotlib. Using the code below, I plotted both cases: <br />
        1. Case 1: m1 = 3,   m2 = 1, d = 1 <br />
        2. Case 2: m1 = 100, m2 = 1, d = 1
 
``` python
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
ax.plot(d[3], 0.0, marker='o', linestyle='none', markersize=10, markerfacecolor='teal', markeredgecolor='k')

########### Calculations ############
g = (gx**2 + gy**2)**0.5
gx *= g**0.25/g
gy *= g**0.25/g

########### Gravitiy Arrows ############
ax.quiver(x, y, gx, gy, scale=190)

########### Save and Plot ############
plt.savefig('case3.jpeg')
plt.show()
```

**Case 1:** This case shows the lagrangian points when the sun's mass is set to '3' and the earth is set to '1' giving us 5 Lagrange points. 3 of which (L1,L2,L3) lie on the x-axis and L4 lies in quadrant 1 and L5 lies in quadrant 4.


 <img src="https://github.com/shivpvtel/Astro401homework1/blob/main/case%201/case1.jpeg" width="400" height="400"/>

**Case 2:** This case shows the lagrangian points when the sun's mass is set to '100' and the earth is set to '1' giving us 5 Lagrange points as well.

 <img src="https://github.com/shivpvtel/Astro401homework1/blob/main/case%202/case2.jpeg" width="400" height="400"/>
  <br />
   <br />
 
 
 3. (20 points) Use a simple root solver to determine the locations of the Lagrange points for both these cases, using your a priori knowledge of their approximate locations (i.e., first do a search along the x axis for the L1, L2, and L3 points, then along the `x = (x1 + x2)/2` axis for the L4 and L5 points. Watch out for singularities!...).
 Using the rtbis function in the python code, we are able to calculate the lagrangian points.

For case 1 the Lagrangian points are:   <br />
     L1:  3.60743540E-01  0.00000000E+00  <br />
     L2:  1.26585852E+00  0.00000000E+00  <br />
     L3: -1.10316674E+00  0.00000000E+00  <br />
     L4:  2.50000000E-01  8.66025061E-01 <br />
     L5:  2.50000000E-01 -8.66025061E-01 <br />
    <br />
 
 For case 2 the Lagrangian points are:   <br />
     L1:  8.48624046E-01  0.00000000E+00 <br />
     L2:  1.14632026E+00  0.00000000E+00 <br />
     L3: -1.00412446E+00  0.00000000E+00 <br />
     L4:  4.90099010E-01  8.66025061E-01 <br />
     L5:  4.90099010E-01 -8.66025061E-01 <br />
    <br />  
    <br />
    <br />

4. (30 points) The James Webb Space Telescope was successfully launched on 12/24/2021. It is the largest and the most complex telescope ever launched into space. It will observe primarily the infrared light from faint and very distant objects. But all objects, including telescopes, also emit infrared light. To avoid swamping the very faint astronomical sig- nals with radiation from the telescope, the telescope and its instruments must be very cold. Therefore, JWST has a large shield that blocks the light from the Sun, Earth, and Moon. To have this work, JWST must be in an orbit where all three of these objects are in about the same direction. Please calculate and plot the ideal location and orbit for JWST.
    
To figure out what was the ideal location for JWST i used the same code used to plot in #2. In which L2 was the only spot that could JWST could be in to its solar shield to project itself from the sun, as well as the far enough to avoid the moons glare. to determine and plot this i used the `Lagrange_Points_Data_Calc.py` to calculate the Grid data as well as the Lagrange points, however for the inputs i converted to solar units so `m1 = 1` `m2 = 1/333000` because in the earth is 0.000003 times the suns mass. and `d=1` becasue it represented 1 AU. The code then outputted three files, `lagrange_grid.dat` which contained the x,y points, the gravitational accelaration and the gravatational potential, `lagrange_root.dat` which contained the lagrangian points in x,y coordinates, `lagrange_mass.dat` which contained the masses and out calculate inputs. I then used the `Lagrange_Plotter_with_Contours.py` but changed it up to better fit my plot and named it `JWST_Orbit_Plotter.py`. I added circles to represent the orbits of the earth (depicted in the teal circle) and I added another circle in light green that depicted the L2 orbit. Again I had no idea how to use python so I had to google how to draw a circle and the reference link I used is commented out in the code it self.  Below you will find the plot depicting the orbit of JWST from a top down view of our solar system.

The Lagrangian points for this plot was:
- L1:  8.48624046E-01  0.00000000E+00
- L2:  1.14632026E+00  0.00000000E+00
- L3: -1.00412446E+00  0.00000000E+00
- L4:  4.90099010E-01  8.66025061E-01
- L5:  4.90099010E-01 -8.66025061E-01


 <img src="https://github.com/shivpvtel/Astro401homework1/blob/main/jwst/JWST's%20Orbit.jpeg" width="400" height="400"/>


I am also going to attach my C++ code (`main.cpp`) because i worked hard on it, and i belive my root finding algorithm was spot on. I worked alongside mohammed and he taught me new functions like the pair function and a couple other things that were neccassary for this Homework assingment. 
