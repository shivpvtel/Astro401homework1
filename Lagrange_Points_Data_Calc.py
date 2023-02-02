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

m = m1 + m2
w2 = m / (d ** 3)
x1 = -m2 * d / m
x2 = m1 * d / m
xc = 0.5 * (x1 + x2)

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
""" 

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

########### Grav y function ############
def gravy(ir13, ir23, x, y):
    return -m1 * y * ir13 - m2 * y * ir23 + w2 * y


def fx(x):
    r12 = (x - x1) ** 2
    r22 = (x - x2) ** 2
    return gravx(1 / r12 ** 1.5, 1 / r22 ** 1.5, x, 0.0)


def fy(y):
    r12 = (xc - x1) ** 2 + y ** 2
    r22 = (xc - x2) ** 2 + y ** 2
    return gravy(1 / r12 ** 1.5, 1 / r22 ** 1.5, xc, y)



########### Writing Data to File ############

    """NEW CODE"""
print('writting mass data to %s'%MASS_FILE)
with open(MASS_FILE, 'w') as f:
    f.write('%12.5E\n'%m1)
    f.write('%12.5E\n'%m2)
    f.write('%12.5E\n'%x1)
    f.write('%12.5E\n'%x2)
    """NEW CODE ENDS"""

print('Writing grid data to %s' % GRID_FILE)
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
            p = - m1 * isqrtr12 - m2 * isqrtr22 - 0.5 * w2 * r2
            f.write('%15.8E %15.8E %15.8E %15.8E %15.8E\n' % (x, y, p, gx, gy))


print('writing lagrange points to %s' % ROOT_FILE)

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
