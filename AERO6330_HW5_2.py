# Zane Grothe
# AERO 6330
# HW 5
# 3/29/22

# CONVERTED TO PYTHON 4/9/24

import sympy as sp
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Problem 2 ~~~~~~~~~~~~~~~~~~~~

# GIVENS

M1 = 5.9722e24  # Mass of Earth
M2 = 7.342e22  # Mass of Moon
M = M1 + M2  # Total mass
tol = 1e-8  # Tolerance

# Nondimensional mass parameters
mu = M2/M

# FIND L1
x = sp.symbols('x')
# For all solutions (5 of them):
# dUdx = x - (1-mu)*(x+mu)/sqrt((x+mu)**2)**3 - mu*(x+mu-1)/sqrt((x+mu-1)**2)**3
# For position of L1, x must be >-mu and <1-mu so:
dUdx = x - (1 - mu)/(x + mu)**2 + mu/(x + mu - 1)**2
d2Udx2 = sp.diff(dUdx, x)

fraction = 1
x0 = 0.5  # initialize loop with position

# Newton's Method
while abs(fraction) > tol:
    fraction = sp.N(dUdx.subs(x, x0))/sp.N(d2Udx2.subs(x, x0))
    x0 = x0 - fraction

L1x = float(x0)  # L1x = scalar
L1y = 0

# FIND Trajectory IC's
x, y = sp.symbols('x y')

d = sp.sqrt((x + mu)**2 + y**2)  # Distance from Earth to Spacecraft
r = sp.sqrt((x + mu - 1)**2 + y**2)  # Distance from Moon to Spacecraft

U = (x**2 + y**2)/2 + (1 - mu)/d + mu/r  # Pseudo-potential

Uxx = sp.diff(U, x, 2)  # Second derivatives
Uyy = sp.diff(U, y, 2)

Uxx = sp.N(Uxx.subs([(x, L1x), (y, 0)]))  # Solve second derivatives
Uyy = sp.N(Uyy.subs([(x, L1x), (y, 0)]))

# Constants
B1 = float(2 - (Uxx + Uyy)/2)
B2 = float(-Uxx*Uyy)
s = np.sqrt(B1 + np.sqrt(B1**2 + B2))
B3 = float((s**2 + Uxx)/(2*s))

# L1 offset
delta = 0.01

# Initial Conditions
ep0 = delta
et0 = 0
ep0_dot = et0*s/B3
et0_dot = -ep0*s*B3
# linear_solution_initial_velocity = et0_dot
xy0 = [L1x + ep0, et0, ep0_dot, et0_dot]


# Numerically integrate
def CR3BP_EOM(xy0, t):
    # Distances
    d = np.sqrt((mu + xy0[0])**2 + (xy0[1])**2)
    r = np.sqrt((1 - mu - xy0[0])**2 + (xy0[1])**2)

    # Outputs
    return [
        xy0[2],
        xy0[3],
        xy0[0] + 2*xy0[3] + (1 - mu)*(-mu - xy0[0])/(d**3) + mu*(1 - mu - xy0[0])/(r**3),
        xy0[1] - 2*xy0[2] - (1 - mu)*(xy0[1])/(d**3) -mu*xy0[1]/(r**3),
    ]


# Period
P = 2*np.pi/s

tspan = np.linspace(0, P, 101)
xy = odeint(CR3BP_EOM, xy0, tspan)


# PLOTS

plt.figure(1) # (follow trajectory no matter how big)
# Planets
plt.plot(-mu, 0, 'o:g', ms=12) # Green Earth
plt.plot(1 - mu, 0, 'o:m', ms=6) # Magenta Moon
# L1
plt.plot(L1x, 0, 's:k', ms=4) # L1 point on black box
# Trajectory
plt.plot(xy[:,0], xy[:,1], 'r')
plt.title(f'Full View Nonlinear Dynamics for {delta} from L1')
plt.xlabel('X')
plt.ylabel('Y')
# plt.set_aspect('equal', 'box')

plt.figure(2) # (now zoom in on L1)
# L1
plt.plot(L1x, 0, 's:k', ms=8) # L1 point on black box
# Trajectory
plt.plot(xy[:,0], xy[:,1], 'r')
plt.xlim([L1x*(1 - 5*abs(delta)/L1x), L1x*(1 + 5*abs(delta)/L1x)])
plt.ylim([et0 - L1x*5*(abs(delta)/L1x), et0 + L1x*5*(abs(delta)/L1x)])
plt.title(f'Scaled View Nonlinear Dynamics for {delta} from L1')
plt.xlabel('X')
plt.ylabel('Y')
# plt.set_aspect('equal', 'box')

plt.figure(3) # Are x and y periodic?
plt.plot(xy[:,0],'r')
plt.plot(xy[:,1],'b')
plt.title(f'Range of Position Movements for {delta} from L1')
plt.xlabel('Nondim. Time')
plt.ylabel('Position Change')
plt.legend(['x', 'y'])

plt.show()

