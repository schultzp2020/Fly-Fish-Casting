# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as la
from scipy.integrate import odeint

data = np.genfromtxt('input_data.csv', delimiter=',')
data_ith = np.genfromtxt('input_data_ith.csv', delimiter=',')

# get the phi array
phi = data_ith[0]

# get the mass array
m = data_ith[1]

# get the length array
L = data_ith[2]

# get the radius array
r = data_ith[3]

# get the initial x position
initial_x_position = data[0]

# get the initial y position
initial_y_position = data[1]

# get the initial x velocity
initial_x_velocity = data[2]

# get the initial y velocity
initial_x_velocity = data[3]


def dYdt(t, Y):
    n = len(Y)
    phi = Y[0:n/2]  # we know n is even since Y = [phi phidot]
    phidot = Y[n/2:n]

    # build RHS vector Q
    Q = build_Q(phi, phidot)

    # matrix [I 0;0 ainv]  #see eqn 32 of wang paper
    M = np.zeros(shape(n, n))
    for i in range(n/2):
        M[i, i] = 1.

    # build little a, then invert and put into M
    # a needs to be n/2 by n/2 np.array/matrix
    a = build_a(phi, phidot)
    ainv = la.inv(a)

    for i in range(n/2):
        M[n/2+i, n/2+i] = ainv[i, i]

    b = np.zeros(shape(n, 1))
    for i in range(n/2):
        b[i] = phidot[i]
        b[n/2+i] = Q[i]
    Ydot = M@b

    return Ydot

# helper routines to build RHS vectors and matrices


def build_Q(phi, phidot, m, L, r, xt, yt, xtdot, ytdot, xtddot, ytddot):
    n = len(phi)
    Q = np.zeros(shape(n, 1))
    g = -9.8
    rhoa = 0.0
    CS = 0.005

    A = build_bigA(phi, phidot, m)

    for i in range(n):
        s = m[i]*0.5*L[i]*(g*np.cos(phi[i])-np.sin(phi[i])
                            * xtddot+np.cos(phi[i])*ytddot)
        for k in range(i, n):
            if k != i:
                s += m[k]*L[i]*(g*np.cos(phi[i])-np.sin(phi[i])
                                * xtddot+np.cos(phi[i])*ytddot)
        Q[i, 1] = -s

        s = 0
        for k in range(i, n):
            if k != i:
                s += A[i, k]*L[i]*L[k]*phidot[k]**2*np.sin(phi[k]-phi[i])
        Q[i, 1] += s

        # cos(phi[k] - phi[i]) = 1 when i = k therefore cos is not needed
        s = getCD(getRe(phi, phidot, i, r, L, xtdot, ytdot)) * \
            rhoa*r[i]*L[i]*getViD(xtdot, ytdot, phi[i]) * \
            np.abs(getViD(xtdot, ytdot, phi[i]))*0.5*L[i]
        for k in range(i, n):
            if k != i:
                s += getCD(getRe(phi, phidot, k, r, L, xtdot, ytdot)) * \
                    rhoa*r[k]*L[k]*getViD(xtdot, ytdot, phi[k]) * \
                    np.abs(getViD(xtdot, ytdot, phi[k])
                           )*L[i]*np.cos(phi[k]-phi[i])
        Q[i, 1] += -s

        # sin(phi[i] - phi[k]) = 0 when i = k therefore s = 0
        s = 0
        for k in range(i, n):
            if k != i:
                s += CS*rhoa*np.pi * \
                    r[k]*L[k]*getViS(xtdot, ytdot, phi[k]) * \
                    np.abs(getViS(xtdot, ytdot, phi[k])
                           )*L[i]*np.sin(phi[i]-phi[k])
        Q[i, 1] += -s

    return Q


def build_a(phi, phidot, A, L, m):
    # build little a matrix from wang paper
    # m is vector of segment masses
    # L is vector of segment lengths
    # A is big A matrix
    n = len(phi)
    a = np.zeros(shape(n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                s = m[k]*L[k]*L[k]*0.25
                for k in range(i, n):
                    s += m[k]*L[k]*L[k]
                a[i, j] = s
            else:
                a[i, j] = A[i, j]*L[i]*L[j]*np.cos(phi[i]-phi[j])
    return a


def build_bigA(phi, phidot, m):
    # build big A matrix from wang paper
    # m is vector of masses of segments
    n = len(phi)
    A = np.zeros(shape(n, n))
    for i in range(n):
        for j in range(n):
            if j > i:
                s = 0
                for k in range(j, n):
                    s += m[k] - 0.5*m[j]
                A[i, j] = s
            else:
                s = 0
                for k in range(i, n):
                    s += m[k] - 0.5*m[i]
                A[i, j] = s
    return A

# Reynolds number, drag, and dissipation calculations


def getRe(phi, phidot, i, r, L, xtdot, ytdot):
    # return the Reynolds number for segment i
    # r is the radius of the segment
    # L is a vector of segment lengths
    Re = 0.
    Re = 1.364e5*r*np.sqrt(getxdot(phi, phidot, i, xtdot, L)
                           ** 2 + getydot(phi, phidot, i, ytdot, L)**2)
    return Re


def getCD(Re):
    # return drag coefficient in normal direction
    # Re is Reynolds number
    CD = 0.
    if Re < 1:
        return 7.16
    if Re < 34:
        return 7.16*Re**-0.42
    if Re < 1580:
        return 3.02*Re**-0.165
    else:
        return 0.9
    return CD

# note CS = 0.005 constant dont need a function here

# routines to get positions and velocities


def getx(phi, i, xt, L):
    # get x at position i
    # rod tip x coordinate at xt
    # L is vector of segment Lengths
    s = xt
    for j in range(i-1):
        s += L[j]*np.cos(phi[j])
    s += L[i]*np.cos(phi[i])*0.5
    return s


def gety(phi, i, yt, L):
    # get y at position i
    # rod tip y coordinate at yt
    # L is vector of segment Lengths
    s = yt
    for j in range(i-1):
        s += L[j]*np.sin(phi[j])
    s += L[i]*np.sin(phi[i])*0.5
    return s


def getxdot(phi, phidot, i, xtdot, L):
    # get xdot at position i
    # rod tip x coordinate movign with velocity xtdot
    # L is vector of segment Lengths
    s = xtdot
    for j in range(i-1):
        s += -L[i]*np.sin(phi[j])*phidot[j]
    s += -L[i]*0.5*np.sin(phi[i])*phidot[i]
    return s


def getydot(phi, phidot, i, ytdot, L):
    # get ydot at position i
    # rod tip y coordinate moving with velocity ytdot
    # L is vector of segment Lengths
    s = ytdot
    for j in range(i-1):
        s += L[j]*np.cos(phi[j])*phidot[j]
    s += L[i]*0.5*np.cos(phi[i])*phidot[i]
    return s


def getViD(xtdot, ytdot, phi_i):
    # get VD at i
    return -xtdot*np.sin(phi_i)+ytdot*np.cos(phi_i)


def getViS(xtdot, ytdot, phi_i):
    # get VS at i
    return xtdot*np.cos(phi_i)+ytdot*np.sin(phi_i)
