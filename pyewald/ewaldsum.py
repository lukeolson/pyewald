from numpy import dot, sqrt
from numpy.linalg import norm
from scipy.special import erfc
from scipy import exp, pi


def energy():
    pass


def potential(i, r, q, cell, area, invcell, alpha=5, nmax=3, mmax=9):
    """
    i : potential location
    r : list of radii
    q : list of charges
    alpha : Ewald parameter
    nmax : real space box cutoff
    mmax : Fourier space box cutoff
    """
    Vr = potential_realsum(i, r, q, cell, alpha, nmax)
    Vf = potential_recipsum(i, r, q, invcell, alpha, mmax, area)
    Vs = potential_selfsum(i, r, q, cell, alpha)
    return Vr+Vf+Vs


def potential_realsum(i, r, q, cell, alpha, nmax):
    Vr = 0
    for j in range(0, len(q)):
        rij = r[i, :] - r[j, :]
        Vrloc = 0
        for n2 in range(-nmax, nmax+1):
            for n1 in range(-nmax, nmax+1):
                for n0 in range(-nmax, nmax+1):
                    if max([abs(n0), abs(n1), abs(n2)]):
                        n = n0*cell[0, :]+n1*cell[1, :]+n2*cell[2, :]
                        rn = norm(rij - n)
                        Vrloc += erfc(alpha*rn)/rn
        Vr += q[j]*Vrloc
    return Vr


def potential_recipsum(i, r, q, cell, alpha, mmax, area):
    Vf = 0
    for j in range(0, len(q)):
        rij = r[i, :] - r[j, :]
        Vfloc = 0
        for n2 in range(-mmax, mmax+1):
            for n1 in range(-mmax, mmax+1):
                for n0 in range(-mmax, mmax+1):
                    if max([abs(n0), abs(n1), abs(n2)]):
                        m = 2*pi*(n0*cell[0, :]+n1*cell[1, :]+n2*cell[2, :])
                        Vfloc += exp(1.j * dot(m, rij) - dot(m, m) /
                                     (alpha**2 * 4))/(dot(m, m)/(4*pi**2))

        Vf += q[j]/(pi*area) * Vfloc.real
    return Vf


def potential_selfsum(i, r, q, cell, alpha):
    Vs = 0
    for j in range(0, len(q)):
        if i == j:
            Vs -= 2*q[j]*alpha/sqrt(pi)
        else:
            rn = norm(r[i, :] - r[j, :])
            Vs += q[j]*erfc(alpha*rn)/rn
    return Vs
