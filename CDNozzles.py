
"""Module containing converging-diverging nozzle equations e.g. property ratios to sonic velocity values, Mach number from area ratio"""

import numpy as np

PRECISION = 1e-7

def SonicTemperatureRatio(M, gamma=1.4):
    '''Returns T*/T from from incident Mach number'''
    top = 1+(gamma-1)/2*M**2
    bottom = (gamma+1)/2
    TstarT = (top/bottom)
    return TstarT

def SonicPressureRatio(M, gamma=1.4):
    '''Returns P*/P from incident Mach number'''
    top = 1+(gamma-1)/2*M**2
    bottom = (gamma+1)/2
    PstarP = (top/bottom)**(gamma/(gamma-1))
    return PstarP

def SonicAreaRatio(M, gamma=1.4):
    '''Returns A*/A from incident Mach number'''
    top = (gamma+1)/2
    bottom = 1+(gamma-1)/2*M**2
    AstarA = M*(top/bottom)**((gamma+1)/(2*(gamma-1)))
    return AstarA

def MachfromAreaRatio(A, gamma=1.4, case=0):
    """Returns Mach number from area ration (A/A*)
    Supersonic solution: case = 0
    Subsonic solution: case = 1"""

    P = 2/(gamma+1)
    Q = 1 - P
    
    if case == 0: # Supersonic case
        R = A**(2*Q/P)
        a = Q**(1/P)

    elif case == 1: # Subsonic case
        R = A**2
        a = P**(1/Q)

    r = (R-1)/(2*a)
    X_new = 1/(1+r+np.sqrt(r*r+2*r))
    X = 0
    
    while abs(X - X_new) > PRECISION: #Newton-raphson
        X = X_new

        if case == 0: # Supersonic case
            F = (P*X+Q)**(1/P) - R*X
            f = (P*X+Q)**((1/P)-1) - R

        elif case == 1: # Subsonic case
            F = (P+Q*X)**(1/Q) - R*X
            f = (P+Q*X)**((1/Q)-1) - R
    
        X_new = X - F/f

    if case == 0:
        M = 1/np.sqrt(X)
    elif case == 1: # Subsonic case
        M = np.sqrt(X)
    
    return M

def find_Astar(P, P0, R, M, A, gamma=1.4):
    a = (1+(gamma-1)/2*M**2)**.5
    b = M/P0
    c = (gamma+1)/(gamma-1)
    d = np.sqrt(gamma/R* (2/(gamma+1))**c)**-1

    Astar = P/R*A*np.sqrt(gamma*R)*a*b*d

    return Astar