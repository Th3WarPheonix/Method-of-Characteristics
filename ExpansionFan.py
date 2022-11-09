
"""Module containing expansion fan equations e.g. property ratios acorss a fan, mach from PM function value"""

import numpy as np

def press_ratio_ExpFan(M1, M2, gamma=1.4):
    """Returns P2/P1"""
    a = 1 + (gamma-1)/2*M2**2
    b = 1 + (gamma-1)/2*M1**2
    p1p2 = (a/b)**(gamma/(gamma-1))
    return p1p2**-1

def temp_ratio_ExpFan(M1, M2, gamma=1.4):
    """Returns T2/T1"""
    a = 1 + (gamma-1)/2*M2**2
    b = 1 + (gamma-1)/2*M1**2
    T1T2 = a/b
    return T1T2**-1

def PranMeyer(M1, gamma=1.4, deg=True):
    """Returns the PM function based on incident mach number"""
    A = np.sqrt((gamma+1)/(gamma-1))
    B = np.atan(np.sqrt((gamma-1)/(gamma+1)*(M1**2-1)))
    C = np.atan(np.sqrt(M1**2-1))
    # Search Prandtl-Meyer Function
    nu = (A*B-C)

    if deg:
        nu = nu*180/np.pi

    return nu

def inv_PranMeyer(v, deg=True, gamma=1.4):
    """Returns the Mach number at which the given value of the Prandtl-Meyer function occurs (given in degrees by default)"""
    a = np.sqrt((gamma+1)/(gamma-1))

    if deg: # Converting to radians for trig functions
        v *= np.pi/180

    if v > np.pi/2*(a-1): # max nu
        print('Given value of nu exceeds max value of nu')
        return None

    M = 1.3 
    M_next = (a**2-1) /((np.pi/2)*(a-1)-v) # Initial guess

    precision = 1e-10
    while abs(M-M_next)>precision:
        
        M = M_next
        A = a*np.atan(a**-1*np.sqrt(M**2-1)) # 1st part PM function
        B = np.atan(np.sqrt(M**2-1)) # 2nd part PM function
        C = (M*(gamma+1)/((gamma*M**2-M**2+2)*np.sqrt(M**2-1))) # 1st part PM function 

        D = 1/(M*np.sqrt(M**2-1)) # Derivative of 2nd part PM function

        M_next = M - (A-B-v)/(C-D) # Newton-Raphson formula
        
    return M

def solve_ExpFan(M1, theta, gamma=1.4):
    """Returns M2 T2/T1 and P2/P1 values as a dictionary from incident Mach number and turn angle

    Dictionary Keys
    ---------------
    'M2'\n
    'T2T1'\n
    'p2p1'"""
    nu = PranMeyer(M1, gamma)
    M2 = inv_PranMeyer(theta + nu, gamma=gamma)
    T2T1 = temp_ratio_ExpFan(M1, M2, gamma)
    p2p1 = press_ratio_ExpFan(M1, M2, gamma)

    results = {}
    results['M2'] = M2
    results['T2T1'] = T2T1
    results['p2p1'] = p2p1

    return results