# Austin Ryan

import numpy as np
import matplotlib.pyplot as plt
import ExpansionFan as EF

def main(M1, Mf, L, num_div, gamma=1.4):
    num_div += 1
    thetat = EF.PranMeyer(M1, gamma)-EF.PranMeyer(Mf, gamma)
    theta = np.linspace(0, thetat, num_div)

    kminus = EF.PranMeyer(M1, gamma) # deg
    nu = np.zeros(num_div)
    Mach = np.zeros(num_div)
    mu = np.zeros(num_div)
    slope = np.zeros(num_div)
    xs = np.zeros(num_div)
    ys = np.zeros(num_div)

    Mach[0] = M1
    nu[0] = kminus # deg
    mu[0] = np.arcsin(1/Mach[0])*180/np.pi # deg
    slope[0] = theta[0] + mu[0] # deg
    for i in range(1, num_div):
        nu[i] = nu[0] - theta[i] # deg
        Mach[i] = EF.inv_PranMeyer(nu[i], gamma=gamma) # deg
        mu[i] = np.arcsin(1/Mach[i])*180/np.pi # deg
        slope[i] = theta[i] + mu[i]
   
    slope *= np.pi/180
    theta *= np.pi/180
    
    for i in range(num_div-1):
        xf = L
        yf = L*np.tan(slope[0])

        a = xf - xs[i]
        b = yf - ys[i]
        top = b-a*np.tan(slope[i+1])
        bottom = np.tan(theta[i])-np.tan(slope[i+1])
        dx = top/bottom
        dy = dx*np.tan(theta[i])
        xs[i+1] = xs[i] + dx
        ys[i+1] = ys[i] + dy
        if i == 0: 
            plt.plot([xs[i], L], [ys[i], L*np.tan(slope[0])], color='red', label='Positive Characteristics')
            plt.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]], color='blue', marker='+', label='Inlet Contour')
        elif i == num_div-1-1:
            plt.plot([xs[i], L], [ys[i], L*np.tan(slope[0])], color='red')
            plt.plot([xs[i+1], L], [ys[i+1], L*np.tan(slope[0])], color='red')
            plt.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]], color='blue', marker='+')
        else:
            plt.plot([xs[i], L], [ys[i], L*np.tan(slope[0])], color='red')
            plt.plot([xs[i], xs[i+1]], [ys[i], ys[i+1]], color='blue', marker='+')

    length = -xs[0]+xs[-1]
    height = -ys[0]+ys[-1]
    plt.plot([xs[0], xs[-1]], [ys[0], ys[0]], color='gray', linestyle='--', label='Length ({})'.format(round(length, 1)))
    plt.plot([xs[-1], xs[-1]], [ys[0], ys[-1]], color='lightblue', linestyle='--', label='Height ({})'.format(round(height, 1)))

    plt.xlabel('X Coordinates')
    plt.ylabel('Y Coordinates')
    plt.title('Supersonic Inlet Contour and Characteristics\n\u03B3={}, M$_1$={}, M$_f$={}, \u03B8$_f$={}'.format(gamma, M1, Mf, round(thetat, 1)))
    plt.legend()
    #plt.savefig('Supersonic Inlet g{}M1{}Mf{}.pdf'.format(gamma, M1, Mf, round(thetat, 1)))
    plt.show()
    #plt.close()

    return length, height, thetat

if __name__ == '__main__':
    M1 = 1.95 # Free stream Mach
    Mf = 1.33 # Desired Ending Mach
    L = 10 # inlet distance
    gamma = [1.4, 1.1, 1.8]
    num_div = 19 # Number of divisions

    M1 = 7
    Mf = 1.5
    main(M1, Mf, L, num_div, gamma[0])