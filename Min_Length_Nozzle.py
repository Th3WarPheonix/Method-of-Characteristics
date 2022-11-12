

import matplotlib.pyplot as plt
import ExpansionFan as expfan
import numpy as np

def MinLength_nozzle(Me, num_waves, angle0=.375, gamma=1.4):
    """Creates a minimum length nozzle contour using the method of charcteristics
    Parameters
    ----------
    Me : float
        Exit mach supersonic
    num_waves : float
        Number of mach waves to simulate"""

    # Loops over plus characteristics of the kernel as a row then moves to the next plus characteristic
    x = np.zeros((num_waves + 2, num_waves)) # x coordinates of waves; the last row contains the wall points
    y = np.zeros((num_waves + 2, num_waves)) # y coordinates of waves; the last row contains the wall points
    theta = np.zeros((num_waves + 1, num_waves)) # Angle (deg) between velcotiy vector and horizontal
    nu = np.zeros((num_waves + 1, num_waves)) # Prandtl-meyer angle (deg)
    Mach = np.zeros((num_waves + 1, num_waves)) # Mach number at point
    mu = np.zeros((num_waves + 1, num_waves)) # Angle (deg) between velocity vector and tangent line to both characteristic lines at the point of analysis
    Kminus = np.zeros((num_waves + 1, num_waves)) # Changes along a row, only a da of unique values are present
    Kplus = np.zeros((num_waves + 1, num_waves)) # Constant along a row
    Cplus = np.zeros((num_waves + 1, num_waves)) # Slope of positive sloped characterisitc lines AKA right running, _plus
    Cminus = np.zeros((num_waves + 1, num_waves)) # Slope of negative sloped characterisitc lines AKA left running, _minus
    
    theta_max = expfan.PranMeyer(Me, gamma=gamma)/2 # Returns max wall angle in degrees

    # Setting initial values
    theta[0] = np.linspace(angle0, theta_max, len(theta[0]))
    nu[0] = np.linspace(angle0, theta_max, len(theta[0]))
    
    for i in range(num_waves):
        Mach[0][i]   = expfan.inv_PranMeyer(nu[0][i], gamma=gamma)
        mu[0][i]     = np.arcsin(1/Mach[0][i])*180/np.pi
        Kminus[0][i] = theta[0][i] + nu[0][i]
        Kplus[0][i]  = theta[0][i] - nu[0][i]
        # print(Mach[0][i], mu[0][i], Kminus[0][i], Kplus[0][i])
        
    x[0] = 0
    y[0] = 1
    #print(theta[0])
    #print(nu[0])
    #print(Mach[0])
    #print(mu[0])

    point = 1
    for j in range(num_waves):
        Kplus[1] = 0
        Kminus[1][j] = Kminus[0][j]
        theta[1][j] = 1/2*(Kminus[1][j] + Kplus[1][j])
        nu[1][j] = 1/2*(Kminus[1][j] - Kplus[1][j])
        Mach[1][j] = expfan.inv_PranMeyer(nu[1][j], gamma=gamma)
        mu[1][j] = np.arcsin(1/Mach[1][j])*180/np.pi

        Cminus[1][j] = .5*(theta[0][j]+theta[1][j]) - .5*(mu[0][j]+mu[1][j])
        Cminus[1][j] = np.tan(Cminus[1][j]*np.pi/180)
        if j != 0:
            Cplus[1][j] = (theta[1][j-1]+theta[1][j])/2 + (mu[1][j-1]+mu[1][j])/2
            Cplus[1][j] = np.tan(Cplus[1][j]*np.pi/180)
            a = -y[1][j-1] + y[0][j]
            b = Cplus[1][j]*x[1][j-1] - Cminus[1][j]*x[0][j]
            x[1][j] = (a+b)/(Cplus[1][j] - Cminus[1][j])
            y[1][j] = Cplus[1][j]*(x[1][j] - x[1][j-1]) + y[1][j-1]
        else:
            x[1][j] = -1/Cminus[1][j]
            y[1][j] = 0
        # print(point, j, Kminus[1][j], Kplus[1][j], theta[1][j], nu[1][j], Mach[1][j], mu[1][j], np.arctan(Cminus[1][j])*180/np.pi, np.arctan(Cplus[1][j])*180/np.pi)
        # print(x[1][j], y[1][j])
        if j == num_waves - 1: # Wall points
            slope = (theta_max+theta[1][j])/2
            slope = np.tan(slope*np.pi/180)
            a = -y[0][0]+y[1][j]
            b = slope*x[0][0] - Cplus[1][j]*x[1][j]
            x[-1][0] = (a+b)/(slope - Cplus[1][j])
            y[-1][0] = Cplus[1][j]*(x[-1][0] - x[1][j]) + y[1][j]
        point += 1
    point += 1
    
    for i in range(2, num_waves+1):
        Kplus[i] = -Kminus[0][i-1]
        for j in range(len(theta[i])-i+1):
            Kminus[i][j] = Kminus[i-1][j+1]
            theta[i][j] = 1/2*(Kminus[i][j] + Kplus[i][j])
            nu[i][j] = 1/2*(Kminus[i][j] - Kplus[i][j])
            Mach[i][j] = expfan.inv_PranMeyer(nu[i][j], gamma=gamma)
            mu[i][j] = np.arcsin(1/Mach[i][j])*180/np.pi
            Cminus[i][j] = (theta[i-1][j+1]+theta[i][j])/2 - (mu[i-1][j+1]+mu[i][j])/2
            Cminus[i][j] = np.tan(Cminus[i][j]*np.pi/180)
            if j != 0:
                Cplus[i][j] = (theta[i][j-1]+theta[i][j])/2 + (mu[i][j-1]+mu[i][j])/2
                Cplus[i][j] = np.tan(Cplus[i][j]*np.pi/180)
                a = -y[i][j-1]+y[i-1][j+1]
                b = Cplus[i][j]*x[i][j-1] - Cminus[i][j]*x[i-1][j+1]
                x[i][j] = (a+b)/(Cplus[i][j] - Cminus[i][j])
                y[i][j] = Cplus[i][j]*(x[i][j] - x[i][j-1]) + y[i][j-1]
            #print(point, j, Kminus[i][j], Kplus[i][j], theta[i][j], nu[i][j], Mach[i][j], mu[i][j], np.arctan(Cminus[i][j])*180/np.pi, np.arctan(Cplus[i][j])*180/np.pi)            
            elif j == 0:
                #print(np.arctan(Cminus[i][j])*180/np.pi)
                y[i][j] = 0
                x[i][j] = -y[i-1][j+1]/Cminus[i][j]+x[i-1][j+1]
            #print(x[i][j], y[i][j])
            if i == num_waves and j == 0:
                slope = (theta[i-1][j+1] + theta[i][j])/2
                slope = np.tan(slope*np.pi/180)
                a = -y[-1][i-2]+y[i][j]
                Cplus[i][j] = np.tan(mu[i][j]*np.pi/180)
                b = slope*x[-1][i-2] - Cplus[i][j]*x[i][j]
                x[-1][i-1] = (a+b)/(slope - Cplus[i][j])
                y[-1][i-1] = Cplus[i][j]*(x[-1][i-1] - x[i][j]) + y[i][j]
            elif j == num_waves - i:
                slope = (theta[i-1][j+1] + theta[i][j])/2
                slope = np.tan(slope*np.pi/180)
                a = -y[-1][i-2]+y[i][j]
                b = slope*x[-1][i-2] - Cplus[i][j]*x[i][j]
                x[-1][i-1] = (a+b)/(slope - Cplus[i][j])
                y[-1][i-1] = Cplus[i][j]*(x[-1][i-1] - x[i][j]) + y[i][j]

            point += 1
        point += 1

    # Initial expansion
    for i in range(num_waves):
        # minus characteristics
        plt.plot([x[0][i], x[1][i]], [y[0][i], y[1][i]], color='orchid')
    for j in range(1, num_waves):
        # plus characteristics
        plt.plot([x[1][j-1], x[1][j]], [y[1][j-1], y[1][j]], color='red')
    
    for i in range(2, num_waves+1):
        for j in range(len(theta[i])-i+1):
            # minus characteristics
            plt.plot([x[i-1][j+1], x[i][j]], [y[i-1][j+1], y[i][j]], color='gold')
            if j != 0:
                # plus characteristics
                plt.plot([x[i][j-1], x[i][j]], [y[i][j-1], y[i][j]], color='teal')
    
    # Throat to first wall point
    plt.plot([x[-1][0], x[0][0]], [y[-1][0], y[0][0]], color='blue')
    plt.plot([x[-1][0], x[1][-1]], [y[-1][0], y[1][-1]], color='green')

    # Wall points
    for i in range(2, num_waves+1):
        k = num_waves - i
        plt.plot([x[-1][i-1], x[i][k]], [y[-1][i-1], y[i][k]], color='pink')
        plt.plot([x[-1][i-1], x[-1][i-2]], [y[-1][i-1], y[-1][i-2]], color='rebeccapurple')
    

    plt.title('Method of Characteristics Nozzle Contour')
    plt.xlabel('Distance from Throat')
    plt.ylabel('Vertical Distance')
    # plt.legend()
    plt.savefig('MOC Nozzle Design.png')
    # plt.close()
    plt.show()
    print(Kminus)
    print(Kplus)
if __name__ == '__main__':

    Me = 3.8486
    num_waves = 15
    r_t = 1
    gamma = 1.21

    MinLength_nozzle(Me, num_waves, gamma=gamma)