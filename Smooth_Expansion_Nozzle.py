

import matplotlib.pyplot as plt
import ExpansionFan as expfan
import numpy as np

def MinLength_nozzle(Me, num_waves, throat_radius, gamma=1.4):
    """Creates a minimum length nozzle contour using the method of charcteristics
    Parameters
    ----------
    Me : float
        Exit mach supersonic
    num_waves : float
        Number of mach waves to simulate"""

    downstream_radius = throat_radius * .4

    # Loops over plus characteristics of the kernel as a row then moves to the next plus characteristic
    x = np.zeros((num_waves + 2, num_waves)) # x coordinates of waves; the last row contains the wall points
    y = np.zeros((num_waves + 2, num_waves)) # y coordinates of waves; the last row contains the wall points
    theta = np.zeros((num_waves + 1, num_waves)) # Angle (deg) between velcotiy vector and horizontal
    nu = np.zeros((num_waves + 1, num_waves)) # Prandtl-meyer angle (deg)
    Mach = np.zeros((num_waves + 1, num_waves)) # Mach number at point
    mu = np.zeros((num_waves + 1, num_waves)) # Angle (deg) between velocity vector and tangent line to both characteristic lines at the point of analysis
    Kminus = np.zeros((num_waves + 1, num_waves)) # Changes along a row, only a da of unique values are present
    Kplus = np.zeros((num_waves + 1, num_waves)) # Constant along a row
    Cplus = np.zeros((num_waves + 1, num_waves)) # Slope of positive sloped characterisitc lines AKA left running, _plus
    Cminus = np.zeros((num_waves + 1, num_waves)) # Slope of negative sloped characterisitc lines AKA right running, _minus
    
    # Last row in x, y lists are the final wall points of converging section
    # First row in x, y lists are the initial wall points of the downstream radius of the throat

    theta_max = expfan.PranMeyer(Me, gamma=gamma, deg=False)/2 # Returns max wall angle in radians

    # Setting initial values, theta and nu are equal
    theta[0] = np.linspace(0, theta_max, len(theta[0])+1)[1:]
    nu[0] = np.linspace(0, theta_max, len(theta[0])+1)[1:]
    
    for i in range(num_waves): # Calculating the conditions at the throat radius poitns
        Mach[0][i]   = expfan.inv_PranMeyer(nu[0][i], gamma=gamma, deg=False)
        mu[0][i]     = np.arcsin(1/Mach[0][i])
        Kminus[0][i] = theta[0][i] + nu[0][i]
        Kplus[0][i]  = theta[0][i] - nu[0][i]
        y[0][i] = throat_radius + downstream_radius*(1-np.cos(theta[0][i]))
        x[0][i] = downstream_radius*np.sin(theta[0][i])

    # Setting compatibility values because they are constant along a characteristic
    Kminus[1] = Kminus[0]
    Kplus[1][:num_waves-1+1] = -Kminus[0][1-1]
    for i in range(2, num_waves+1):
        Kplus[i][:num_waves-i+1] = -Kminus[0][i-1]
        Kminus[i][0:-(i-1)] = Kminus[0][i-1:]

    # First centerline point
    theta[1][0] = (Kplus[1][0] + Kminus[1][0]) / 2
    nu[1][0] = (Kminus[1][0] - Kplus[1][0]) / 2
    Mach[1][0] = expfan.inv_PranMeyer(nu[1][0], gamma=gamma, deg=False)
    mu[1][0] = np.arcsin(1/Mach[1][0])
    Cminus[1][0] = (theta[0][0] + theta[1][0] - mu[0][0] - mu[1][0]) / 2
    # Solving for coordiantes with new point y coord known
    y[1][0] = 0
    x[1][0] = -y[0][0]/np.tan(Cminus[1][0]) + x[0][0]
    
    k = 1
    for j in range(1, num_waves): # Moving along first left running characteristic
        theta[k][j] = (Kplus[k][j] + Kminus[k][j]) / 2
        nu[k][j] = (Kminus[k][j] - Kplus[k][j]) / 2
        Mach[k][j] = expfan.inv_PranMeyer(nu[k][j], gamma=gamma, deg=False)
        mu[k][j] = np.arcsin(1/Mach[k][j])
        Cplus[k][j]  = (theta[k][j-1] + mu[k][j-1] + theta[k][j] + mu[k][j]) / 2 # Slope of the line going into the new point
        Cminus[k][j] = (theta[k-1][j] - mu[k-1][j] + theta[k][j] - mu[k][j]) / 2 # Slope of the line going into the new point
        # Solving for coordinates
        a = np.tan(Cplus[k][j])*x[k][j-1] - np.tan(Cminus[k][j])*x[k-1][j]
        b = np.tan(Cplus[k][j]) - np.tan(Cminus[k][j])
        x[k][j] = (a + y[k-1][j] - y[k][j-1])/b
        y[k][j] = np.tan(Cminus[k][j])*(x[k][j] - x[k-1][j]) + y[k-1][j]

    # First Wall Point
    a = np.tan((theta[0][-1]+theta[k][j])/2)
    b = np.tan(Cplus[k][j])
    x[-1][0] = (a*x[0][-1] - b*x[k][j] + y[k][j] - y[0][-1]) / (a-b)
    y[-1][0] = a*(x[-1][0]-x[0][-1]) + y[0][-1]

    for k in range(2,num_waves+1):
        for j in range(num_waves-k+1):
            if j == 0:
                theta[k][j] = (Kplus[k][j] + Kminus[k][j]) / 2
                nu[k][j] = (Kminus[k][j] - Kplus[k][j]) / 2
                Mach[k][j] = expfan.inv_PranMeyer(nu[k][j], gamma=gamma, deg=False)
                mu[k][j] = np.arcsin(1/Mach[k][j])
                Cminus[k][j] = (theta[k-1][j+1] - mu[k-1][j+1] + theta[k][j] - mu[k][j]) / 2 # Slope of the line going into the new point
                Cplus[k][j]  = (theta[k][j-1] + mu[k][j-1] + theta[k][j] + mu[k][j]) / 2 # Slope of the line going into the new point
                # Solving for coordinates
                y[k][j] = 0
                x[k][j] = -y[k-1][j+1]/np.tan(Cminus[k][j]) + x[k-1][j+1]
                if j == num_waves-k: # Wall points
                    a = np.tan((theta[k-1][j+1] + theta[k][j])/2)
                    b = np.tan(Cplus[k][j])
                    x[-1][k-1] = (a*x[-1][k-2] - b*x[k][j] + y[k][j] - y[-1][k-2]) / (a-b)
                    y[-1][k-1] = a*(x[-1][k-1]-x[-1][k-2]) + y[-1][k-2]
            else:
                theta[k][j] = (Kplus[k][j] + Kminus[k][j]) / 2
                nu[k][j] = (Kminus[k][j] - Kplus[k][j]) / 2
                Mach[k][j] = expfan.inv_PranMeyer(nu[k][j], gamma=gamma, deg=False)
                mu[k][j] = np.arcsin(1/Mach[k][j])
                Cplus[k][j]  = (theta[k][j-1] + mu[k][j-1] + theta[k][j] + mu[k][j]) / 2 # Slope of the line going into the new point
                Cminus[k][j] = (theta[k-1][j+1] - mu[k-1][j+1] + theta[k][j] - mu[k][j]) / 2 # Slope of the line going into the new point
                # Solving for coordinates
                a = np.tan(Cplus[k][j])*x[k][j-1] - np.tan(Cminus[k][j])*x[k-1][j+1]
                b = np.tan(Cplus[k][j]) - np.tan(Cminus[k][j])
                x[k][j] = (a + y[k-1][j+1] - y[k][j-1])/b
                y[k][j] = np.tan(Cminus[k][j])*(x[k][j] - x[k-1][j+1]) + y[k-1][j+1]
                if j == num_waves-k: # Wall points, always calculated last in a left characteristic
                    a = np.tan((theta[k-1][j+1] + theta[k][j])/2)
                    b = np.tan(Cplus[k][j])
                    x[-1][k-1] = (a*x[-1][k-2] - b*x[k][j] + y[k][j] - y[-1][k-2]) / (a-b)
                    y[-1][k-1] = a*(x[-1][k-1]-x[-1][k-2]) + y[-1][k-2]

    # Throat radius to first left running
    for i in range(num_waves):
        # minus/left characteristics
        if i == 0:
            plt.plot([x[0][i], x[1][i]], [y[0][i], y[1][i]], marker='.', color='#FFD700', label='Left Characteristics')
        else:
            plt.plot([x[0][i], x[1][i]], [y[0][i], y[1][i]], marker='.', color='#FFD700')
    for j in range(1, num_waves):
        # plus/right characteristics
        if j == 1:
            plt.plot([x[1][j-1], x[1][j]], [y[1][j-1], y[1][j]], marker='.', color='#0028FF', label='Right Characteristics')
        else:
            plt.plot([x[1][j-1], x[1][j]], [y[1][j-1], y[1][j]], marker='.', color='#0028FF')

    # Remaining expansion region plotting
    for i in range(2, num_waves+1):
        for j in range(len(theta[i])-i+1):
            # minus characteristics
            plt.plot([x[i-1][j+1], x[i][j]], [y[i-1][j+1], y[i][j]], marker='.', color='#FFD700')
            if j != 0:
                # plus characteristics
                plt.plot([x[i][j-1], x[i][j]], [y[i][j-1], y[i][j]], marker='.', color='#0028FF')

    # Wall points
    plt.plot(x[0], y[0], label='Initial Radius', color='#8000FF') # Throat radius
    plt.plot(x[-1], y[-1], marker='o')

    plt.xlim([0, np.max([np.max(x), np.max(y)])])
    plt.ylim([0, np.max([np.max(x), np.max(y)])])
    plt.legend()
    plt.show()
    plt.savefig('Smooth Expanion Nozzle.png')

if __name__ == '__main__':

    Me = 3.8486
    num_waves = 15
    throat_area = .042; #m^2
    throat_radius = np.sqrt(throat_area/np.pi) #m
    throat_radius = 1
    downstream_radius = throat_radius * .4 #m
    gamma = 1.21

    MinLength_nozzle(Me, num_waves, throat_radius, gamma=gamma)