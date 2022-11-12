# Austin Ryan

import matplotlib.pyplot as plt
import numpy as np
import ExpansionFan as EF
import CDNozzles as CD
import IsentropicFlow as IF

def main(angle, start_pts, length, M1, radius=1, clines=True, contour=False,):
    if start_pts % 2 == 1:
        start_pts += 1
    act_pts = int((start_pts-1)/2+1)

    tot_cols = 10
    Mach = np.zeros([tot_cols, act_pts])
    theta = np.zeros([tot_cols, act_pts])
    nu = np.zeros([tot_cols, act_pts])
    kminus = np.zeros([tot_cols, act_pts]) # theta + nu
    kplus = np.zeros([tot_cols, act_pts]) # theta - nu
    slopep = np.zeros([tot_cols, act_pts])
    slopem = np.zeros([tot_cols, act_pts])
    mu = np.zeros([tot_cols, act_pts])

    xs = np.zeros([tot_cols, act_pts])
    ys = np.zeros([tot_cols, act_pts])

    theta[0] = np.linspace(angle, 0, act_pts)

    ys[0] = np.linspace(radius, 0, act_pts)
    xs[0] = np.zeros(act_pts)
    
    for i in range(act_pts):
        Mach[0][i] = M1
        nu[0][i] = EF.PranMeyer(M1)
        kminus[0][i] = theta[0][i] + nu[0][i]
        kplus[0][i] = theta[0][i] - nu[0][i]
        mu[0][i] = np.arcsin(1/Mach[0][i])*180/np.pi

    length_check = False
    i = 1
    while length_check == False:
        if i > len(Mach)-1:
            # Adds more space for additional values if initial column count isn't high enough
            Mach = np.append(Mach, [np.zeros(act_pts)], axis=0)
            theta = np.append(theta, [np.zeros(act_pts)], axis=0)
            nu = np.append(nu, [np.zeros(act_pts)], axis=0)
            kminus = np.append(kminus, [np.zeros(act_pts)], axis=0) # theta + nu
            kplus = np.append(kplus, [np.zeros(act_pts)], axis=0) # theta - nu
            slopep = np.append(slopep, [np.zeros(act_pts)], axis=0)
            slopem = np.append(slopem, [np.zeros(act_pts)], axis=0)
            mu = np.append(mu, [np.zeros(act_pts)], axis=0)
            xs = np.append(xs, [np.zeros(act_pts)], axis=0)
            ys = np.append(ys, [np.zeros(act_pts)], axis=0)

        for j in range(len(Mach[i])):
            if i % 2 == 0: # Even columns
                if j == 0: # Wall point
                    theta[i][j] = angle
                    kplus[i][j] = kplus[i-1][j]
                    nu[i][j] = theta[i][j] - kplus[i][j]
                    Mach[i][j] = EF.inv_PranMeyer(nu[i][j])
                    kminus[i][j] = theta[i][j] + nu[i][j]
                    mu[i][j] = np.arcsin(1/Mach[i][j])*180/np.pi
                    slopep[i-1][j] = np.tan((.5*(theta[i-1][j]+theta[i][j]) + .5*(mu[i-1][j]+mu[i][j]))*np.pi/180)
                    x1 = 0
                    y1 = radius

                    slope1 = angle*np.pi/180

                    xs[i][j] = (ys[i-1][j]-y1-slopep[i-1][j]*xs[i-1][j]+slope1*x1)/(slope1-slopep[i-1][j])
                    ys[i][j] = slopep[i-1][j]*(xs[i][j]-xs[i-1][j])+ys[i-1][j]
                elif j == len(Mach[i])-1: # Axis point
                    theta[i][j] = 0
                    kminus[i][j] = kminus[i-1][j-1]
                    nu[i][j] = kminus[i][j] - theta[i][j]
                    Mach[i][j] = EF.inv_PranMeyer(nu[i][j])
                    kplus[i][j] = theta[i][j] - nu[i][j]
                    mu[i][j] = np.arcsin(1/Mach[i][j])*180/np.pi

                    slopem[i-1][j-1] = np.tan((.5*(theta[i-1][j-1]+theta[i][j]) - .5*(mu[i-1][j-1]+mu[i][j]))*np.pi/180)

                    ys[i][j] = 0
                    xs[i][j] = (ys[i][j]-ys[i-1][j-1])/(slopem[i-1][j-1])+xs[i-1][j-1]
                else: # Interior point
                    theta[i][j] = .5*(kminus[i-1][j-1] + kplus[i-1][j])
                    nu[i][j] = .5*(kminus[i-1][j-1] - kplus[i-1][j])
                    Mach[i][j] = EF.inv_PranMeyer(nu[i][j])
                    kminus[i][j] = kminus[i-1][j-1]
                    kplus[i][j] = kplus[i-1][j]
                    mu[i][j] = np.arcsin(1/Mach[i][j])*180/np.pi

                    slopem[i-1][j-1] = np.tan((.5*(theta[i-1][j-1]+theta[i][j]) - .5*(mu[i-1][j-1]+mu[i][j]))*np.pi/180)
                    slopep[i-1][j] = np.tan((.5*(theta[i-1][j]+theta[i][j]) + .5*(mu[i-1][j]+mu[i][j]))*np.pi/180)

                    xs[i][j] = (ys[i-1][j]-ys[i-1][j-1]-slopep[i-1][j]*xs[i-1][j]+slopem[i-1][j-1]*xs[i-1][j-1])/(slopem[i-1][j-1]-slopep[i-1][j])
                    ys[i][j] = slopep[i-1][j]*(xs[i][j]-xs[i-1][j])+ys[i-1][j]

            elif i % 2 == 1: # Odd columns
                if j == len(Mach[i])-1: # Axis point
                    continue
                else:# Interior point
                    theta[i][j] = .5*(kminus[i-1][j] + kplus[i-1][j+1])
                    nu[i][j] = .5*(kminus[i-1][j] - kplus[i-1][j+1])
                    Mach[i][j] = EF.inv_PranMeyer(nu[i][j])
                    kminus[i][j] = kminus[i-1][j]
                    kplus[i][j] = kplus[i-1][j+1]
                    mu[i][j] = np.arcsin(1/Mach[i][j])*180/np.pi

                    slopem[i-1][j] = np.tan((.5*(theta[i-1][j]+theta[i][j]) - .5*(mu[i-1][j]+mu[i][j]))*np.pi/180)
                    slopep[i-1][j+1] = np.tan((.5*(theta[i-1][j+1]+theta[i][j]) + .5*(mu[i-1][j+1]+mu[i][j]))*np.pi/180)

                    xs[i][j] = (ys[i-1][j+1]-ys[i-1][j]-slopep[i-1][j+1]*xs[i-1][j+1]+slopem[i-1][j]*xs[i-1][j])/(slopem[i-1][j]-slopep[i-1][j+1])
                    ys[i][j] = slopep[i-1][j+1]*(xs[i][j]-xs[i-1][j+1])+ys[i-1][j+1]
        # Stop computing once characteristics pass design length
        if i % 2 == 0: # Even columns
            if xs[i][-1] >= length:
                length_check = True
        elif i % 2 == 1: # Odd columns
            if xs[i][-2] >= length:
                length_check = True
        i += 1

    if clines == True:
        COLOR = 'blue'
        for i in range(1, len(Mach)):
            for j in range(len(Mach[i])):
                if Mach[i][j] == 0:
                    continue
                elif i % 2 == 0: # Even column
                    if j == len(Mach[i])-1: # Axis point
                        plt.plot([xs[i-1][j-1], xs[i][j]], [ys[i-1][j-1], ys[i][j]], color=COLOR)
                    elif j == 0: # Wall point
                        plt.plot([xs[i-1][j], xs[i][j]], [ys[i-1][j], ys[i][j]], color=COLOR)
                    else: # Interior point
                        plt.plot([xs[i-1][j-1], xs[i][j]], [ys[i-1][j-1], ys[i][j]], color=COLOR)
                        plt.plot([xs[i-1][j], xs[i][j]], [ys[i-1][j], ys[i][j]], color=COLOR)

                if i % 2 == 1: # Odd column
                    if j == len(Mach[i])-1: # Axis point
                        continue
                    else: # Interior point
                        plt.plot([xs[i-1][j+1], xs[i][j]], [ys[i-1][j+1], ys[i][j]], color=COLOR)
                        plt.plot([xs[i-1][j], xs[i][j]], [ys[i-1][j], ys[i][j]], color=COLOR)
    
    maxx = 0
    maxy = 0
    maxM = 0
    for i in xs:
        for j in i:
            if j > maxx:
                maxx = j
    for i in ys:
        for j in i:
            if j > maxy:
                maxy = j
    for i in Mach:
        for j in i:
            if j > maxM:
                maxM = j
    md = angle*np.pi/180
    maxx2 = maxx+.5
    yd = np.tan(md)*maxx2+radius

    plt.plot([0, maxx2], [0, 0], linestyle='--', color='gray', label='Centerline')
    plt.plot([0, maxx2], [radius, yd], color='red', label='Duct Contour')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Diverging Duct Method of Characteristics')
    #plt.ylim([-.1, maxx2])
    plt.legend()

    if contour == True:
        indices = []
        for i in range(1, len(xs)):
            for j in range(len(xs[i])):
                if xs[i][0] == 0:
                    indices.append(i)
                elif xs[i][j] == 0:
                    xs[i][j] = xs[i, j-1]
                    Mach[i][j] = Mach[i, j-1]
        xs = np.delete(xs, indices, axis=0)
        ys = np.delete(ys, indices, axis=0)
        Mach = np.delete(Mach, indices, axis=0)

        plt.contourf(xs, ys, Mach, levels=np.around(np.linspace(M1, maxM, 30), 3), cmap='inferno')
        plt.colorbar(orientation='vertical').set_label('Mach Number')
    #plt.ylim([-.1, maxx2])
    #plt.savefig('MOC Diverging Duct {}.pdf'.format(int(start_pts-1)))
    plt.close()
    plt.show()
    
    return xs, ys, Mach

def analytic_solution(half_angle, length, M1, dist, MachA, start_pts):
    for i in range(len(dist)):
        indices = []
        for j in range(len(dist[i])):
            if dist[i][j] == 0 and j != 0:
                indices.append(j)
        dist[i] = np.delete(dist[i], indices)
        MachA[i] = np.delete(MachA[i], indices)
        plt.plot(dist[i], MachA[i], label='MOC {} Points'.format(start_pts[i]))

    xs = np.linspace(0, length, 100)
    ys = np.tan(half_angle*np.pi/180)*xs+1
    Astar = CD.SonicAreaRatio(M1)
    Mach = np.zeros(len(ys))
    for i in range(len(ys)):
        Mach[i] = CD.MachfromAreaRatio(ys[i]/Astar)

    plt.plot(xs, Mach, label='Analytic Result')
    
    #plt.ylim([1, max(Mach)+.5])
    plt.ylabel('Mach Number')
    plt.xlabel('Horiztonal Distance')
    plt.title('Mach Number Along Nozzle')
    plt.legend()
    # plt.savefig('Diverging Duct Mach Number.pdf')
    #plt.close()
    plt.show() 

    diff = np.zeros(len(MachA))
    for i in range(len(MachA)):
        diff[i] = abs(MachA[i][-1] - Mach[-1])

    plt.plot(start_pts, diff, label='$L_1$ norm')
    plt.title('$L_1$ Norm of Final Mach Number')
    plt.ylabel('$L_1$')
    plt.xlabel('Start Points')
    plt.legend()
    # plt.savefig('L1 Norm.pdf')
    #plt.close()
    plt.show()

def presstemp(xs, Mach, start_points, length, GAMMA=1.4):
    for i in range(len(xs)):
        indices = []
        for j in range(len(xs[i])):
            if xs[i][j] == 0 and j != 0:
                indices.append(j)
        xs[i] = np.delete(xs[i], indices)
        Mach[i] = np.delete(Mach[i], indices)
 
    for run in range(len(xs)):
        temp = np.zeros(len(xs[run]))
        for i in range(len(xs[run])):
            temp[i] = IF.TRatio(Mach[run][i], GAMMA)
        plt.plot(xs[run], temp, label='T/T0 {}'.format(start_points[run]))

    xsa = np.linspace(0, length, 100)
    ysa = np.tan(half_angle*np.pi/180)*xsa+1
    Astar = CD.SonicAreaRatio(M1)
    tempa = np.zeros(len(ysa))
    for i in range(len(ysa)):
        tempa[i] = IF.TRatio(CD.MachfromAreaRatio(ysa[i]/Astar))
    plt.plot(xsa, tempa, label='Analytic Result')

    plt.title('Temperature Ratio Along the Nozzle')
    plt.xlabel('X Direction')
    plt.ylabel('Temperature Ratio ($T/T_0$)')
    plt.legend()
    # plt.savefig('Diverging Duct Temperature Ratio.pdf')
    plt.close()
    plt.show()

    for run in range(len(xs)):
        press = np.zeros(len(xs[run]))
        for i in range(len(xs[run])):
            press[i] = IF.pRatio(Mach[run][i], GAMMA)
        plt.plot(xs[run], press, label='P/P0 {}'.format(start_points[run]))

    xsa = np.linspace(0, length, 100)
    ysa = np.tan(half_angle*np.pi/180)*xsa+1
    Astar = CD.SonicAreaRatio(M1)
    pressa = np.zeros(len(ysa))
    for i in range(len(ysa)):
        pressa[i] = IF.pRatio(CD.MachfromAreaRatio(ysa[i]/Astar))
    plt.plot(xsa, pressa, label='Analytic Result')

    plt.title('Pressure Ratio Along the Nozzle')
    plt.xlabel('X Direction')
    plt.ylabel('Pressure Ratio ($P/P_0$)')
    plt.legend()
    #plt.savefig('Diverging Duct Pressure Ratio.pdf')
    plt.close()
    plt.show()

if __name__ == '__main__':
    # STARTING "RADIUS" IS 1 
    full_angle = 18 # degrees
    half_angle = full_angle/2
    start_points = [3,5,7]
    M1 = 1.4349
    length = 10
    Mach1 = []
    xs1 = []
    ys1 = []

    for start_pts in start_points:
        xs, ys, Mach = main(half_angle, start_pts, length, M1)
        Mach1.append(Mach[:,0])
        ys1.append(ys[:, 0])
        xs1.append(xs[:, 0])

    presstemp(xs1, Mach1, start_points, length)
    analytic_solution(half_angle, length, M1, xs1, Mach1, start_points)