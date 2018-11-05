# Isabel Cristina López Giraldo
# 4 bar mechanism (Natural coordinates)

import math
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():

    #Lenghts
    # la=sqrt((x1-xa)^2+(y1-ya)^2);
    a=2
    b=4
    c=5
    d=6

    #Initial variables
    phi1=0  #Angle from ground
    phi2initial=0
    omega2initial=1
    counter=0

    #Coordinates of the bars
    x1=2
    x2=3
    xa=0
    xb=d*math.cos(phi1)+xa

    y1=0.3
    y2=5
    ya=0
    yb=d*math.sin(phi1)+ya

    xi=[[x1],[y1],[x2],[y2]]

    i=2*math.pi
    #Newton-Raphson
    fig=plt.figure()
    plt.ion()

    for phi2 in np.arange(0,float(i), 0.05):

        counter+=1
        residue=1

        while residue > 0.00001:
            
            #Time
            t=(phi2-phi2initial)/omega2initial
            #tp=((-omega2initial)+(sqrt(omega2initial^2-2*alpha2initial*(phi2initial-phi2)))/alpha2initial);
            #tn=((-omega2initial)-(sqrt(omega2initial^2-2*alpha2initial*(phi2initial-phi2)))/alpha2initial);

            x1=xi[0][0]
            y1=xi[1][0]
            x2=xi[2][0]
            y2=xi[3][0]

            #Function
            #x1-xa-a*cos((0.5*alpha2initial*t^2)+(omega2initial*t)+phi2initial)

            e1=(x1-xa)**2+(y1-ya)**2-a**2
            e2=(x2-x1)**2+(y2-y1)**2-b**2
            e3=(x2-xb)**2+(y2-yb)**2-c**2
            e4=x1-xa-a*math.cos((omega2initial*t)+phi2initial)
            fx =[[e1],[e2],[e3],[e4]]

            #The Jacobian
            Jx=[[2*(x1-xa),2*(y1-ya),0,0],[-2*(x2-x1),-2*(y2-y1),2*(x2-x1),2*(y2-y1)],[0,0,2*(x2-xb),2*(y2-yb)],[1,0,0,0]]

            xf = xi - np.dot(np.linalg.pinv(Jx),fx)
            xi = xf
            residue = xf - xi
            residue = np.linalg.norm(fx)

        plt.clf()
        plt.plot(xa,ya,'ok')
        plt.plot(xb,yb,'ok')
        plt.plot(xf[0][0],xf[1][0],'or')
        plt.plot(xf[2][0],xf[3][0],'ob')
        plt.plot([xa,xf[0][0]],[ya,xf[1][0]],'m')
        plt.plot([xf[0][0],xf[2][0]],[xf[1][0],xf[3][0]],'m')
        plt.plot([xf[2][0],xb],[xf[3][0],yb],'m')
        fig.suptitle('4 bar linkage')
        plt.xlabel('x (arbitrary units)')
        plt.ylabel('y (arbitrary units)')
        plt.axis((-3,8,-4,6))
        plt.pause(0.005)
        plt.ioff()
    plt.show()


if __name__ == "__main__":
    main()
