# Isabel Cristina López Giraldo
# 4 bar mechanism (kinematic analysis)

import math
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    
    #Dimensions of the bars
    a=2.23
    b=4.47
    c=6
    d=8

    #Initial variables
    phi2initial=0.0
    omega2initial=1.0
    alpha2initial=1.0
    counter=0 #counter
    i=2*math.pi
    #Newton-Raphson
    fig=plt.figure()
    plt.ion()

    for phi2 in np.arange(0,float(i), 0.05):
        counter+=1
        xi=[[1],[2.2]]
        residue=1

        while residue> 0.00001:
            #funtion
            e1= a*math.cos(phi2)+b*math.cos(xi[0][0])-c*math.cos(xi[1][0])-d
            e2= a*math.sin(phi2)+b*math.sin(xi[0][0])-c*math.sin(xi[1][0])
            fx = [[e1],[e2]]

            #The Jacobian
            j1=-b*math.sin(xi[0][0])
            j2=c*math.sin(xi[1][0])
            j3= b*math.cos(xi[0][0])
            j4= -c*math.cos(xi[1][0])
            Jx =[[j1,j2],[j3,j4]]

            xf = xi - np.linalg.pinv(Jx)*fx
            xi = xf
            residue = xf - xi
            residue = np.linalg.norm(fx)

        #Aceleration
        alpha2=alpha2initial
        
        #Angular velocity as a function of rotation and not time
        omega2=math.sqrt(2*alpha2*(phi2-phi2initial)+omega2initial**2)
        
        #Velocity
        vfinal=-np.linalg.pinv(Jx)*[[-a*omega2*math.sin(phi2)],[a*omega2*math.cos(phi2)]]

        dj1=-b*vfinal[0][0]*math.cos(xf[0][0])
        dj2=c*vfinal[1][0]*math.cos(xf[1][0])
        dj3=-b*vfinal[0][0]*math.sin(xf[0][0])
        dj4=c*vfinal[1][0]*math.sin(xf[1][0])
        dJx=[[dj1,dj2],[dj3,dj4]]

        atan=[[-a*alpha2*math.sin(phi2)],[a*alpha2*math.cos(phi2)]]
        
        anorm=[[-a*omega2**2*math.cos(phi2)],[-a*omega2**2*math.sin(phi2)]]
        afinal=-np.linalg.pinv(Jx)*(dJx*vfinal+atan+anorm)

        X=[0,a*math.cos(phi2),d+c*math.cos(xf[1][0]),d]
        Y=[0,a*math.sin(phi2),c*math.sin(xf[1][0]),0]

        plt.plot(X,Y,'b')
        fig.suptitle('4 bar linkage')
        plt.xlabel('x (arbitrary units)')
        plt.ylabel('y (arbitrary units)')
        plt.axis((-2,8,-4,7))
        plt.draw()
        plt.pause(0.00001)
        plt.clf()
    plt.show()


if __name__ == "__main__":
    main()
