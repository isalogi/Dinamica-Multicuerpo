# Isabel Cristina López Giraldo
# 2 bar mechanism SCARA(Natural coordinates)

import math
import sys
import numpy as np
import matplotlib.pyplot as plt

def main():

    #Chasis position
    xa = 0
    ya = 0

    ##Dimensions of the bars
    a = 12
    b = 6

    #Initial variables
    x1 = 4
    y1 = 0
    x2 = 4
    y2 = 4
    
    i=6.28319
    counter=0 #counter

    #Newton-Raphson
    fig=plt.figure()
    plt.ion()

    for t in np.arange(0,float(i), 0.01):

        counter+=1
        xi=np.array([[x1], [y1], [x2], [y2]])
        residue=1

        while residue > 0.000001:

            #funtion
            e1= (x1-xa)**2+(y1-ya)**2-a**2
            e2= (x2-x1)**2+(y2-y1)**2-b**2
            e3= x2-6-math.sin(math.pi*2*t)
            e4=y2-5-math.sin(math.pi*(4*t-2))
            fx = np.array([[e1],[e2],[e3],[e4]])

            #The Jacobian
            Jx= np.array([[2*(x1-xa), 2*(y1-y2), 0, 0],
                         [-2*(x2-x1), -2*(y2-y1), 2*(x2-x1), 2*(y2-y1)],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        
            Jin=np.linalg.inv(Jx)

            xf = xi - np.dot(Jin,fx)
            xi = xf
            x1 = xi[0,0]
            y1 = xi[1,0]
            x2 = xi[2,0]
            y2 = xi[3,0]
            residue = xf - xi
            residue = np.linalg.norm(fx)

        plt.clf()
        plt.plot([xa,x1,x2],[ya,y1,y2],color='blue') 
        plt.plot([6+math.sin(math.pi*2*t)],[5+math.sin(math.pi*(4*t-2))],'ro')
        fig.suptitle('SCARA MECHANISM')
        plt.xlabel('x (arbitrary units)')
        plt.ylabel('y (arbitrary units)')
        plt.axis((-2,15,-2,15))
        plt.draw()
        plt.pause(0.0001)
        
    plt.show()


if __name__ == "__main__":
    main()