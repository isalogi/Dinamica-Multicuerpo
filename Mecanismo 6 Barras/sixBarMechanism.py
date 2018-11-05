# Isabel Cristina López Giraldo
# 6 bar mechanism

import math

import matplotlib.pyplot as plt
import numpy as np


def getJ(a, bars):
    J = np.array([[1,  0,  bars[0]*math.sin(a[0])/2, -1,  0, -bars[1]*math.sin(a[1])/2,  0,  0,              0,  0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  1, -bars[0]*math.cos(a[0])/2,  0, -1,  bars[1]*math.cos(a[1])/2,  0,  0,              0, 0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  1,  0, -bars[1]*math.sin(a[1])/2, -1,  0, -bars[2]*math.sin(a[2])/2,  0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  0,  1, bars[1]*math.cos(a[1])/2,  0, -1,  bars[2]*math.cos(a[2])/2,  0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  1,  0, -bars[2]*math.sin(a[2])/2, -1,  0, -bars[3]*math.sin(a[3])/2,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  0,  1,  bars[2]*math.cos(a[2])/2,  0, -1,  bars[3]*math.cos(a[3])/2,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  0,  0,              0,  1,0, -bars[3]*math.sin(a[3])/2, -1,  0, -bars[4]*math.sin(a[4])/2,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  0,  0,              0,  0,1,  bars[3]*math.cos(a[3])/2,  0, -1,  bars[4]*math.cos(a[4])/2,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  0,  0,              0,  0,0,              0,  1,  0, -bars[4]*math.sin(a[4])/2, -1,  0, -bars[5]*math.sin(a[5])/2],
                  [0,  0,              0,  0,  0,              0,  0,  0,              0,  0,  0,              0,  0,  1,  bars[4]* math.cos(a[4])/2,  0, -1,  bars[5]*math.cos(a[5])/2],
                  [-1,  0, bars[0]*math.sin(a[0])/2,  0,  0,              0,  0,  0,              0,0,  0,              0,  0,  0,              0,  1,  0, -bars[5]*math.sin(a[5])/2],
                  [0, -1, -bars[0]*math.cos(a[0])/2,  0,  0,              0,  0,  0,              0,0,  0,  bars[5]*math.cos(a[3])/2,  0,  0,              0,  0,  1,              0],
                  [1,  0,  bars[0]*math.sin(a[0])/2,  0,  0,              0,  0,  0,              0,0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  1, -bars[0]*math.cos(a[0])/2,  0,  0,              0,  0,  0,              0,0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  0,              1,  0,  0,              0,  0,  0,              0, 0,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  0,  0,              0,1,  0,              0,  0,  0,              0,  0,  0,              0],
                  [0,  0,              0,  0,  0,              0,  0,  0,              0,0,  1,              0,  0,  0,              0,  0,  0,              0],[0,  0,              0,  0,  0,              0,  0,  0,              0,  0,  0,              1,  0,  0, 0, 0, 0, 0]])
                            
    return J

def getFx(x, y, a, bars, xValue, Yvalue) :
    fx = np.array([[x[0] + (-bars[0]/2)*math.cos(a[0]) - (x[1] + (-bars[1]/2)*math.cos(a[1]))],
                   [y[0] + (-bars[0]/2)*math.sin(a[0]) - (y[1] + (-bars[1]/2)*math.sin(a[1]))],
                   [x[1] + (bars[1]/2)*math.cos(a[1]) - (x[2] + (-bars[2]/2)*math.cos(a[2]))],
                   [y[1] + (bars[1]/2)*math.sin(a[1]) - (y[2] + (-bars[2]/2)*math.sin(a[2]))],
                   [x[2] + (bars[2]/2)*math.cos(a[2]) - (x[3] + (-bars[3]/2)*math.cos(a[3]))],
                   [y[2] + (bars[2]/2)*math.sin(a[2]) - (y[3] + (-bars[3]/2)*math.sin(a[3]))],
                   [x[3] + (bars[3]/2)*math.cos(a[3]) - (x[4] + (-bars[4]/2)*math.cos(a[4]))],
                   [y[3] + (bars[3]/2)*math.sin(a[3]) - (y[4] + (-bars[4]/2)*math.sin(a[4]))],
                   [x[4] + (bars[4]/2)*math.cos(a[4]) - (x[5] + (-bars[5]/2)*math.cos(a[5]))],
                   [y[4] + (bars[4]/2)*math.sin(a[4]) - (y[5] + (-bars[5]/2)*math.sin(a[5]))],
                   [x[5] + (bars[5]/2)*math.cos(a[5]) - (x[0] + (bars[0]/2)*math.cos(a[0]))],
                   [y[5] + (bars[5]/2)*math.sin(a[5]) - (y[0] + (bars[0]/2)*math.sin(a[0]))],
                   [x[0]-bars[0]/2*math.cos(a[0])],
                   [y[0]-bars[0]/2*math.sin(a[0])],
                   [a[0]],
                   [x[3]-xValue],
                   [y[3]-Yvalue],
                   [a[3]]])

    return fx

def graphic(x, y, a, bars, xValue, Yvalue) :
    # x1=xa-a/2*math.cos(aa)
    # x2=xa+a/2*math.cos(aa)
    # y1=ya-a/2*math.sin(aa)
    # y2=ya+a/2*math.sin(aa)

    x3 = x[1]-bars[1]/2*math.cos(a[1])
    x4 = x[1]+bars[1]/2*math.cos(a[1])
    x5 = x[2]-bars[2]/2*math.cos(a[2])
    x6 = x[2]+bars[2]/2*math.cos(a[2])
    x7 = x[3]-bars[3]/2*math.cos(a[3])
    x8 = x[3]+bars[3]/2*math.cos(a[3])
    x9 = x[4]-bars[4]/2*math.cos(a[4])
    x10 = x[4]+bars[4]/2*math.cos(a[4])
    x11 = x[5]-bars[5]/2*math.cos(a[5])
    x12 = x[5]+bars[5]/2*math.cos(a[5])

    y3 = y[1]-bars[1]/2*math.sin(a[1])
    y4 = y[1]+bars[1]/2*math.sin(a[1])
    y5 = y[2]-bars[2]/2*math.sin(a[2])
    y6 = y[2]+bars[2]/2*math.sin(a[2])
    y7 = y[3]-bars[3]/2*math.sin(a[3])
    y8 = y[3]+bars[3]/2*math.sin(a[3])
    y9 = y[4]-bars[4]/2*math.sin(a[4])
    y10 = y[4]+bars[4]/2*math.sin(a[4])
    y11 = y[5]-bars[5]/2*math.sin(a[5])
    y12 = y[5]+bars[5]/2*math.sin(a[5])

    # Graficar barras y punto actual
    plt.clf()
    plt.axis((-10, bars[0]+10, -40, 40))
    plt.grid()
    plt.plot([x3, x4, x5, x6, x7, x8, x9, x10, x11, x12], [
                y3, y4, y5, y6, y7, y8, y9, y10, y11, y12], color='red')
    plt.scatter([xValue], [Yvalue], 2, 'b')
    plt.pause(0.0001)

def calculateAngles(a, angles) :
    # Crear vector con angulos actuales
    ang = np.array([a[0], a[1], a[2], a[3], a[4], a[5]])
    ang = ((ang*180)/math.pi)
    # Definir cantidad de decimales a tener en cuenta en matriz de angulos
    ang = np.around(ang, decimals=3)

    # Guardar angulos de cada iteracion en matriz "angulos"
    return np.vstack((angles, ang))

def main():
 
    # Dimensions of the bars
    bar1 = 30
    bar2 = 20
    bar3 = 10
    bar4 = 0.6*(bar1+bar2+bar3+bar2+bar1)  
    bars = [bar4, bar1, bar2, bar3, bar2, bar1]

    # initial conditions
    x = [0, 0, 0, 0, 0, 0]
    y = [0, 0, 0, 0, 0, 0]

    a = [0, 0, 0, 0, 0, math.pi]

    # angle of efector
    pointsx = np.arange(27, 39, 0.5)  # Points x
    pointsy = (pointsx-33)**2  # Points y

    n = len(pointsx)
    angles = np.array(["aa  ", "ab     ", "ac     ", "ad ", "ae    ", "af"])

    xi = np.array([[x[0]], [y[0]], [a[0]], [x[1]], [y[1]], [a[1]], [x[2]], [y[2]], [a[2]],[x[3]], [y[3]], [a[3]], [x[4]], [y[4]], [a[4]], [x[5]], [y[5]], [a[5]]])

    # max iterations
    ei = 100
    itmax = 2000
    iti = 0
    i = 0

    while i < n:
        error = ei
        it = iti
        while (error > 0.00001 and it < itmax):  # Desarrollo de metodo numerico

            fx = getFx(x,y,a,bars,pointsx[i], pointsy[i])
            J = getJ(a, bars)
            JI = np.linalg.pinv(J)
            Xf = xi-np.dot(JI, fx)
            xi = Xf

            x[0] = xi[0, 0]
            x[1] = xi[3, 0]
            x[2] = xi[6, 0]
            x[3] = xi[9, 0]
            x[4] = xi[12, 0]
            x[5] = xi[15, 0]

            y[0] = xi[1, 0]
            y[1] = xi[4, 0]
            y[2] = xi[7, 0]
            y[3] = xi[10, 0]
            y[4] = xi[13, 0]
            y[5] = xi[16, 0]

            a[0] = xi[2, 0]
            a[1] = xi[5, 0]
            a[2] = xi[8, 0]
            a[3] = xi[11, 0]
            a[4] = xi[14, 0]
            a[5] = xi[17, 0]

            error = np.linalg.norm(fx)
            it = it+1

        if (it >= itmax):
            print("Error, cantidad maxima de iteaciones sobrepasada")
            break

        # Pintar la gráfica  
        graphic(x, y, a, bars, pointsx[i], pointsy[i])
        angles = calculateAngles(a, angles)
        i = i+1

    # Crear o sobreescribir archivo .txt con angulos
    np.savetxt("angles.txt", angles, delimiter="   ", fmt="%s")
    # Mostrar ultima grafica
    plt.show()


if __name__ == "__main__":
    main()
