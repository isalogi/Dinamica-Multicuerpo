# Isabel Cristina López Giraldo
# 2 bar mechanism SCARA

import math
import sys
import numpy as np
import matplotlib.pyplot as plt


def main():
    # Dimensions of the bars
    a = 12
    b = 6

    x0 = 0
    y0 = 0

    i = 6.28319

    counter = 0  # counter

    phi1 = math.pi/2
    phi2 = 0.0
    phi3 = math.pi/4
    c = 0.0
    xt = []
    yt = []
    # Newton-Raphson
    fig = plt.figure()
    plt.ion()

    for t in np.arange(0, float(i), 0.01):

        counter += 1
        xi = [[phi1], [phi2], [phi3], [c]]

        residue = 1

        while residue > 0.000001:

            # funtion
            e1 = a*math.cos(phi1)+b*math.cos(phi2)-c*math.cos(phi3)
            e2 = a*math.sin(phi1)+b*math.sin(phi2)-c*math.sin(phi3)
            e3 = c*math.cos(phi3)-3-math.sin(math.pi*2*t)
            e4 = c*math.sin(phi3)-2-math.sin(math.pi*(4*t-2))
            fx = np.array([[e1], [e2], [e3], [e4]])

            # The Jacobian
            Jx = np.array([[-math.sin(xi[0][0])*a, -math.sin(xi[1][0])*b, math.sin(xi[2][0])*c, -math.cos(xi[2][0])],
                           [math.cos(xi[0][0])*a, math.cos(xi[1][0])*b, -
                            math.cos(xi[2][0])*c, -math.sin(xi[2][0])],
                           [0, 0, -math.sin(xi[2][0])*c, math.cos(xi[2][0])],
                           [0, 0, math.cos(xi[2][0]), math.sin(xi[2][0])]])

            Jin = np.linalg.pinv(Jx)

            xf = xi - np.dot(Jin, fx)
            xi = xf
            
            residue = np.linalg.norm(fx)

            x1 = x0 + a*math.cos(xi[1][0])
            y1 = y0 + a*math.sin(xi[1][0])
            x2 = x1 + b*math.cos(xi[2][0])
            y2 = y1 + b*math.sin(xi[2][0])

        plt.plot([x0, x1], [y0, y1], 'r-o')
        plt.plot([x1, x2], [y1, y2], 'r-o')
        plt.plot([x0, x2], [y0, y2], 'k--o')
        plt.plot(xt, yt, 'b')
        fig.suptitle('SCARA MECHANISM')
        plt.xlabel('x (arbitrary units)')
        plt.ylabel('y (arbitrary units)')
        plt.axis((-2, 15, -2, 15))
        plt.draw()
        plt.pause(0.0001)
        plt.clf()
        plt.show()


if __name__ == "__main__":
    main()
