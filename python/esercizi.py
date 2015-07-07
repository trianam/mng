import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def deCasteljau(V, n, t):
    """
    deCasteljau calculate the value of the Bezier curve in a point
    @type V: numpy.Array
    @param V: the vector of the n+1 control points (matrix n+1 X 2: column 1 = x; column 2 = y).
    @type n: number
    @param n: the grade of the curve.
    @type t: number
    @param t: the parameter on wich I want to evaluate the curve.
    """
    
    #Q = np.copy(V)
    Q = V.astype(float)
    for k in range(1, n):
        for i in range(0, n-k):
            Q[i] = (1.0 - t)*Q[i] + t*Q[i+1]
    return Q[0]

def isBidimensional(V):
    return len(V[0,:]) == 2

def isTridimensional(V):
    return len(V[0,:]) == 3

def drawBezier(V):
    """
    drawBezier draw the Bezier curve
    @type V: numpy.array
    $param V: the vector of control vertex (matrix n X 2: column 1 = x; column 2 = y)
    """
    n = len(V)

    C = [];
    for t in np.arange(0,1,0.01):
        C.append(deCasteljau(V, n, t))

    C = np.array(C)

    if (isBidimensional(C)):
        plt.plot(C[:,0], C[:,1])
    elif (isTridimensional(C)):
        plt.plot(C[:,0], C[:,1], C[:,2])
    else:
        print('ERR')

def drawControlVertexes(V):
    if (isBidimensional(V)):
        plt.plot(V[:,0], V[:,1])
    elif (isTridimensional(V)):
        plt.plot(V[:,0], V[:,1], V[:,2])
    else:
        print('ERR')


V = np.array([[0,0,0],[1,2,0],[3,2,3]])
plt.figure(figsize=(13, 8));
if (isTridimensional(V)):
    plt.gca(projection='3d')
drawBezier(V)
drawControlVertexes(V)
plt.show()

