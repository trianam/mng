import numpy as np
import scipy as sp
import scipy.misc
import sympy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def deCasteljau(V, n, t):
    """
    calculate the value of the Bezier curve in a point.

    Arguments:
    - V (numpy.Array): the vector of the n+1 control points (matrix n+1 X 2: column 1 = x; column 2 = y);
    - n (number): the grade of the curve;
    - t (number): the parameter on wich I want to evaluate the curve.

    Return: the point of the curve in t.
    """
    
    #Q = np.copy(V)
    Q = V.astype(float)
    for k in range(1, n+1):
        for i in range(0, n-k+1):
            Q[i] = (1.0 - t)*Q[i] + t*Q[i+1]
    return Q[0]

def isBidimensional(V):
    """
    evaluate length of points for checking if they are bidimensional
    """
    return len(V[0,:]) == 2

def isTridimensional(V):
    """
    evaluate length of points for checking if they are tridimensional
    """
    return len(V[0,:]) == 3

def bezier(V, step):
    """
    calculate points of the Bezier curve.

    Arguments:
    - V (numpy.array): the vector of control vertex (matrix n X 2: column 1 = x; column 2 = y);
    - step (number): the granularity of the points for t in [0,1].

    Return: all the points of the curve, the number of them depends of the step.
    """
    n = len(V)-1

    C = [];
    for t in np.arange(0,1+step,step):
        C.append(deCasteljau(V, n, t))

    C = np.array(C)

    return C

def drawVertexes(V):
    """
    draw the points.

    Arguments:
    - V (numpy.array): the vector of points to draw.
    """
    if (isBidimensional(V)):
        plt.plot(V[:,0], V[:,1])
    elif (isTridimensional(V)):
        plt.plot(V[:,0], V[:,1], V[:,2])
    else:
        print('ERR')

def exp2bernstein(P):
    """
    Transform a series of points given in the exponential base,
    in the equivalent series of points in the Bernstein base.

    Arguments:
    - P (numpy.array): the vector of points (matrix n X 2: column 1 = x; column 2 = y) in the exponential base.

    Return: the correspondent vector of points in the Bernstein base.
    """
    n = len(P)-1
    C = []
    for i in range(0,n+1):
        C.append(sum(P[j]*sp.misc.comb(n-j,i-j)/sp.misc.comb(n,i) for j in range(0,i+1)))
    C = np.array(C)
    return C    

def mergePoly(X, Y):
    """
    Merge two polynomies in one matrix of coefficients.
    The result is an array of points in the exponential base for
    the curve expressed by the parametric representation of the two
    polynomes.

    Arguments:
    - X (sympy.Poly): the first polynome;
    - Y (sympy.Poly): the second polynome;

    Return: the array of points in the exponential base. 
    """
    degX = X.degree()
    degY = Y.degree()
    deg = max(degX, degY)
    coeffX = np.zeros(deg+1)
    coeffY = np.zeros(deg+1)
    coeffX[-degX-1:] = X.all_coeffs()[-degX-1:]
    coeffY[-degY-1:] = Y.all_coeffs()[-degY-1:]
    return np.column_stack((coeffX[::-1], coeffY[::-1]))

V = np.array([[0,0,0],[1,2,0],[3,2,3]])
plt.figure(figsize=(13, 8));
if (isTridimensional(V)):
    plt.gca(projection='3d')
drawVertexes(bezier(V, 0.01))
drawVertexes(V)
plt.show()

X = sympy.Poly('1+t+t**2')
Y = sympy.Poly('t**3')
P = mergePoly(X,Y)
#P = np.array([[1,0],[1,0],[1,0],[0,1]])
plt.figure(figsize=(13, 8));
V = exp2bernstein(P)
drawVertexes(bezier(V, 0.01))
drawVertexes(V)
plt.show()

