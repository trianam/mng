import numpy as np
import scipy as sp
import scipy.misc
import sympy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

def deCasteljau(V, t, n = None):
    """
    calculate the value of the Bezier curve in a point.

    Arguments:
    - V (numpy.Array): the vector of the n+1 control points (matrix n+1 X 2: column 1 = x; column 2 = y);
    - t (number): the parameter on wich I want to evaluate the curve.
    - n (number): the grade of the curve;

    Return: a tuple containing
    - the point of the curve in t;
    - the control points of the splitted curve from 0 to t;
    - the control points of the splitted curve from t to 1;
    """

    if n is None:
        n = len(V)-1
    
    C2 = V.astype(float)
    C1 = np.empty_like(C2)
    for k in range(1, n+1):
        C1[k-1] = C2[0] 
        for i in range(0, n-k+1):
            C2[i] = (1.0 - t)*C2[i] + t*C2[i+1]

    C1[n] = C2[0]

    return (C2[0], C1, C2[::-1])

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

def bezier(V, tMin = 0., tMax = 1., step = 0.01):
    """
    calculate points of the Bezier curve.

    Arguments:
    - V (numpy.array): the vector of control vertex (matrix n X 2: column 1 = x; column 2 = y);
    - tMin (number): the min t of the evaluation of the curve;
    - tMax (number): the max t of the evaluation of the curve;
    - step (number): the granularity of the points for t in [0,1].

    Return: all the points of the curve, the number of them depends of the step.
    """
    C = [];
    for t in np.arange(tMin,tMax+step,step):
        C.append(deCasteljau(V, t)[0])

    C = np.array(C)

    return C

def splitControlPoints(V, tSplit):
    """
    calcutate points of the two Bezier curves obtained splitting
    the curve defined by the passed control vertex in tSplit

    Arguments:
    - V
    - tSplit

    Return:
    """

    return deCasteljau(V, tSplit)[1:3]


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

def exer2():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    plt.figure(figsize=(13, 8));
    if (isTridimensional(V)):
        plt.gca(projection='3d')
    drawVertexes(bezier(V))
    drawVertexes(V)
    plt.show()

def exer3():
    X = sympy.Poly('1+t+t**2')
    Y = sympy.Poly('t**3')
    P = mergePoly(X,Y)
    #P = np.array([[1,0],[1,0],[1,0],[0,1]])
    plt.figure(figsize=(13, 8));
    V = exp2bernstein(P)
    drawVertexes(bezier(V))
    drawVertexes(V)
    plt.show()

def exer4():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    plt.figure(figsize=(13, 8));
    if (isTridimensional(V)):
        plt.gca(projection='3d')
    C1,C2 = splitControlPoints(V, 0.6)
    drawVertexes(bezier(C1))
    drawVertexes(C1)
    drawVertexes(bezier(C2))
    drawVertexes(C2)
    plt.show()

def exer5():
    V1 = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    V2 = np.array([[0,0,0],[1,2,0],[3,2,3],[3,2,3],[2,6,2]])
    V3 = np.array([[0,0,0],[1,2,0],[3,2,3],[3,2,3],[3,2,3],[2,6,2]])
    V4 = np.array([[0,0,0],[1,2,0],[3,2,3],[3,2,3],[3,2,3],[3,2,3],[2,6,2]])
    plt.figure(figsize=(13, 8));
    if (isTridimensional(V1)):
        plt.gca(projection='3d')
    drawVertexes(V1)
    drawVertexes(bezier(V1))
    drawVertexes(bezier(V2))
    drawVertexes(bezier(V3))
    drawVertexes(bezier(V4))
    plt.show()


menu = {
    2 : exer2,
    3 : exer3,
    4 : exer4,
    5 : exer5,
    }
    
exer = 0
while((exer < 2) or (exer > 5)):
    exer = input('Execute exercise number: ')

menu[exer]()
