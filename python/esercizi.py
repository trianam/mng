#!/usr/bin/python

# Copyright (C) 2015  Stefano Martina

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

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
    - V: the control points of the original curva
    - tSplit: the point in the parametric domain where to split the curve

    Return: a tuple with the control points of the first
            and the control points of the second part of the curve
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

def increaseGrade(V):
    """
    Increase by one the grade of the curve represented by passed control vertexes.
    """
    n = len(V) - 1
    W = np.zeros(tuple(map(lambda a,b: a+b, V.shape, (2,0))))
    W[1:-1] = V
    for i in range(n+2,0,-1):
        W[i] = ((i-1.)/(n+1.))*W[i-1] + (1. - (i-1.)/(n+1.))*W[i]

    return W[1:]

def attachWithC(Vorig, Vextra, C=2, hOrig=1., hAttach=1., start=True):
    """
    Attach to the curve defined by control vertexes Vorig
    another curve with continuity C and extra vertexes Vextra.
    To the start or end of Vorig. The dimension of returned will be
    the dimension of Vextra + C+1.

    Arguments:
    - Vorig: the control vertexes of the original curve, they are not changed;
    - Vextra: attach those vertex adding two further vertexes for keeping C continuity;
    - C: the continuity, can be 0, 1 or 2;
    - hOrig: the length of the parametrization domain of vOrig in the total curve;
    - hAttach: the length of the parametrization domain of the attached curve in the total curve;
    - start: if true attach to the start of Vorig, if false to the end.

    Returns: The control vertexes of the curve attached to Vorig
    """

    W = np.zeros(tuple(map(lambda a,b: a+b, Vextra.shape, (C+1,0))))
    if start:
        W[:-C-1] = Vextra
        W[-1] = Vorig[0]
        if(C >= 1):
            W[-2] = ((hAttach+hOrig)/hOrig)*Vorig[0] - (hAttach/hOrig)*Vorig[1]
            if(C >= 2):
                W[-3] = (hAttach/hOrig)*W[-2]+W[-2]-(hAttach/hOrig)*Vorig[1]-((hAttach/hOrig)**2)*(Vorig[1]-Vorig[2])
    else:
        W[C+1:] = Vextra
        W[0] = Vorig[-1]
        if(C >= 1):
            W[1] = ((hAttach+hOrig)/hOrig)*Vorig[-1] - (hAttach/hOrig)*Vorig[-2]
            if(C >= 2):
                W[2] = (hAttach/hOrig)*W[1]+W[1]-(hAttach/hOrig)*Vorig[-2]-((hAttach/hOrig)**2)*(Vorig[-2]-Vorig[-3])

    return W

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

def exer6():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    V1 = increaseGrade(V)
    V2 = increaseGrade(V1)
    V3 = increaseGrade(V2)
    plt.figure(figsize=(13, 8));
    if (isTridimensional(V)):
        plt.gca(projection='3d')
    drawVertexes(bezier(V))
    drawVertexes(V)
    drawVertexes(V1)
    drawVertexes(V2)
    drawVertexes(V3)
    plt.show()

def exer7():
    V = np.array([[0,0,0],[1,2,0],[3,2,0],[6,-1,0]])
    W = attachWithC(V, np.array([[2,-1,3]]),C=2,hOrig=1.,hAttach=1.)
    plt.figure(figsize=(13, 8));
    if (isTridimensional(V)):
        plt.gca(projection='3d')
    drawVertexes(bezier(V))
    drawVertexes(V)
    drawVertexes(bezier(W))
    drawVertexes(W)
    plt.show()
    
menu = {
    2 : exer2,
    3 : exer3,
    4 : exer4,
    5 : exer5,
    6 : exer6,
    7 : exer7,
    }


def checkIn(v):
    try:
        return int(v)
    except:
        return 0

exer = 0
while(exer not in menu.keys()):
    exer = checkIn(raw_input('Execute exercise number: '))

menu[exer]()
