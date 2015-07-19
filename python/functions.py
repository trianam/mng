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

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import scipy as sp
import scipy.misc
import math
import matplotlib.pyplot as plt

#Bezier

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

def rational2Poly(V, W):
    Vcomp = np.column_stack((V,np.ones(len(V))))
    Vcomp = np.multiply(Vcomp,W.reshape((len(W),1)))
    return Vcomp

def poly2Rational(P):
    return (np.divide(P[:,:-1], P[:,-1].reshape((len(P[:,-1]),1))), P[:,-1])

def bezierRational(V, W, tMin = 0., tMax = 1., step = 0.01):
    """
    calculate points of the rational Bezier curve.

    Arguments:
    - V (numpy.array): the vector of control vertex (matrix n X 2: column 1 = x; column 2 = y);
    - W (numpy.array): the vector of weigths (the dimension must be the same as the rows of V
    - tMin (number): the min t of the evaluation of the curve;
    - tMax (number): the max t of the evaluation of the curve;
    - step (number): the granularity of the points for t in [0,1].

    Return: all the points of the curve, the number of them depends of the step.
    """
    Vcomp = rational2Poly(V, W)
    P = bezier(Vcomp, tMin, tMax, step)
    return poly2Rational(P)
    
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


def drawVertexes(V, opt = ''):
    """
    draw the points.

    Arguments:
    - V (numpy.array): the vector of points to draw.
    - opt: the option string to pass to plot
    """
    if (isBidimensional(V)):
        plt.gca(projection='rectilinear')
        plt.plot(V[:,0], V[:,1], opt)
    elif (isTridimensional(V)):
        plt.gca(projection='3d')
        plt.plot(V[:,0], V[:,1], V[:,2], opt)
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
    Vi = np.zeros(tuple(map(lambda a,b: a+b, V.shape, (2,0))))
    Vi[1:-1] = V
    for i in range(n+2,0,-1):
        Vi[i] = ((i-1.)/(n+1.))*Vi[i-1] + (1. - (i-1.)/(n+1.))*Vi[i]

    return Vi[1:]

def increaseGradeRational(V, W):
    """
    Increase by one the grade of the rational curve represented by passed control vertexes and weigths.
    """
    return poly2Rational(increaseGrade(rational2Poly(V, W)))

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

    Vattach = np.zeros(tuple(map(lambda a,b: a+b, Vextra.shape, (C+1,0))))
    if start:
        Vattach[:-C-1] = Vextra
        Vattach[-1] = Vorig[0]
        if(C >= 1):
            Vattach[-2] = ((hAttach+hOrig)/hOrig)*Vorig[0] - (hAttach/hOrig)*Vorig[1]
            if(C >= 2):
                Vattach[-3] = (hAttach/hOrig)*Vattach[-2]+Vattach[-2]-(hAttach/hOrig)*Vorig[1]-((hAttach/hOrig)**2)*(Vorig[1]-Vorig[2])
    else:
        Vattach[C+1:] = Vextra
        Vattach[0] = Vorig[-1]
        if(C >= 1):
            Vattach[1] = ((hAttach+hOrig)/hOrig)*Vorig[-1] - (hAttach/hOrig)*Vorig[-2]
            if(C >= 2):
                Vattach[2] = (hAttach/hOrig)*Vattach[1]+Vattach[1]-(hAttach/hOrig)*Vorig[-2]-((hAttach/hOrig)**2)*(Vorig[-2]-Vorig[-3])

    return Vattach

def findArcWeigth(a):
    return math.sin(math.acos((-2.*a*(1.-a))/(a**2+(1.-a)**2))/2.)

def findCircle(a):
    centre = ([(1.-a)/(1.-2.*a), (1.-a)/(1.-2.*a)])
    radius = math.sqrt(((1.-a)**2 + a**2)/(1-2.*a)**2)
    return (centre, radius)

def circum(c, r, n=100):
    return np.array([(math.cos(2*math.pi/n*x)*r + c[0], math.sin(2*math.pi/n*x)*r + c[1]) for x in xrange(0,n+1)])


#B-splines

def bSplineBaseFactory(part,k):
    """
    Define the B-splines base for the given extended partition.

    Arguments:
    - part: the extended partition of the parametric domain.
    - k: the order

    Returns: the function that give a base for a given i and h.
    """
    if (part[len(part)-k+1 : len(part)] == part[len(part)-k : len(part)-1]).all():
        iEnd = len(part)-k
    else:
        iEnd = len(part) -1

    def N(i,h=k):
        """
        Define a single base of B-splines for a given i and h.
        
        Arguments:
        - i: the index of the base
        - h: the order of the base

        Returns:
        - the function in t of the base
        """
        def bSpline(t):
            if (h > 1):
                n1 = N(i,h-1)
                n2 = N(i+1, h-1)

                a1 = ((t - part[i]) * n1(t))
                a2 = (part[i+h-1] - part[i])
                if(a1 == 0.):
                    a = 0.
                else:
                    a = (a1 / a2)

                b1 = ((part[i+h] - t) * n2(t))
                b2 = (part[i+h] - part[i+1])
                if(b1 == 0.):
                    b = 0.
                else:
                    b = (b1 / b2)

                return a + b
            else:
                #if ((part[i] <= t) and (t < part[i+1])):
                #if ((part[i] <= t) and ((t < part[i+1]) or ((i+1 == len(part)-k) and (t == part[i+1])))):
                #if ((part[i] <= t) and ((t < part[i+1]) or ((i+1 == len(part)-1) and (t == part[i+1])))):
                if ((part[i] <= t) and ((t < part[i+1]) or ((i+1 == iEnd) and (t == part[iEnd])))):

                    return 1.
                else:
                    return 0.
        return bSpline
    return N

def bSplineFun(V, k, dom=(0.,1.), part=None, m=None):
    """
    Define a B-spline function using the given control vectors, of
    the given order, and using either the given extended partition
    or creating one uniform and clamped in the given domain and
    molteplicity.

    Arguments:
    - V: the control points vector;
    - k: the order;
    - dom: the parametric domain, mandatory if a partition is not
    given, but with default value (0,1);
    - part: the extended partition of the parametric domain,
    mandatory if the multiplicity and domain are not given;
    - m: the multiplicity, mandatory if the partition is not given.

    Returns: the function in t of the curve
    """
    if (part == None):
        part = np.concatenate((np.zeros(k)+dom[0], np.repeat(np.linspace(dom[0], dom[1], len(m)+2)[1:-1], m), np.zeros(k)+dom[1]))
    else:
        dom = (part[k-1], part[-k])

    N = bSplineBaseFactory(part,k)
    return lambda t: sum(V[i]*N(i)(t) for i in range(len(V)))

def bSpline(V, k, part, step = 0.01):
    """
    calculate the points of the B-spline curve defined by the
    given control points, the given extended partition, and of
    the given order.

    Arguments:
    - V: the control points vector;
    - k: the order;
    - part: the extended partition of the parametric domain;
    - step: the steps in the domain for creating the points.

    Returns:
    the array of points in the space of the curve.
    """
    C = bSplineFun(V,k,part=part)
    
    P = np.empty((0,V.shape[1]), dtype=float);
    for t in np.arange(part[k-1],part[-k]+step,step):
        P = np.append(P, [C(t)], axis=0)

    return P

def bSplineClosed(V, k, step = 0.01):
    """
    calculate the points of the closed B-spline curve defined by the
    given control points and of the given order.

    Arguments:
    - V: the control points vector, if the poligon is "open" that will be "closed";
    - k: the order;
    - step: the steps in the domain for creating the points.

    Returns:
    the array of points in the space of the curve.
    """
    if not (V[0] == V[-1]).all():
        V = np.vstack((V,V[0]))
    return bSpline(
        np.concatenate((V,V[1:k+1])),
        k,
        np.arange(float(-k)/float(len(V)-1), float(k+len(V))/float(len(V)-1), 1./float(len(V)-1)),
        step)

def bSplineBases(part, k, step = 0.01):
    """
    calculate the points of the B-spline bases defined by the
    given extended partition, and of
    the given order.

    Arguments:
    - part: the extended partition of the parametric domain;
    - k: the order;
    - step: the steps in the domain for creating the points.

    Returns:
    the array of array of points of the curves.
    """

    for i in range(len(part)-k):
        N = bSplineBaseFactory(part,k)(i)
        P = np.empty((0,2), dtype=float);
        #for t in np.arange(part[k-1],part[-k]+step,step):
        for t in np.arange(part[0],part[-1]+step,step):
            P = np.vstack((P, np.array([t,N(t)]))) 

        if i == 0:
            A = P
        else:
            A = np.dstack((A,P))

    A = np.swapaxes(A,0,2)
    A = np.swapaxes(A,1,2)
    return A
