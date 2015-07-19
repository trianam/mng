#!/usr/bin/env python

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

import functions as fn
import sympy
import numpy as np
import math

def exer1_2():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    fn.drawVertexes(fn.bezier(V))
    fn.drawVertexes(V)

def exer1_3():
    X = sympy.Poly('1+t+t**2')
    Y = sympy.Poly('t**3')
    P = fn.mergePoly(X,Y)
    #P = np.array([[1,0],[1,0],[1,0],[0,1]])
    V = fn.exp2bernstein(P)
    fn.drawVertexes(fn.bezier(V))
    fn.drawVertexes(V)

def exer1_4():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    C1,C2 = fn.splitControlPoints(V, 0.6)
    fn.drawVertexes(fn.bezier(C1))
    fn.drawVertexes(C1)
    fn.drawVertexes(fn.bezier(C2))
    fn.drawVertexes(C2)

def exer1_5():
    V1 = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    V2 = np.array([[0,0,0],[1,2,0],[3,2,3],[3,2,3],[2,6,2]])
    V3 = np.array([[0,0,0],[1,2,0],[3,2,3],[3,2,3],[3,2,3],[2,6,2]])
    V4 = np.array([[0,0,0],[1,2,0],[3,2,3],[3,2,3],[3,2,3],[3,2,3],[2,6,2]])
    fn.drawVertexes(V1)
    fn.drawVertexes(fn.bezier(V1))
    fn.drawVertexes(fn.bezier(V2))
    fn.drawVertexes(fn.bezier(V3))
    fn.drawVertexes(fn.bezier(V4))

def exer1_6():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    V1 = fn.increaseGrade(V)
    V2 = fn.increaseGrade(V1)
    V3 = fn.increaseGrade(V2)
    fn.drawVertexes(fn.bezier(V))
    fn.drawVertexes(V)
    fn.drawVertexes(V1)
    fn.drawVertexes(V2)
    fn.drawVertexes(V3)

def exer1_7():
    V = np.array([[0,0,0],[1,2,0],[3,2,0],[6,-1,0]])
    Va = fn.attachWithC(V, np.array([[2,-1,3]]),C=2,hOrig=1.,hAttach=1.)
    fn.drawVertexes(fn.bezier(V))
    fn.drawVertexes(V)
    fn.drawVertexes(fn.bezier(Va))
    fn.drawVertexes(Va)

def exer2_1():
    V = np.array([[0,0],[1,2],[3,1]])
    W1 = np.array([1,1,1])
    W2 = np.array([1,2,1])

    fn.drawVertexes(V)
    P1 = fn.bezierRational(V,W1)[0]
    P2 = fn.bezierRational(V,W2)[0]
    fn.drawVertexes(P1)
    fn.drawVertexes(P2)

def exer2_2():
    V = np.array([[0,0,0],[1,2,0],[3,2,3],[2,6,2]])
    W = np.array([1,1,1,1])
    V1,W1 = fn.increaseGradeRational(V, W)
    V2,W2 = fn.increaseGradeRational(V1, W1)
    V3,W3 = fn.increaseGradeRational(V2, W2)
    fn.drawVertexes(fn.bezierRational(V,W)[0])
    fn.drawVertexes(V)
    fn.drawVertexes(V1)
    fn.drawVertexes(V2)
    fn.drawVertexes(V3)

def exer2_3():
    V = np.array([[1,0],[1,1],[0,1]])
    W = np.array([1,math.sqrt(2)/2,1])

    fn.drawVertexes(V)
    fn.drawVertexes(fn.bezierRational(V,W)[0])

def exer2_4():
    for a in [0.7, 1., 2., 3.]:
        V = np.array([[1,0],[a,a],[0,1]])
        W = np.array([1,fn.findArcWeigth(a),1])
        centre, radius = fn.findCircle(a)

        fn.drawVertexes(fn.circum(centre, radius))
        fn.drawVertexes(V)
        fn.drawVertexes(fn.bezierRational(V,W)[0])
        
    

def test1():
    V = np.array([[0.5,0.5], [0.4,1.], [1.,1.2], [0.9,0.], [2., 0.1], [1.5, 0.6]])
    k = 2
    p = np.array([0., 0., 1./5., 2./5., 3./5., 4./5., 1., 1.])
    #k = 3
    #p = np.array([0., 0., 0., 1./4., 1./2., 3./4., 1., 1., 1.])

    #print('pippo1 ' + str(fn.bSplineBaseFactory(p,k)(5,1)(0.9)))
    #print('pippo2 ' + str(fn.bSplineBaseFactory(p,k)(5,1)(1)))
    #return


    for i in range(6):
        N = fn.bSplineBaseFactory(p,k)(i)
        P = np.empty((0,2), dtype=float);
        for t in np.arange(0.,1.01,0.01):
            P = np.vstack((P, np.array([t,N(t)]))) 
  
        fn.drawVertexes(P)


def example1():
    p = np.array([-2., -1., 0., 1., 2., 3., 4., 5., 6., 7., 8.])
    k = 3

    for base in fn.bSplineBases(p, k):
        fn.drawVertexes(base)
    

def example2():
    p = np.array([0., 0., 0., 1., 2., 3., 4., 4., 5., 5., 5.])
    k = 3

    for base in fn.bSplineBases(p, k):
        fn.drawVertexes(base)


def example3():
    V = np.array([[-0.,0.], [-0.,6.], [-1.,5.], [-3.,8.], [-1.,14.],[2.,14.],[4.,8.], [2.,5.], [1.,6.], [1.,0.]])
    k = 4
    pc = np.array([0., 0., 0., 0., 1./7., 2./7., 3./7., 4./7., 5./7., 6./7., 1., 1., 1., 1.])
    pnc = np.array([-3./7., -2./7., -1./7., 0., 1./7., 2./7., 3./7., 4./7., 5./7., 6./7., 1., 8./7., 9./7., 10./7.])

    
    for base in fn.bSplineBases(pc, k):
        fn.drawVertexes(base)
    
    fn.drawVertexes(V)
    fn.drawVertexes(fn.bSpline(V, k, pc))

    for base in fn.bSplineBases(pnc, k):
        fn.drawVertexes(base)

    fn.drawVertexes(V)
    fn.drawVertexes(fn.bSpline(V, k, pnc))


def example4():
    V = np.array([[0.,0.], [0.,6.], [-1.,5.], [-3.,8.], [-1.,14.],[2.,14.],[4.,8.], [2.,5.], [1.,6.], [1.,0.]])
    k = 4


    fn.drawVertexes(V, 'o')
    fn.drawVertexes(V)
    fn.drawVertexes(fn.bSplineClosed(V, k))



def exerB_1():
    V = np.array([[0.5,0.5], [0.4,1.], [1.,1.2], [0.9,0.], [2., 0.1], [1.5, 0.6]])

    ks = [2, 3, 4, 5, 6]
    ps = [
        np.array([0., 0., 1./5., 2./5., 3./5., 4./5., 1., 1.]),
        np.array([0., 0., 0., 1./4., 1./2., 3./4., 1., 1., 1.]),
        np.array([0., 0., 0., 0., 1./3., 2./3., 1., 1., 1., 1.]),
        np.array([0., 0., 0., 0., 0., 1./2., 1., 1., 1., 1., 1.]),
        np.array([0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1.]),
    ]


    for p,k in zip(ps, ks):
        fn.drawVertexes(fn.bSpline(V, k, p))

    fn.drawVertexes(V, 'o')

def exerB_2():
    V = np.array([[0.,0.], [-0.4,-0.1], [-1.,0.2], [-0.1,1.], [1., 1.1], [1.1, 0.6], [1.2, 0.7], [1.3, 1.2], [2., 0.8], [2.5, 0.7]])

    ks = [4, 4, 6, 6, 8]
    ps = [
        np.array([0., 0., 0., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]),
        np.array([0., 0., 0., 0., 1., 1., 2., 2., 3., 3., 4., 4., 4., 4.]),
        np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 4., 5., 5., 5., 5., 5., 5.]),
        np.array([0., 0., 0., 0., 0., 0., 1., 2., 3., 3., 4., 4., 4., 4., 4., 4.]),
        np.array([0., 0., 0., 0., 0., 0., 0., 0., 1., 2., 3., 3., 3., 3., 3., 3., 3., 3.]),
    ]


    for p,k in zip(ps, ks):
        fn.drawVertexes(fn.bSpline(V, k, p))

    fn.drawVertexes(V, 'o')
    fn.drawVertexes(V)
    
def exerB_3():
    V = np.array([[0.,2.], [1.,0.], [2.,1.], [2.,1.], [3., 0.], [4., 2.]])

    k = 3
    p= np.array([0., 0., 0., 1./4., 1./2., 3./4., 1., 1., 1.])


    fn.drawVertexes(V, 'o')
    fn.drawVertexes(V)
    fn.drawVertexes(fn.bSpline(V, k, p))


def exerB_4():
    V = np.array([[1.,0.], [0.,1.], [2.,1.5], [2.,1.5], [4., 1.], [3., 0.]])

    k = 4
    ps = [
        np.array([0., 0., 0., 0., 1./4., 3./4., 1., 1., 1., 1.]),
        np.array([0., 0., 0., 0., 1./2., 1./2., 1., 1., 1., 1.]),
    ]


    fn.drawVertexes(V, 'o')
    fn.drawVertexes(V)

    for p in ps:
        fn.drawVertexes(fn.bSpline(V, k, p))


def exerB_5():
    Vs = [
        np.array([[0.,1.], [1.,0.], [2.,1.], [3., 0.], [4., 1.]]),
        np.array([[0.,1.], [1.,0.], [2.,1.], [2.,1.], [3., 0.], [4., 1.]]),
        np.array([[0.,1.], [1.,0.], [2.,1.], [2.,1.], [2.,1.], [3., 0.], [4., 1.]])
    ]

    ks = [4,4,4]
    ps = [
        np.array([0., 0., 0., 0.,1./2., 1., 1., 1., 1.]),
        np.array([0., 0., 0., 0., 1./4., 3./4.,1., 1., 1., 1.]),
        np.array([0., 0., 0., 0., 1./2., 1./2., 1./2.,1., 1., 1., 1.])
    ]
    

    fn.drawVertexes(Vs[0], 'o')
    fn.drawVertexes(Vs[0])
    for V,k,p in zip(Vs,ks,ps):
        fn.drawVertexes(fn.bSpline(V, k, p))


def exerB_6a():
    V = np.array([[1.,0.], [0.,1.], [2.,2.], [3., 0.]])
    k = 4


    fn.drawVertexes(V, 'o')
    fn.drawVertexes(V)
    fn.drawVertexes(fn.bSplineClosed(V, k))


def exerB_6b():
    V = np.array([[0.,0.], [0.,3.], [1.,2.], [2., 3.], [2.,0.], [1.,1.], [0.,0.]])
    k = 4


    fn.drawVertexes(V, 'o')
    fn.drawVertexes(V)
    fn.drawVertexes(fn.bSplineClosed(V, k))

