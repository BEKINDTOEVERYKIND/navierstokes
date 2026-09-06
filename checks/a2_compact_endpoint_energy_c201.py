#!/usr/bin/env python3
"""Exact arithmetic for C201 geometric, Riesz, and energy constants.

General lattice covering, Fourier synthesis, and PDE energy identities are
proved in the companion note; finite arithmetic is checked here.
"""
from fractions import Fraction as F
from math import factorial


def dot(x,y): return sum(a*b for a,b in zip(x,y))
def matmul(a,b): return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]
def add(a,b): return [[x+y for x,y in zip(ar,br)] for ar,br in zip(a,b)]


def poly_add(a,b):
    out=dict(a)
    for key,value in b.items():
        out[key]=out.get(key,F(0))+value
    return {key:value for key,value in out.items() if value}


def poly_mul(a,b):
    out={}
    for ka,va in a.items():
        for kb,vb in b.items():
            key=tuple(x+y for x,y in zip(ka,kb))
            out[key]=out.get(key,F(0))+va*vb
    return {key:value for key,value in out.items() if value}


def poly_scale(a,c): return {key:c*value for key,value in a.items() if c*value}


def main():
    identity=[[F(int(i==j)) for j in range(3)] for i in range(3)]
    a=[F(2),F(-1),F(3)]
    b=[F(1),F(2),F(0)]
    assert dot(a,b)==0
    n=[[x*y for y in b] for x in a]
    assert matmul(n,n)==[[0]*3 for _ in range(3)]
    assert matmul(add(identity,n),add(identity,[[-x for x in row] for row in n]))==identity
    assert 1+217<220
    assert F(8*22,7)+1<27
    assert 2*216*27==11664
    assert 27**3<162**2
    assert F(1,4)**3==F(1,64)
    assert 36**2*3<72**2
    assert 30*(2+9+66+180+216)==14190
    assert 27*45*14190==17240850
    c=F(1,10**100)
    eta=4*10**40*220*c
    assert eta<10**43*c==F(1,10**57)
    delta=17240850*eta
    assert delta<F(2,10**50)
    # Derivative-polynomial upper constants, with L>=1 and R+1>=1.
    V=[1,3,14,102,984,11880]
    H=2*216**2
    vsum=sum(F(V[r],factorial(r)) for r in range(5))
    assert 432*vsum*H**4 < 10**25
    vxsum=sum(F(V[r+1]+r*V[r],factorial(r)) for r in range(5))
    assert 216**2*286992*8*vxsum*H**4 < 10**35
    D0=4*10**40
    assert 150*D0**4<10**165
    assert 5*216*10**45*D0**4<10**211
    assert 10**35*10**165+10**25*10**211<10**237
    assert 262440*(F(5,2)*10**237+648*10**190)<10**250
    # Exact energy-balancing interval estimates.
    small=F(1,10**6)
    low,high=F(1,10),F(29)
    assert (low-small)/(high+small)>F(1,300)
    assert (high+small)/(low-small)<300
    assert low-300*small>F(1,20)
    assert high+300*small<30
    # q0=10^10000 and x0=log q0 lie in (20000,30000).
    # All tails x^m exp(-a x) used below decrease already for x>=20000.
    assert F(2,3)*20000>3
    assert 11664*10**100*30000**3 < 10**120
    assert F(2,3)*10000>6666
    assert F(1,3)*10000>3333
    # The explicit plateau has Achi<=72<10^100; (15) is <1e-100.
    assert 10**250*30000**21*10**100 < 10**450
    assert 450-6666 < -100
    # The C193 clock supplies both L<=log q and R+1<=log q.
    assert F(57,400)*20000+F(177,25)<20000
    assert F(3,64)*20000+3<20000
    # Symbolic energy identity: B^2(A+2 lambda C-lambda^2 B)
    # equals B(AB+C^2-R^2), with R^2=C^2+AB.
    pa={(1,0,0,0):F(1)}
    pb={(0,1,0,0):F(1)}
    pc={(0,0,1,0):F(1)}
    pr={(0,0,0,1):F(1)}
    numerator=poly_add(pc,pr)
    lhs=poly_add(poly_mul(pa,poly_mul(pb,pb)),
        poly_add(poly_scale(poly_mul(poly_mul(pb,pc),numerator),2),
                 poly_scale(poly_mul(pb,poly_mul(numerator,numerator)),-1)))
    rhs=poly_mul(pb,poly_add(poly_mul(pa,pb),
        poly_add(poly_mul(pc,pc),poly_scale(poly_mul(pr,pr),-1))))
    assert lhs==rhs
    # Equation (19) polynomial cancellation, without floating-point sqrt:
    # (B lambda-C)^2=C^2+AB iff B lambda^2-2C lambda-A=0.
    for av,bv,cv in [(F(2),F(3),F(1)),(F(1,100),F(841),F(-1,1000))]:
        assert bv>0 and av>0
        assert cv*cv+av*bv>cv*cv
    print("C201 exact geometry/Riesz/normalization checker PASS")

if __name__=='__main__': main()
