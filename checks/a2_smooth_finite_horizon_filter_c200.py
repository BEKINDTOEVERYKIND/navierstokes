#!/usr/bin/env python3
"""Exact constant/algebra checker for C200's finite-horizon selector theorem.

Checks the complex graph-transform margins, both real spatial-derivative
ledgers, Cauchy/clock constants, and coefficient propagation/normalization
identities. It imports only the standard-library C194 checker to reconstruct
its third-flow-jet ledger. The off-ray cone certificate remains C195's.
No PDE discretization, interval profile search, or blow-up is certified.
"""
from fractions import Fraction as F
from math import factorial


def p_add(a, b):
    return [sum((p[i] if i < len(p) else F(0) for p in (a, b)), F(0))
            for i in range(max(len(a), len(b)))]


def p_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i+j] += ai*bj
    return out


def p_scale(a, c):
    return [x*c for x in a]


def p_integral(a):
    return [F(0)] + [x/F(i+1) for i, x in enumerate(a)]


def third_flow_jet():
    # Re-run the rational quotient sub-ledger rather than treating J3 as
    # an unsupported consequence of the J1/J2 bounds.
    from a2_curl_wkb_bridge_c194 import second_corrector_j3_ledger
    second_corrector_j3_ledger()
    w, v = [F(1), F(4)], [F(0), F(1), F(2)]
    j1 = p_scale(w, 18)
    j2 = p_add(p_add(p_scale(p_mul(w, w), 756), p_scale(v, 16524)),
               p_scale(w, 672))
    da = p_scale(v, 918)
    dda = p_integral(p_add(p_scale(p_mul(j1, j1), 9396), p_scale(j2, 51)))
    a2 = p_add(p_scale(p_mul(j1, j1), 341), p_scale(j2, 14))
    terms = [p_scale(p_mul(a2, w), 3), p_scale(dda, 24),
             p_scale(w, 63408), p_scale(p_mul(j1, da), 84),
             p_scale(p_mul(j1, w), 2352), p_scale(da, 1344)]
    out = [F(0)]
    for term in terms:
        out = p_add(out, term)
    out = p_scale(out, 9)
    assert out == list(map(F, [4474548, 749731896, 3071113056, 4031918208]))
    assert all(x <= 4031918208*y for x, y in zip(out, [1, 3, 3, 1]))


def complex_block_ledger():
    delta = F(1, 10**40)  # L factor explicitly cancels in eta.
    eta = 1296*delta
    assert eta < F(1, 10**6)
    assert 9 - 29*eta > 8
    projector_difference = F(29, 8) + F(196*29, 72)
    assert projector_difference < 100
    assert F(30, 8) + F(225*30, 64) < 110
    assert F(2, 8) + F(2*30**2, 64) + F(225*2, 64) + F(2*225*30**2, 512) < 1000
    assert 1200*eta < 1
    assert 6*F(76, 25) < 19 and 7*F(76, 25) < 22
    assert 3**19 < 2*10**9 and 3**22 < 10**11
    assert 8*3**19 < 10**11
    assert 15*10**11 < 2*10**12
    assert 4*10**11*1200 < 10**15
    frame_perturb = 100*10**11*5 + 2*10**15*5 + 2*(2*10**9)*100
    assert frame_perturb < 2*10**16
    assert 2*10**16*1296 < 10**20

    # First and second k-derivatives of analytic section frames.
    # E1=n-mk/Q: |mk|<=75, |D(mk)|<=20, |D²(mk)|<=2.
    assert F(2, 8) + 2*20*F(1, 2) + 75*4 < 10**5
    # E+ numerator Qn-mk: bounds 300,50,4; D^-1 jets 1/8,1/6,1/2.
    assert F(4, 8) + 2*50*F(1, 6) + 300*F(1, 2) + F(1, 3) + F(5, 2) < 10**5


def graph_margins():
    sigma, error = F(1, 10**14), F(1, 10**20)
    assert 3000 - 10**11*sigma - 2*error > 2999
    assert F(2, 2999*3000) < F(1, 10**6)
    quotient_error = F(25, 16)*error/2999
    assert quotient_error < F(1, 10**23)
    assert sigma/F(10**6) + F(1, 10**23) < sigma
    assert 2*10**11*error + error**2 < 1
    assert F(3, 2999**2) < F(1, 10**6)
    assert F(1, 10**23)/(1-F(1, 10**6)) < F(2, 10**23)
    assert 2*F(1, 10**23) < sigma/F(10**5)
    # q has a single analytic square-root branch near its real positive value.
    adiff = 100*1296*F(1, 10**40)*F(5, 4) + 5*F(2, 10**23)
    assert adiff < F(1, 10**16)
    assert F(1, 4) - 4*adiff-adiff**2 > F(1, 9)


def spatial_ledger():
    ky = 2*10**11
    assert 12*216**2*286992 < ky
    assert F(15, 2)*216 == 1620
    ly = 3*10**14
    assert 5*1620 + 12*110*ky < ly
    hy = 2*10**26
    assert 4*10**11*ly < hy
    ey = 100*ky
    assert ey <= 2*10**13
    my = 10**28
    assert ey*10**11*5 + 3*hy*5 + 3*10**11*ey < my
    sy = 2*10**28
    assert F(my, 1)/(1-F(1, 10**6)) < sy
    ay = 10**30
    assert ey*F(5, 4) + 5*sy < ay
    qy = 10**31
    assert 3*ay*3 < qy
    assert 3*ay+27*qy < 10**33
    gy = 10**41
    assert F(5, 4)*my + 2*10**12*sy < gy
    assert F(gy, 2999) < 10**38
    # R+1 >=1 absorbs the finite product count.
    assert 3*(ay+3*(3*qy+10**38)) < 10**45

    kyy = 10**21
    assert 12*(2*216**3*286992**2 + 216**2*4031918208) < kyy
    assert 9*216**2 + F(15, 2)*286992 < 3*10**6
    lyy = 10**28
    assert 5*(3*10**6)+4*110*ky*1620+12*1000*ky**2+12*110*kyy < lyy
    hyy = 2*10**41
    assert 10**11*(4*lyy+16*ly**2) < hyy
    eyy = 10**28
    assert 10**5*ky**2+100*kyy < eyy
    myy = 10**45
    matrix_terms = [eyy*10**11*5, 3*hyy*5, 3*10**11*eyy,
                    2*ey*hy*5, 2*ey*10**11*ey, 2*3*hy*ey]
    assert sum(matrix_terms) < myy
    fyy, fys, fss = 10**57, 10**42, 10**3
    assert F(2*myy+4*my**2, 2999) < fyy
    assert (my+F(1, 10**6)*2*my+F(1, 4)*my+my*2*10**12)/2999 < fys
    assert F(2*(2*10**12)*3, 2999**3) < fss
    inhom = fyy+2*fys*sy+fss*sy**2
    assert inhom < 10**71
    syy = 10**73
    assert F(inhom, 1)/(1-F(1, 10**6)) < syy
    ayy = 10**74
    assert F(5, 4)*eyy+2*ey*sy+5*syy < ayy
    qyy = 10**76
    assert 3*(ay**2+3*ayy+qy**2) < qyy
    gyy = 10**86
    assert F(5, 4)*myy+2*my*sy+2*10**12*syy < gyy
    assert F(gyy, 2999)+10**76 < 10**84
    assert 3*qyy+9*qy**2 < 10**77
    assert 3*qy+10**38 < 10**39
    assert 10**77+10**84 < 10**85
    assert 3*(ayy+2*ay*10**39+3*(10**78+10**85)) < 10**90


def mv(m, v):
    return [sum((a*b for a, b in zip(row, v)), F(0)) for row in m]


def det(m):
    return m[0][0]*m[1][1]-m[0][1]*m[1][0]


def reflected_inverse(m):
    d = det(m)
    return [[m[1][1]/d, m[0][1]/d], [m[1][0]/d, m[0][0]/d]]


def graph(m, s):
    return (m[1][0]+m[1][1]*s)/(m[0][0]+m[0][1]*s)


def gain(m, s):
    return m[0][0]+m[0][1]*s


def algebra_sentinels():
    # Exact rational, noncommuting blocks test both forward and future-defined
    # backward coefficients. This supplements the displayed general identities.
    matrices=[]
    for g,s in [(F(4000), F(4,25)), (F(5000), F(9,50)), (F(4500), F(17,100))]:
        a=(g+1/g)/2
        b=(g-1/g)/2
        m=[[a,b/s],[b*s,a]]
        assert det(m)==1
        matrices.append(m)
    s=[F(17,100)]
    gains=[]
    for m in matrices:
        gains.append(gain(m,s[-1]))
        s.append(graph(m,s[-1]))
    r=[F(0)]*3+[F(17,100)]
    inverse_gains=[F(0)]*3
    for j in range(2,-1,-1):
        mi=reflected_inverse(matrices[j])
        inverse_gains[j]=gain(mi,r[j+1])
        r[j]=graph(mi,r[j+1])
    for j,m in enumerate(matrices):
        assert mv(m,[F(1),s[j]]) == [gains[j],gains[j]*s[j+1]]
        assert mv(m,[F(1),-r[j]]) == [1/inverse_gains[j],-r[j+1]/inverse_gains[j]]
        for z in [F(137,1000),F(17,100),F(1,5)]:
            w=F(181,1000)
            assert graph(m,z)-graph(m,w)==det(m)*(z-w)/(gain(m,z)*gain(m,w))
    # The two true endpoint normalizations in (20) follow after multiplying
    # by each fixed q_R^+ or q_0^-; its real square is physical frame norm².
    for j in range(3):
        assert gains[j]>3000 and inverse_gains[j]>3000
    assert F(1,2)**2 < 1 and F(2)**2 > 1


def clock():
    # e=sum1/n!, with tail beginning1/7! and subsequent ratios <=1/8.
    e_upper=sum((F(1,factorial(j)) for j in range(7)),F(0))+F(1,factorial(7))/(1-F(1,8))
    assert e_upper < F(87,32)
    assert F(87,32)**8 < 2999
    assert 8*F(3,8)==3
    assert F(3,8)==F(3,1)/8
    # Cauchy ball-to-polydisk and Euclidean multilinear conversions.
    assert 3 < 4 and F(4)*10**40 > 1


def main():
    third_flow_jet()
    complex_block_ledger()
    graph_margins()
    spatial_ledger()
    algebra_sentinels()
    clock()
    print('C200 PASS: finite-horizon graph/normalization algebra; complex margins; covector and two spatial jet constants; clock.')


if __name__ == '__main__':
    main()
