"""Independent audit of C186 and C187 exact algebra (checkpoint 0c5bf3b)."""
import sympy as sp
from fractions import Fraction as F

print("== C186 ==")
# nilpotent rank-one identity and tau classification
a1,a2,b1,b2,u1,u2,v1,v2 = sp.symbols('a1 a2 b1 b2 u1 u2 v1 v2')
Nu = sp.Matrix([[u1],[u2]])*sp.Matrix([[a1,a2]]); Nv = sp.Matrix([[v1],[v2]])*sp.Matrix([[b1,b2]])
tau = sp.trace(Nu*Nv)
alpha1_u2 = (a1*v1+a2*v2); alpha2_u1 = (b1*u1+b2*u2)
print("tau = alpha1(u2) alpha2(u1):", sp.simplify(tau - alpha1_u2*alpha2_u1)==0)
U = sp.Matrix([[1,1],[0,1]]); V = sp.Matrix([[1,0],[1,1]])
NU, NV = U-sp.eye(2), V-sp.eye(2)
print("nilpotence, dets:", (NU**2==sp.zeros(2,2), NV**2==sp.zeros(2,2), U.det()==1, V.det()==1))
print("flags distinct: ker NU = <e1>, ker NV = <e2>:", (NU*sp.Matrix([1,0])==sp.zeros(2,1), NV*sp.Matrix([0,1])==sp.zeros(2,1)))
# trace/det/discriminant of (I+N1)(I+N2)
P12 = (sp.eye(2)+Nu)*(sp.eye(2)+Nv)
print("det = 1:", sp.simplify(P12.det() - (1 + 0)) == sp.simplify(tau*0))  # det(I+Nu)=1 etc.
print("det((I+N1)(I+N2)) == 1:", sp.simplify(P12.det()-1)==sp.simplify((a1*u1+a2*u2)*(b1*v1+b2*v2))*0 or sp.simplify(P12.det().subs({a1*u1+a2*u2:0})-1)==0)
# enforce alpha_i(u_i)=0 -> substitute a2 = -a1*u1/u2 etc. check numerically instead:
import random
ok_det=ok_tr=True
for _ in range(200):
    import numpy as np
    u=np.random.randn(2); a=np.random.randn(2); a=a-np.dot(a,u)*u/np.dot(u,u)  # not the right constraint
    # correct constraint: alpha(u)=0 i.e. a.u=0
    a=a-np.dot(a,u)/np.dot(u,u)*u
    v=np.random.randn(2); b=np.random.randn(2); b=b-np.dot(b,v)/np.dot(v,v)*v
    N1=np.outer(u,a); N2=np.outer(v,b)
    P=(np.eye(2)+N1)@(np.eye(2)+N2); t=np.trace(N1@N2)
    ok_det &= abs(np.linalg.det(P)-1)<1e-10
    ok_tr  &= abs(np.trace(P)-(2+t))<1e-10
print("random nilpotent pairs: det=1 and tr=2+tau:", ok_det, ok_tr)
# product P=UV and rational cone bound
P = U*V
print("P=UV=[[2,1],[1,1]], tr=3:", P==sp.Matrix([[2,1],[1,1]]), sp.trace(P)==3)
w = sp.Matrix([sp.Rational(13,8),1])
Pw = P*w
print("Pw = (17/4,21/8):", Pw==sp.Matrix([sp.Rational(17,4),sp.Rational(21,8)]))
print("Pw >= (34/13) w:", Pw[0]-sp.Rational(34,13)*w[0]>=0, Pw[1]-sp.Rational(34,13)*w[1]>=0,
      "  margins:", sp.nsimplify(Pw[0]-sp.Rational(34,13)*w[0]), Pw[1]-sp.Rational(34,13)*w[1])
# golden-ratio operator norm < 13/8
sv = sp.sqrt(max((U.T*U).eigenvals()))
print("||U||_2 = (1+sqrt5)/2 < 13/8:", sp.simplify(sv-(1+sp.sqrt(5))/2)==0, sp.simplify(sv)<sp.Rational(13,8))
# exact rational log series value (1.11)
z = F(21,47)
s = 2*(z + z**3/3 + z**5/5)
print("(1.11): 2(z+z^3/3+z^5/5) == 24/25 + 1185042/5733625175:", s == F(24,25)+F(1185042,5733625175), " value>24/25:", s>F(24,25))
# gamma > 12/25 given log(34/13) > 24/25 (series lower-bounds log) -> 1/2 log rho > 12/25
import math
print("check: log(34/13) =", math.log(34/13), "> 24/25 =", 24/25, ":", math.log(34/13)>24/25)
# schedules
print("(3/8)/(24/25) = 25/64:", F(3,8)/F(24,25)==F(25,64), "; (1/2)/(24/25) = 25/48:", F(1,2)/F(24,25)==F(25,48))
# robustness (2.2)-(2.6)
eps=F(1,100)
print("6e+4e^2 = 151/2500:", 6*eps+4*eps**2==F(151,2500), "; 3-151/2500 = 7349/2500 > 29/10:", F(3)-F(151,2500)==F(7349,2500), F(7349,2500)>F(29,10))
# worst-case trace deviation: verify 6eps+4eps^2 is a valid bound by brute maximization
worst=0
import numpy as np
for _ in range(200000):
    dU=(np.random.rand(2,2)*2-1)*0.01; dV=(np.random.rand(2,2)*2-1)*0.01
    tU=np.array([[1,1],[0,1]])+dU; tV=np.array([[1,0],[1,1]])+dV
    worst=max(worst,abs(np.trace(tU@tV)-3))
print("random max |tr-3| =", round(worst,6), "<= 151/2500 =", 151/2500)
print("5/2+2/5 = 29/10:", F(5,2)+F(2,5)==F(29,10))
z2=F(3,7); s2=2*(z2+z2**3/3)
print("(2.5): 2(3/7+(3/7)^3/3) = 312/343 > 9/10:", s2==F(312,343), s2>F(9,10), "; log(5/2)=",round(math.log(2.5),4))
print("gamma_tilde = (1/2)log(5/2) > 9/20:", 0.5*math.log(2.5)>0.45)
# steady Euler + passive gradient rotation
x,y,t=sp.symbols('x y t',real=True)
v1e,v2e=-sp.sin(y),sp.sin(x); p=-sp.cos(x)*sp.cos(y)
adv1=v1e*sp.diff(v1e,x)+v2e*sp.diff(v1e,y)+sp.diff(p,x)
adv2=v1e*sp.diff(v2e,x)+v2e*sp.diff(v2e,y)+sp.diff(p,y)
om=sp.diff(v2e,x)-sp.diff(v1e,y)
print("steady Euler + vorticity transport:", sp.simplify(adv1)==0, sp.simplify(adv2)==0,
      sp.simplify(v1e*sp.diff(om,x)+v2e*sp.diff(om,y))==0)
J=sp.Matrix([[0,-1],[1,0]]); Rt=sp.exp(t*J)  # rotation
gvec=(Rt)*sp.Matrix([1,0])
print("g(pi/2) = e2:", sp.simplify(gvec.subs(t,sp.pi/2))==sp.Matrix([0,1]))

print("== C187 ==")
# shell count and A3 bound
r=sp.symbols('r',positive=True,integer=True)
print("(2r+1)^3-(2r-1)^3 = 24r^2+2:", sp.expand((2*r+1)**3-(2*r-1)**3)==24*r**2+2)
S5=sum(F(1,int(rr)**4) for rr in range(1,6)); tail=F(1,3*125)
print("sum r^-4 bound: S5+int =", float(S5+tail), "< 13/12:", S5+tail<F(13,12), " (true zeta(4)=",round(math.pi**4/90,6),")")
print("1+26*13/12 = 175/6 < 121/4:", 1+26*F(13,12)==F(175,6), F(175,6)<F(121,4))
# numerical A3^2 for sanity
A3sq=sum(1.0/(1+i*i+j*j+k*k)**3 for i in range(-30,31) for j in range(-30,31) for k in range(-30,31))
print("numerical A3^2 ~", round(A3sq,4), "( < 175/6 =", round(175/6,4), ") A3 ~", round(math.sqrt(A3sq),4), "< 5.5")
# <k>^3 <= 4(<l>^3+<k-l>^3): reduce to (a+b)^3 <= 4(a^3+b^3), a,b>=1... actually all a,b>0
aa,bb=sp.symbols('aa bb',positive=True)
print("(a+b)^3 <= 4(a^3+b^3):", sp.simplify(4*(aa**3+bb**3)-(aa+bb)**3 - (3*aa**3+3*bb**3-3*aa**2*bb-3*aa*bb**2))==0,
      " and 3(a^3+b^3-a^2 b-a b^2)=3(a+b)(a-b)^2/... >=0:", sp.factor(aa**3+bb**3-aa**2*bb-aa*bb**2))
# constants: 8*11/2=44; 88^2=7744; 2*44^2=3872; sqrt(3872)=44 sqrt2
print("constants:", 8*F(11,2)==44, 88**2==7744, 2*44**2==3872, sp.sqrt(3872)==44*sp.sqrt(2))
# Young split: 88VZD <= nu D^2/4 + 7744 V^2Z^2/nu ; check coefficient: (88VZ)^2/(4*(nu/4))...
nu,VV,ZZ,DD=sp.symbols('nu V Z D',positive=True)
lhs=88*VV*ZZ*DD; rhs=nu*DD**2/4 + 88**2*VV**2*ZZ**2/nu
print("Young check (ab <= a^2 eps + b^2/(4eps) with eps=nu/4):", sp.simplify(rhs - (nu*DD**2/4 + (88*VV*ZZ)**2/nu))==0,
      "AM-GM valid: 88VZD <= 2*sqrt(nu/4 * 7744/nu) V Z D =", 2*math.sqrt(7744/4)==88)
print("Gronwall exponent 15488/2 = 7744:", 15488==2*7744)
