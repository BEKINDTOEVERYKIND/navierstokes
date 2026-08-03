"""
Attempt part 3: wake floor + joint exponent window.

Exponent conventions (all as powers of N; polylog factors G = (log Lambda)^2
tracked separately):
  parent frequency N, minor radius r_N = N^-1, strain g = N^beta
  carrier Lambda = N^{beta/2}, child Q = N^b, eps=b-1, Delta=alpha-beta, mu=beta-2b
  child block: minor r_Q = N^-b, major R_Q = N^{b(4-2alpha)}
  m12 := b(2 alpha - 4) - 1  (margin of nested-torus condition (12); R_Q = N^{-(1+m12)})
Wake floor: any transverse collar of length L inside the strain-coherent flank
(L <= r_N/G^2) has DC wake coefficient >= psi_inf / L, so
  wake/child >= R_Q/(r_N/G^2) = G^2 N^{-m12}.
Conditions:
  (1)   b>1, 2b<beta<alpha<5/2
  (12)  m12 > 0
  (17a) 5-2alpha > eps*Delta        (17b) mu > 2 eps*Delta
  (L2)  m12 > 2 b (b-1) Delta       [wake vs rho_Q^2 malignant tolerance]
  (L3)  beta/2 > 1 + m12            [Lambda R_Q > 1: carrier oscillates in core]
  (L2') m12 > (b-1)(2 alpha - beta - 1)  [paranoid: wake strain vs THIS parent at rho_N^2]
Claims to verify numerically:
  A. (L2') infeasible on the whole window (sup margin < 0)  -> the fork is real
  B. (1)+(12)+(17)+(L2)+(L3) nonempty; find max-min-margin point
  C. finite-N crossover N* where N^{m_min} >= G^2 = ((beta/2) ln N)^2 ... use G = (ln Lambda)^2
  D. DC-floor demo: oscillatory collars cannot reduce the zero-mode of d(chi)/dtau
"""
import numpy as np

rng = np.random.default_rng(42)

def margins(al, bt, bv):
    e, D, m = bv-1, al-bt, bt-2*bv
    m12 = bv*(2*al-4)-1
    return {
        'b>1'   : e,
        'beta>2b': m,
        'alpha>beta': D,
        '5/2>alpha': 2.5-al,
        '(12)'  : m12,
        '(17a)' : (5-2*al)-e*D,
        '(17b)' : m-2*e*D,
        '(L2)'  : m12-2*bv*e*D,
        '(L3)'  : bt/2-1-m12,
    }, m12, e, D

# ---- A. paranoid branch infeasibility
supA = -9
for _ in range(6_000_000):
    bv=rng.uniform(1.0,1.25); bt=rng.uniform(2*bv,2.5); al=rng.uniform(bt,2.5)
    if not (2*bv<bt<al<2.5): continue
    m12=bv*(2*al-4)-1
    if m12<=0: continue
    val = m12-(bv-1)*(2*al-bt-1)
    supA=max(supA,val)
print(f"A. paranoid (L2') sup margin over (1)+(12): {supA:+.5f}  (must be < 0)")
# analytic bound: m12 <= b-1 (alpha<=5/2); (b-1)(2a-bt-1) > (b-1)(a-1) > b-1 since a>2
print("   analytic: m12 <= b-1 < (b-1)(alpha-1) <= (b-1)(2alpha-beta-1)  [alpha>2, beta<alpha]")

# ---- B. joint window max-min margin
best=None
for _ in range(12_000_000):
    bv=rng.uniform(1.0,1.25); bt=rng.uniform(2*bv,2.5); al=rng.uniform(bt,2.5)
    mg,_,_,_ = margins(al,bt,bv)
    mm=min(mg.values())
    if mm<=0: continue
    if best is None or mm>best[0]: best=(mm,al,bt,bv)
if best:
    mm,al,bt,bv=best
    mg,m12,e,D = margins(al,bt,bv)
    print(f"B. joint window NONEMPTY. max-min-margin = {mm:.4f} at alpha={al:.4f}, beta={bt:.4f}, b={bv:.4f}")
    for k,v in mg.items(): print(f"     {k:11s} {v:+.4f}")
else:
    print("B. joint window EMPTY")

# candidates from the review
for (al,bt,bv,tag) in [(2.45,2.35,1.15,"cand A"),(2.49,2.45,1.21,"cand B")]:
    mg,_,_,_=margins(al,bt,bv)
    mm=min(mg.values())
    print(f"   {tag}: min-margin {mm:+.4f} " + " ".join(f"{k}:{v:+.3f}" for k,v in mg.items() if v==mm or v<0.05))

# ---- C. finite-N crossover for the best point
if best:
    mm,al,bt,bv=best
    m_min=mm
    import math
    # need N^{m_min} >= (ln Lambda)^4 with Lambda=N^{beta/2}; G=(ln Lambda)^2, floor factor G^2=(ln L)^4
    f=lambda lnN: m_min*lnN - 4*math.log(max((bt/2)*lnN,1.000001))
    lo,hi=1.0,1e9
    for _ in range(200):
        mid=(lo+hi)/2
        if f(mid)>0: hi=mid
        else: lo=mid
    print(f"C. crossover: N* = e^{hi:.1f} = 10^{hi/math.log(10):.1f}  (wake floor beats polylog above this)")

# ---- D. DC floor demo
print("D. DC floor: collar profiles chi from 1 to 0 over length L; zero-mode of dchi/dtau")
L=1.0; n=4096; tau=np.linspace(0,L,n)
for name,chi in [
    ("smooth ramp", 0.5*(1+np.cos(np.pi*tau/L))),
    ("oscillatory ramp", 0.5*(1+np.cos(np.pi*tau/L)) + 0.02*np.sin(200*np.pi*tau/L)*np.sin(np.pi*tau/L)),
    ("staircase", np.clip(1-np.floor(8*tau/L)/8,0,1)),
]:
    dchi=np.gradient(chi,tau)
    print(f"   {name:18s}: integral of dchi/dtau = {np.trapezoid(dchi,tau):+.4f}  (forced = -1)")
