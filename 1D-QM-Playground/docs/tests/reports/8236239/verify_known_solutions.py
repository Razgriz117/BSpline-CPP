"""Verify the four TISE test outputs against the known analytic solutions."""
import numpy as np, json, os
from pathlib import Path
from scipy.optimize import brentq
from scipy.special import eval_hermite, factorial, genlaguerre
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = Path("/tmp/BSpline-CPP/1D-QM-Playground/data")
OUT = Path("./verify_out"); OUT.mkdir(exist_ok=True)
report = {}

def load_eigs(t):
    a = np.loadtxt(D/t/"tise/eigenvalues.dat"); return a[:,1]
def load_state(t, j):
    a = np.loadtxt(D/t/f"tise/eigenstate_{j:03d}.dat"); return a[:,0], a[:,1]
def load_cont(t, j):
    a = np.loadtxt(D/t/f"tise/continuum_state_{j:03d}.dat"); return a[:,0], a[:,1]
def load_ps(t):
    return np.loadtxt(D/t/"tise/phase_shifts.dat")
def norm(x, y): return np.trapezoid(y**2, x)
def wrap_pi(d): return (d + np.pi/2) % np.pi - np.pi/2   # map to (-pi/2, pi/2]
def match_state(x, y_num, y_exact):
    """relative L2 error after fixing overall sign (eigenvectors are sign-arbitrary)."""
    s = np.sign(np.dot(y_num, y_exact)) or 1.0
    return np.sqrt(np.trapezoid((s*y_num - y_exact)**2, x) / np.trapezoid(y_exact**2, x)), s
def code_delta(k, psi, dpsi, R):
    return np.arctan(k*psi/dpsi) - k*R

def overlay(t, pairs, fname, title, ylab=r"$\phi_n(x)$"):
    n = len(pairs); fig, axes = plt.subplots(1, n, figsize=(4.2*n, 3.4))
    if n == 1: axes = [axes]
    for ax, (lab, x, ynum, yex) in zip(axes, pairs):
        ax.plot(x, yex, "k-", lw=2.5, alpha=.35, label="analytic")
        ax.plot(x, ynum, "C0-", lw=1, label="B-spline")
        ax.set_title(lab, fontsize=10); ax.set_xlabel("x")
    axes[0].set_ylabel(ylab); axes[0].legend(fontsize=8)
    fig.suptitle(title); fig.tight_layout(); fig.savefig(OUT/fname, dpi=110); plt.close(fig)

# ───────────────────────── 1. FREE PARTICLE (box [0,100]) ─────────────────────────
t = "free_particle"; L = 100.0
E = load_eigs(t); n = np.arange(1, len(E)+1)
E_ex = n**2 * np.pi**2 / (2*L**2)
rel = np.abs(E-E_ex)/E_ex
rows = [(int(i), float(E[i-1]), float(E_ex[i-1]), float(rel[i-1])) for i in list(range(1,11))+[15,20,30,40,50,len(E)]]
pairs=[]; wf=[]
for j in (1,2,3,5,10):
    x, y = load_state(t, j); yex = np.sqrt(2/L)*np.sin(j*np.pi*x/L)
    err, s = match_state(x, y, yex); wf.append((j, float(norm(x,y)), float(err)))
    if j in (1,2,3,10): pairs.append((f"n={j}", x, s*y, yex))
overlay(t, pairs, "free_particle_states.png", "Free particle in a box: $\\phi_n$ vs $\\sqrt{2/L}\\,\\sin(n\\pi x/L)$")
ps = load_ps(t); cont=[]; cpairs=[]
for i,(Ei, d, dd) in enumerate(ps, start=1):
    k = np.sqrt(2*Ei); x, y = load_cont(t, i); yex = np.sqrt(2/(np.pi*k))*np.sin(k*x)
    err, s = match_state(x, y, yex)
    cont.append((float(Ei), float(d), float(wrap_pi(d)), 0.0, float(dd), 0.0, float(err), float(np.abs(y).max()), float(np.sqrt(2/(np.pi*k)))))
    if i in (1,3): cpairs.append((f"E={Ei:.1f}", x[:150], s*y[:150], yex[:150]))
overlay(t, cpairs, "free_particle_continuum.png", "Free particle continuum: $\\psi_E$ vs $\\sqrt{2/\\pi k}\\,\\sin(kx)$", r"$\psi_E(x)$")
report[t] = dict(energies=rows, max_rel_err_first10=float(rel[:10].max()), max_rel_err_all=float(rel.max()),
                 argmax=int(np.argmax(rel)+1), wavefunctions=wf, continuum=cont)

# ───────────────────────── 2. FINITE SQUARE WELL ─────────────────────────
# V=-V0 on [0,a), 0 on [a,Lbox]; psi(0)=0 -> half-well = odd states of symmetric well of half-width a.
t = "finite_square_well"; V0=1.0; a=10.0; Lbox=100.0
E = load_eigs(t)
def f_inf(En):   # matching condition with decaying tail: K cot(Ka) = -kappa
    K = np.sqrt(2*(En+V0)); kap = np.sqrt(-2*En); return K/np.tan(K*a) + kap
def f_box(En):   # exact with hard wall at Lbox: outside psi = sinh(kappa (Lbox-x))
    K = np.sqrt(2*(En+V0)); kap = np.sqrt(-2*En)
    return K/np.tan(K*a) + kap/np.tanh(kap*(Lbox-a))
# bracket roots: Ka in ((m-1/2)pi, m pi), m=1..  (cot from 0 to -inf)
bound_inf=[]; bound_box=[]
m=1
while True:
    Klo=(m-0.5)*np.pi/a; Khi=m*np.pi/a
    Elo=Klo**2/2-V0; Ehi=min(Khi**2/2-V0, -1e-12)
    if Elo >= 0: break
    lo=Elo+1e-9; hi=Ehi-1e-12 if Ehi<0 else -1e-12
    try:
        bound_inf.append(brentq(f_inf, lo, hi))
    except ValueError: bound_inf.append(None)
    try:
        bound_box.append(brentq(f_box, lo, hi))
    except ValueError: bound_box.append(None)
    m+=1
nb = len(bound_inf)
Ka_max = np.sqrt(2*V0)*a
fsw_rows=[]
for i in range(nb):
    fsw_rows.append((i+1, float(E[i]), bound_inf[i], bound_box[i]))
# 5th state? Ka_max vs 4.5 pi
fsw_note = dict(Ka_max=float(Ka_max), four_half_pi=float(4.5*np.pi), n_bound_analytic=nb,
                first_positive_numeric=[float(e) for e in E[4:7]],
                box_states_expected=[float(j**2*np.pi**2/(2*Lbox**2)) for j in (1,2,3)])
# wavefunction overlay
pairs=[]; wf=[]
for j in range(1, nb+1):
    En = bound_inf[j-1] if bound_inf[j-1] is not None else E[j-1]
    K=np.sqrt(2*(En+V0)); kap=np.sqrt(-2*En)
    x,y = load_state(t,j)
    yex = np.where(x<a, np.sin(K*x), np.sin(K*a)*np.exp(-kap*(x-a)))
    yex /= np.sqrt(norm(x,yex))
    err,s = match_state(x,y,yex); wf.append((j,float(norm(x,y)),float(err)))
    pairs.append((f"n={j}, E={En:.4f}", x[:60], s*y[:60], yex[:60]))
overlay(t, pairs, "fsw_states.png", "Finite well bound states: $\\phi_n$ vs piecewise $\\sin(Kx)$ / $e^{-\\kappa x}$")
# continuum: delta = atan((k/K) tan(Ka)) - ka  (mod pi); exact dDelta/dE by numerical derivative
def delta_an(En):
    k=np.sqrt(2*En); K=np.sqrt(2*(En+V0)); return np.arctan((k/K)*np.tan(K*a)) - k*a
ps=load_ps(t); cont=[]; cpairs=[]
for i,(Ei,d,dd) in enumerate(ps, start=1):
    k=np.sqrt(2*Ei); K=np.sqrt(2*(Ei+V0))
    dan = delta_an(Ei); h=1e-5
    ddan = (delta_an(Ei+h)-delta_an(Ei-h))/(2*h)
    # exact continuum wavefunction, energy-normalised: outside sqrt(2/pi k) sin(kx+delta)
    x,y = load_cont(t,i)
    A = np.sqrt(2/(np.pi*k))
    inside = np.sin(K*x); outside = np.sin(k*x+dan)
    # continuity at a: scale inside by sin(ka+dan)/sin(Ka)
    yex = A*np.where(x<a, inside*np.sin(k*a+dan)/np.sin(K*a), outside)
    err,s = match_state(x,y,yex)
    cont.append((float(Ei), float(d), float(wrap_pi(d)), float(wrap_pi(dan)), float(dd), float(ddan), float(err), float(np.abs(y[x>a]).max()), float(A)))
    if i in (1,5): cpairs.append((f"E={Ei:.1f}", x[:120], s*y[:120], yex[:120]))
overlay(t, cpairs, "fsw_continuum.png", "Finite well continuum: B-spline vs analytic $\\sin(Kx)$ / $\\sin(kx+\\delta)$", r"$\psi_E(x)$")
report[t]=dict(bound=fsw_rows, note=fsw_note, wavefunctions=wf, continuum=cont)

# ───────────────────────── 3. HARMONIC OSCILLATOR ─────────────────────────
t="harmonic_oscillator"; E=load_eigs(t); n=np.arange(len(E)); E_ex=n+0.5
abserr=np.abs(E-E_ex)
rows=[(int(i), float(E[i]), float(E_ex[i]), float(abserr[i])) for i in list(range(10))+[15,20,30,40,50,60,70,80,88]]
def ho_state(nn, x):
    return (np.pi**-0.25)/np.sqrt(2.0**nn*factorial(nn))*eval_hermite(nn, x)*np.exp(-x**2/2)
pairs=[]; wf=[]
for j in (1,2,3,4,5,6,10,20,40):
    x,y=load_state(t,j); yex=ho_state(j-1,x); err,s=match_state(x,y,yex); wf.append((j,float(norm(x,y)),float(err)))
    if j in (1,2,3,6): pairs.append((f"n={j-1}", x[(x>-8)&(x<8)], (s*y)[(x>-8)&(x<8)], yex[(x>-8)&(x<8)]))
overlay(t, pairs, "ho_states.png", "Harmonic oscillator: $\\phi_n$ vs Hermite functions")
# where does 1e-6 accuracy end? turning point sqrt(2E) vs node spacing 0.5
report[t]=dict(energies=rows, wavefunctions=wf,
               n_below_1e6=int(np.sum(abserr<1e-6)), n_below_1e3=int(np.sum(abserr<1e-3)))

# ───────────────────────── 4. HYDROGEN l=1 ─────────────────────────
t="hydrogen"; E=load_eigs(t); l=1
nb_num = int(np.sum(E<0))
rows=[]
for i in range(nb_num):
    nn=i+2; rows.append((i+1, nn, float(E[i]), float(-0.5/nn**2), float(abs(E[i]+0.5/nn**2)/(0.5/nn**2))))
def u_nl(nn, ll, r):  # reduced radial function r*R_nl, Z=1, a.u.
    rho=2*r/nn
    N=np.sqrt((2/nn)**3*factorial(nn-ll-1)/(2*nn*factorial(nn+ll)))
    return r*N*np.exp(-rho/2)*rho**ll*genlaguerre(nn-ll-1, 2*ll+1)(rho)
pairs=[]; wf=[]
for j in (1,2,3,4,5,6):
    nn=j+1; x,y=load_state(t,j); yex=u_nl(nn,l,x); err,s=match_state(x,y,yex); wf.append((j,nn,float(norm(x,y)),float(err)))
    if j in (1,2,3,5): pairs.append((f"{nn}p  (n={nn}, l=1)", x, s*y, yex))
overlay(t, pairs, "hydrogen_states.png", "Hydrogen $\\ell=1$: $\\phi_n$ vs $u_{n1}(r)=rR_{n1}(r)$", r"$u(r)$")
# continuum: regular Coulomb function F_l(eta, kr), eta = -1/k (attractive Z=1)
mp.mp.dps=30
ps=load_ps(t); cont=[]; cpairs=[]; R=100.0
for i,(Ei,d,dd) in enumerate(ps, start=1):
    k=np.sqrt(2*Ei); eta=-1.0/k
    F=lambda r: float(mp.coulombf(l, eta, k*r))
    x,y=load_cont(t,i)
    yex=np.array([F(r) if r>0 else 0.0 for r in x])*np.sqrt(2/(np.pi*k))  # F ~ sin(...) asymptotically -> same energy normalisation
    err,s=match_state(x,y,yex)
    # apply the code's flat-asymptote delta formula to the exact Coulomb function at R
    psiR=F(R); dpsiR=float(mp.diff(lambda r: mp.coulombf(l, eta, k*r), R))
    d_flat_exact=code_delta(k, psiR, dpsiR, R)
    # what a true Coulomb wave would give: F ~ sin(kr - l pi/2 - eta ln(2kr) + sigma_l)
    sigma=float(mp.arg(mp.gamma(l+1+1j*eta)))
    d_coulomb_at_R = -l*np.pi/2 - eta*np.log(2*k*R) + sigma
    cont.append((float(Ei), float(d), float(wrap_pi(d)), float(wrap_pi(d_flat_exact)), float(wrap_pi(d_coulomb_at_R)), float(dd), float(err), float(np.abs(y[x>50]).max()), float(np.sqrt(2/(np.pi*k)))))
    if i in (1,5): cpairs.append((f"E={Ei:.1f}", x[:200], s*y[:200], yex[:200]))
overlay(t, cpairs, "hydrogen_continuum.png", "Hydrogen $\\ell=1$ continuum: B-spline vs Coulomb function $F_1(\\eta,kr)$", r"$\psi_E(r)$")
report[t]=dict(bound=rows, n_bound_numeric=nb_num, wavefunctions=wf, continuum=cont)

json.dump(report, open(OUT/"report.json","w"), indent=1)
print(json.dumps(report, indent=1))
