"""References for the three new tests (interior_singularity, right_edge_singularity, case3_irregular_tail)."""
import numpy as np, json
from pathlib import Path
from scipy.optimize import brentq
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
from scipy.interpolate import BSpline
import mpmath as mp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
mp.mp.dps = 25
D = Path("/tmp/BSpline-CPP/1D-QM-Playground/data"); OUT = Path("./verify_out"); OUT.mkdir(exist_ok=True)
rep = {}
def eig(t): return np.loadtxt(D/t/"tise/eigenvalues.dat")[:,1]
def state(t,j): a=np.loadtxt(D/t/f"tise/eigenstate_{j:03d}.dat"); return a[:,0],a[:,1]
def cont(t,j): a=np.loadtxt(D/t/f"tise/continuum_state_{j:03d}.dat"); return a[:,0],a[:,1]

# ---------- generic FD reference solver: -1/2 psi'' + V psi = E psi, psi=0 at both ends ----------
def fd_eigs(V, a, b, N=60000, k=12):
    x = np.linspace(a, b, N+2)[1:-1]; h = x[1]-x[0]
    main = 1.0/h**2 + V(x); off = -0.5/h**2*np.ones(N-1)
    H = diags([off, main, off], [-1,0,1], format="csc")
    w = eigsh(H, k=k, sigma=0.0, which="LM", return_eigenvectors=False)
    return np.sort(w)
def fd_richardson(V,a,b,k=12):
    e1 = fd_eigs(V,a,b,N=40000,k=k); e2 = fd_eigs(V,a,b,N=80000,k=k)
    return (4*e2 - e1)/3   # h^2 extrapolation

# ---------- generic B-spline solver (to demonstrate the proper singular-join treatment) ----------
def bspline_eigs(V, knots_interior, a, b, order, drop_x=None, nq=None, k_return=12):
    """order = number of coefficients per interval (B-spline 'order' k, degree k-1). knots_interior: list of interior breakpoints
    with multiplicity already expanded. drop_x: point where basis functions nonzero there are dropped (psi(drop_x)=0)."""
    deg = order-1
    t = np.r_[[a]*order, knots_interior, [b]*order]
    nb = len(t)-order
    # quadrature: Gauss-Legendre on each distinct interval
    br = np.unique(t)
    gx, gw = np.polynomial.legendre.leggauss(order+2)
    X=[];W=[]
    for lo,hi in zip(br[:-1],br[1:]):
        X.append(0.5*(hi-lo)*gx+0.5*(hi+lo)); W.append(0.5*(hi-lo)*gw)
    X=np.concatenate(X); W=np.concatenate(W)
    B  = BSpline.design_matrix(X, t, deg).toarray()
    dB = np.zeros_like(B)
    for i in range(nb):
        c=np.zeros(nb); c[i]=1; dB[:,i]=BSpline(t,c,deg).derivative()(X)
    S = (B*W[:,None]).T @ B
    H = 0.5*(dB*W[:,None]).T @ dB + (B*(W*V(X))[:,None]).T @ B
    keep = np.ones(nb,bool); keep[0]=keep[-1]=False
    if drop_x is not None:
        vals = BSpline.design_matrix(np.array([drop_x]), t, deg).toarray()[0]
        keep &= (np.abs(vals) < 1e-14)
    idx=np.where(keep)[0]
    from scipy.linalg import eigh
    w = eigh(H[np.ix_(idx,idx)], S[np.ix_(idx,idx)], eigvals_only=True)
    return np.sort(w)[:k_return], int(keep.sum())

# =====================================================================================
# 1. interior_singularity: V=0 on [0,20), 1/(x-20) on (20,40]. Exact: psi(20)=0 splits the domain.
# =====================================================================================
t="interior_singularity"; E=eig(t)
box = np.array([n**2*np.pi**2/(2*20**2) for n in range(1,8)])
# Coulomb side: u(s)=F_0(eta=1/k, k s), s=x-20, u(20)=0 -> F_0(1/k, 20k)=0
def F0(Ei): k=np.sqrt(2*Ei); return float(mp.coulombf(0, 1.0/k, 20*k))
Es=np.linspace(0.02,1.2,6000); v=np.array([F0(e) for e in Es]); idx=np.where(np.diff(np.sign(v))!=0)[0]
coul=np.array([brentq(F0,Es[i],Es[i+1]) for i in idx])
exact=np.sort(np.r_[box,coul]); side=['box' if e in box else 'coulomb' for e in exact]
# where does each solver state live?
loc=[]
for j in range(1,13):
    x,y=state(t,j); left=np.trapezoid(y[x<20]**2,x[x<20]); right=np.trapezoid(y[x>20]**2,x[x>20])
    dead=x[(np.abs(y)<1e-13)&(x>18)&(x<22)]
    loc.append((j,float(E[j-1]),'left' if left>right else 'right', float(dead.min()) if len(dead) else None, float(dead.max()) if len(dead) else None))
# proper treatment demo: uniform h=1 grid, order 8, knot at 20 with multiplicity 7 (=order-1), drop the one spline nonzero at 20
V_is = lambda x: np.where(x>20, 1.0/np.maximum(x-20,1e-300), 0.0)
kn = [float(i) for i in range(1,40)]
kn_mult = sorted(kn + [20.0]*6)           # 20 already present once -> total multiplicity 7
proper,nkeep = bspline_eigs(V_is, kn_mult, 0.0, 40.0, 8, drop_x=20.0)
naive,_  = bspline_eigs(V_is, kn, 0.0, 40.0, 8, drop_x=None)   # simple knots, nothing dropped (pre-8236239 behaviour, approx.)
rep[t]=dict(solver=[float(e) for e in E[:12]], exact=[(float(e),s) for e,s in zip(exact[:12],side[:12])],
            box=[float(b) for b in box], coulomb=[float(c) for c in coul[:8]], location=loc,
            proper_bspline=[float(p) for p in proper[:12]], nkeep_proper=nkeep, naive_bspline=[float(p) for p in naive[:12]])
# figure: solver states 1,3 vs exact
fig,ax=plt.subplots(1,3,figsize=(13,3.4))
x,y=state(t,1); ax[0].plot(x,np.sqrt(2/20)*np.sin(np.pi*x/20)*(x<20),'k-',lw=2.5,alpha=.35,label='exact  n=1 box state'); ax[0].plot(x,np.sign(y[50])*y,'C0-',lw=1,label='solver state 1'); ax[0].axvspan(19,21,color='C3',alpha=.15); ax[0].set_title(f'state 1: E={E[0]:.5f}  (exact {box[0]:.5f})',fontsize=9); ax[0].legend(fontsize=7)
x,y=state(t,3); k=np.sqrt(2*coul[0]); fe=np.array([float(mp.coulombf(0,1/k,k*(xx-20))) if xx>20 else 0 for xx in x]); fe/=np.sqrt(np.trapezoid(fe**2,x))
ax[1].plot(x,fe,'k-',lw=2.5,alpha=.35,label='exact  Coulomb-side state 1'); ax[1].plot(x,np.sign(y[np.argmax(np.abs(y))])*np.sign(fe[np.argmax(np.abs(y))])*y,'C0-',lw=1,label='solver state 3'); ax[1].axvspan(19,21,color='C3',alpha=.15); ax[1].set_title(f'state 3: E={E[2]:.5f}  (exact {coul[0]:.5f})',fontsize=9); ax[1].legend(fontsize=7)
ax[2].plot(range(1,13),E[:12],'o',label='solver 8236239'); ax[2].plot(range(1,13),exact[:12],'k_',ms=14,label='exact'); ax[2].plot(range(1,13),proper[:12],'x',label='B-splines, mult-7 knot at 20 + drop 1'); ax[2].set_xlabel('state'); ax[2].set_ylabel('E'); ax[2].legend(fontsize=7); ax[2].set_title('spectrum',fontsize=9)
fig.suptitle('interior_singularity: V=0 | 1/(x−20)'); fig.tight_layout(); fig.savefig(OUT/'interior_singularity.png',dpi=110); plt.close(fig)

# =====================================================================================
# 2. right_edge_singularity: V=1/(100-x) on [0,100). Exact bound-box states: F_0(1/k, 100k)=0
# =====================================================================================
t="right_edge_singularity"; E=eig(t)
def F0R(Ei): k=np.sqrt(2*Ei); return float(mp.coulombf(0, 1.0/k, 100*k))
Es=np.linspace(0.005,0.2,8000); v=np.array([F0R(e) for e in Es]); idx=np.where(np.diff(np.sign(v))!=0)[0]
exR=np.array([brentq(F0R,Es[i],Es[i+1]) for i in idx])
rows=[(j+1,float(E[j]),float(exR[j]),float((E[j]-exR[j])/exR[j])) for j in range(min(12,len(exR)))]
wf=[]
for j in (1,2,5):
    x,y=state(t,j); k=np.sqrt(2*exR[j-1]); fe=np.array([float(mp.coulombf(0,1/k,k*(100-xx))) for xx in x]); fe/=np.sqrt(np.trapezoid(fe**2,x))
    s=np.sign(np.dot(y,fe)); wf.append((j,float(np.sqrt(np.trapezoid((s*y-fe)**2,x)))))
# continuum: what does the code produce at the wall?
ps=np.loadtxt(D/t/"tise/phase_shifts.dat"); cinfo=[]
for i,(Ei,d,dd) in enumerate(ps,1):
    x,y=cont(t,i); cinfo.append((float(Ei),float(d),float(y[-1]),float(np.abs(y[x<90]).max())))
rep[t]=dict(rows=rows,wf=wf,continuum=cinfo, n_exact_found=len(exR))
fig,ax=plt.subplots(1,2,figsize=(9,3.4))
x,y=state(t,1); k=np.sqrt(2*exR[0]); fe=np.array([float(mp.coulombf(0,1/k,k*(100-xx))) for xx in x]); fe/=np.sqrt(np.trapezoid(fe**2,x)); s=np.sign(np.dot(y,fe))
ax[0].plot(x,fe,'k-',lw=2.5,alpha=.35,label='exact F₀(1/k, k(100−x))'); ax[0].plot(x,s*y,'C0-',lw=1,label='solver state 1'); ax[0].legend(fontsize=7); ax[0].set_title(f'state 1: E={E[0]:.6f} (exact {exR[0]:.6f})',fontsize=9)
x,y=cont(t,1); ax[1].plot(x,y,'C0-',lw=1); ax[1].axvline(100,color='C3'); ax[1].set_title(f'continuum state 1 (E={ps[0][0]:.3f}) — nonzero at the singular wall',fontsize=9)
fig.suptitle('right_edge_singularity: V=1/(100−x)'); fig.tight_layout(); fig.savefig(OUT/'right_edge_singularity.png',dpi=110); plt.close(fig)

# =====================================================================================
# 3. case3_irregular_tail: V=x^-1.5 on [0.1,50], solver tapers V to 0 over [R-4.99, R]
# =====================================================================================
t="case3_irregular_tail"; E=eig(t); R=50.0; dl=4.99
def W(x): d=x-R; return np.where(d<=-dl,1.0,np.where(d>=0,0.0,np.sin(np.pi/2*d/dl)**2))
Vraw=lambda x: x**-1.5; Vtap=lambda x: x**-1.5*W(x)
e_raw=fd_richardson(Vraw,0.1,50.0); e_tap=fd_richardson(Vtap,0.1,50.0)
rows=[(j+1,float(E[j]),float(e_tap[j]),float(E[j]-e_tap[j]),float(e_raw[j]),float(E[j]-e_raw[j])) for j in range(12)]
rep[t]=dict(rows=rows, delta=dl, V_at_R=float(R**-1.5), V_at_R_minus_delta=float((R-dl)**-1.5))
fig,ax=plt.subplots(1,2,figsize=(9,3.4))
xx=np.linspace(0.1,50,2000); ax[0].plot(xx,Vraw(xx),'k-',alpha=.4,label='V = x^-1.5'); ax[0].plot(xx,Vtap(xx),'C0-',label='tapered (as solved)'); ax[0].set_ylim(0,0.05); ax[0].set_xlim(30,50); ax[0].legend(fontsize=7); ax[0].set_title('the taper: last 4.99 bohr',fontsize=9)
ax[1].plot(range(1,13),E[:12]-e_tap[:12],'o-',label='solver − FD(tapered V)'); ax[1].plot(range(1,13),E[:12]-e_raw[:12],'s-',label='solver − FD(raw V)'); ax[1].axhline(0,color='k',lw=.5); ax[1].set_yscale('symlog',linthresh=1e-7); ax[1].legend(fontsize=7); ax[1].set_xlabel('state'); ax[1].set_title('energy differences',fontsize=9)
fig.suptitle('case3_irregular_tail: V = x^{-3/2} on [0.1, 50]'); fig.tight_layout(); fig.savefig(OUT/'case3_irregular_tail.png',dpi=110); plt.close(fig)

json.dump(rep,open(OUT/'report_new.json','w'),indent=1); print(json.dumps(rep,indent=1))
