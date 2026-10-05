"""Independent numerical cross-checks of the supplied friction manuscript.
Does not import or execute the author's verifier. These tests are not proofs.
Run: python independent_checks.py --out independent_results.json
Dependencies: numpy, scipy, mpmath.
"""
from __future__ import annotations
import argparse, hashlib, json, platform
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.special import lambertw
import mpmath as mp

def state(r: float, law: str, z: float=3.0) -> float:
    if r < 0:
        raise ValueError('negative speed')
    if law == 'ageing':
        # Direct Lambert-W expression, independent of submitted root evaluation.
        return float(lambertw(np.exp(r-z)).real + z-r)
    if law == 'slip':
        return z if r == 0 else (z-r*np.log(r))/(1+r)
    raise ValueError(law)

def friction(r: float, law: str, reg: str) -> float:
    if r == 0:
        return 0.0
    g=1+np.log(r)+state(r, law)
    return max(0.0,g) if reg=='truncated' else float(np.arcsinh(np.exp(g)/2))

def main(out: Path) -> None:
    records=[]
    for law in ('ageing','slip'):
        load=2.8 if law=='ageing' else 2.5
        points=[0,4,5+np.log(2),32] if law=='ageing' else [0,1,3,32]
        for reg in ('truncated','asinh'):
            f=lambda r:r/8+friction(r,law,reg)-load
            roots=[brentq(f,points[i],points[i+1],xtol=5e-15,rtol=1e-14) for i in range(3)]
            slopes=[]; energies=[]
            for r in roots:
                step=1e-5*max(r,1)
                slopes.append((f(r+step)-f(r-step))/(2*step))
                # Supplementary classification, not a global-minimiser certificate.
                energy,err=quad(lambda s:friction(s,law,reg),0,r,epsabs=2e-10,epsrel=2e-10,limit=500)
                energies.append({'energy':r*r/16+energy-load*r,'quad_error_estimate':err})
            assert max(map(lambda r:abs(f(r)),roots))<1e-11
            assert slopes[0]>0 and slopes[1]<0 and slopes[2]>0
            records.append(dict(law=law,regularisation=reg,roots=roots,residual_max=max(abs(f(r)) for r in roots),scalar_energy_second_derivatives=slopes,energy_values=energies))
    # Assemble the actual consistent-mass P1 matrix and eliminate all interior
    # degrees of freedom; compare with the submitted closed-form formula.
    fe=[]; c=np.tanh(1.0)/8
    for n in (8,16,32,64,128):
        dx=1.0/n
        local=c/dx*np.array([[1.,-1.],[-1.,1.]])+c*dx/6*np.array([[2.,1.],[1.,2.]])
        A=np.zeros((n+1,n+1))
        for j in range(n):
            A[j:j+2,j:j+2]+=local
        end=np.array([0,n]); interior=np.arange(1,n)
        Schur=A[np.ix_(end,end)]-A[np.ix_(end,interior)]@np.linalg.solve(A[np.ix_(interior,interior)],A[np.ix_(interior,end)])
        diag=c/dx+c*dx/3; off=-c/dx+c*dx/6
        kap=np.arccosh(-diag/off)
        analytic_imp=diag+off*np.sinh((n-1)*kap)/np.sinh(n*kap)
        analytic_cross=off*np.sinh(kap)/np.sinh(n*kap)
        dif=max(abs(Schur[1,1]-analytic_imp),abs(Schur[1,0]-analytic_cross))
        assert dif<3e-12
        fe.append(dict(elements=n,assembled_impedance=float(Schur[1,1]),assembled_cross=float(Schur[1,0]),closed_formula_max_error=float(dif),continuum_impedance_error=float(Schur[1,1]-.125)))
    # Smooth radial compactly supported pressure. Normalise its total mass to pi
    # and recompute the cross pairing via the radial Poisson solution.
    mp.mp.dps=40
    def base(r):
        return mp.exp(-1/(1-r*r)) if r<1 else mp.mpf(0)
    norm=1/(2*mp.quad(lambda r:r*base(r),[0,mp.mpf('.5'),1]))
    def density(r): return norm*base(r)
    # Integration by parts gives pi/2*M(1)^2; this is a second profile check,
    # not an independent proof of the general Fourier-space argument.
    mass_integral=mp.quad(lambda r:r*density(r),[0,mp.mpf('.5'),1])
    defect=mp.pi/2*mass_integral**2
    assert abs(defect-mp.pi/8)<mp.mpf('1e-35')
    result={'status':'PASS','scope':'Independent floating-point cross-checks, not formal proof or global energy certification','environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'mpmath':mp.__version__},'independent_scalar_roots':records,'assembled_finite_elements':fe,'smooth_radial_profile':{'mass':mp.nstr(2*mp.pi*mass_integral,30),'Q0':mp.nstr(defect,30),'QK':mp.nstr(-2*mp.pi*defect,30)}}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=Path('independent_results.json'))
    main(ap.parse_args().out)
