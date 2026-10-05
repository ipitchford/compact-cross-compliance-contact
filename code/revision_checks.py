#!/usr/bin/env python3
"""Supporting checks added in v2. Analytic contact theorems are not formalised.

Exact SymPy identities check the coordinate shift and elastic symbols. Numerical
quadrature and spectral samples are illustrations, not interval certificates.
"""
from __future__ import annotations
import csv
import json
from pathlib import Path
import mpmath as mp
import sympy as sp


def exact_coordinate_shift() -> dict:
    ag, zg = sp.symbols('alpha_G z_G', real=True)
    h, L, V, r = sp.symbols('h L V_star r', positive=True)
    alpha, z = ag + sp.log(V/L), zg + sp.log(V/L)
    age_g = ag-zg-h*(sp.exp(-ag)-r/L)
    age_here = alpha-z-h*(V*sp.exp(-alpha)-r)/L
    slip_g = ag-zg+h*r/L*(ag+sp.log(r/L))
    slip_here = alpha-z+h*r/L*(alpha+sp.log(r/V))
    if not (sp.simplify(age_g-age_here) == 0):
        raise AssertionError('check failed: code/revision_checks.py:23')
    if not (sp.simplify(sp.expand_log(slip_g-slip_here, force=True)) == 0):
        raise AssertionError('check failed: code/revision_checks.py:24')
    return {'ageing_residual': '0', 'slip_residual': '0',
            'shift': 'alpha = alpha_G + log(V_star/L)',
            'meaning': 'Same logarithmic implicit-Euler state equations; integrator not new.'}


def exact_elastic_checks() -> dict:
    z, k, mu = sp.symbols('depth k mu', positive=True)
    la = sp.symbols('lambda', real=True)
    a, b = sp.symbols('a b')
    d = (la+mu)/(la+3*mu)
    Nm = sp.Matrix([[-d,-d],[d,d]])
    Np = sp.Matrix([[d,-d],[d,-d]])
    if not (Nm*Nm == sp.zeros(2) and Np*Np == sp.zeros(2)):
        raise AssertionError('check failed: code/revision_checks.py:37')
    Hm = 2*k*mu/(la+3*mu)*sp.Matrix([[la+2*mu,-mu],[-mu,la+2*mu]])
    Hp = -2*k*mu/(la+3*mu)*sp.Matrix([[la+2*mu,mu],[mu,la+2*mu]])
    for sign, N, H in [(-1, Nm, Hm), (1, Np, Hp)]:
        X, Z = sp.exp(sign*k*z)*(sp.eye(2)+k*z*N)*sp.Matrix([a,b])
        div = -k*X + sp.diff(Z,z)
        res = [mu*(sp.diff(X,z,2)-k*k*X)+(la+mu)*k*div,
               mu*(sp.diff(Z,z,2)-k*k*Z)+(la+mu)*sp.diff(div,z)]
        if not (all(sp.simplify(x)==0 for x in res)):
            raise AssertionError('check failed: code/revision_checks.py:45')
        traction = sp.Matrix([-mu*(sp.diff(X,z)+k*Z), -la*div-2*mu*sp.diff(Z,z)])
        actual = traction.subs(z,0).jacobian([a,b])
        if not ((actual-H).applyfunc(sp.simplify) == sp.zeros(2)):
            raise AssertionError('check failed: code/revision_checks.py:48')
    H = sp.symbols('H', positive=True)
    q = k*H
    # Verify the top boundary without expanding exponentially growing factors.
    R0 = (sp.eye(2)-q*Np)*(sp.eye(2)+q*Nm)
    if not (((sp.eye(2)+q*Np)*R0-(sp.eye(2)+q*Nm)).applyfunc(sp.simplify)==sp.zeros(2)):
        raise AssertionError('check failed: code/revision_checks.py:53')
    Cinf = Hm.inv().applyfunc(sp.factor)
    if not (sp.simplify(Cinf[1,0] - 1/(2*k*(la+mu)))==0):
        raise AssertionError('check failed: code/revision_checks.py:55')
    # Reflection reverses only the normal/tangent cross sign in a common frame.
    J = sp.diag(1,-1)
    if not ((J*Cinf*J)[1,0] == -Cinf[1,0]):
        raise AssertionError('check failed: code/revision_checks.py:58')
    cp = sp.Rational(1,2)/(2+1)
    cm = sp.Rational(1,2)/(1+2)
    if not (cp == cm == sp.Rational(1,6)):
        raise AssertionError('check failed: code/revision_checks.py:61')
    # Low-frequency boundary ODE, independently integrated rather than expanded
    # from the compliance inverse (which is poorly conditioned at tiny k).
    T = sp.symbols('T', real=True)
    X0 = T*(H-z)/mu
    Z1prime0 = la*H*T/(mu*(la+2*mu))
    Z1second = -(la+mu)*T/(mu*(la+2*mu))
    Z1zero = -Z1prime0*H-Z1second*H**2/2
    Z1 = Z1zero+Z1prime0*z+Z1second*z*z/2
    if not (sp.simplify((la+2*mu)*sp.diff(Z1,z,2)-(la+mu)*sp.diff(X0,z))==0):
        raise AssertionError('check failed: code/revision_checks.py:70')
    if not (sp.simplify((la+2*mu)*sp.diff(Z1,z).subs(z,0)-la*X0.subs(z,0))==0):
        raise AssertionError('check failed: code/revision_checks.py:71')
    if not (sp.simplify(Z1.subs(z,H)) == 0):
        raise AssertionError('check failed: code/revision_checks.py:72')
    coeff = sp.factor(Z1zero/T)
    if not (sp.simplify(coeff-(mu-la)*H**2/(2*mu*(la+2*mu)))==0):
        raise AssertionError('check failed: code/revision_checks.py:74')
    delta = coeff.subs({la:2,mu:1,H:1})-coeff.subs({la:1,mu:2,H:1})
    if not (delta == -sp.Rational(7,40)):
        raise AssertionError('check failed: code/revision_checks.py:76')
    return {'Navier_modes_residuals': [0,0,0,0], 'traction_matrix_residuals': 'zero',
            'nilpotent_mode_matrices': 'N_minus^2=N_plus^2=0',
            'top_clamp_reflection_residual': 'zero',
            'outward_normal': '-e_3 for upper layer',
            'cartesian_u3_from_t1': '-i*xi_1/[2*(lambda+mu)*|xi|^2]',
            'outward_un_from_t1': '+i*xi_1/[2*(lambda+mu)*|xi|^2]',
            'relative_outward_normal_symbol': 'i*(chi_plus-chi_minus)*(xi dot s)/|xi|^2',
            'matched_chi': '1/6', 'one_layer_small_k_cross_coefficient': str(coeff),
            'matched_layer_difference_small_k_coefficient': str(delta),
            'scope': 'Exact identities only; compactness and contact existence require analytic proofs.'}


def layer_compliance(k, la, mu, H=1):
    k, la, mu, H = map(mp.mpf, (k,la,mu,H))
    if min(k,mu,H)<=0 or 3*la+2*mu<=0:
        raise ValueError('Require k, mu, H positive and positive strain energy')
    d=(la+mu)/(la+3*mu)
    Nm=mp.matrix([[-d,-d],[d,d]])
    Np=mp.matrix([[d,-d],[d,-d]])
    Hm=2*k*mu/(la+3*mu)*mp.matrix([[la+2*mu,-mu],[-mu,la+2*mu]])
    Hp=-2*k*mu/(la+3*mu)*mp.matrix([[la+2*mu,mu],[mu,la+2*mu]])
    q=k*H
    R=mp.exp(-2*q)*(mp.eye(2)-q*Np)*(mp.eye(2)+q*Nm)
    return (mp.eye(2)-R)*(Hm-Hp*R)**-1


def spectral_samples() -> list[dict]:
    rows=[]
    with mp.workdps(100):
        for k0 in ['0.0001','0.01','0.1','1','2','4','8','16','32','64']:
            k=mp.mpf(k0)
            A=layer_compliance(k,2,1)
            B=layer_compliance(k,1,2)
            D=layer_compliance(k,1,1)
            for C in (A,B,D):
                if not (abs(C[0,1]-C[1,0]) < mp.mpf('1e-75')):
                    raise AssertionError('check failed: code/revision_checks.py:112')
                if not (C[0,0]>0 and mp.det(C)>0):
                    raise AssertionError('check failed: code/revision_checks.py:113')
            matched=A[1,0]-B[1,0]
            mismatched=A[1,0]-D[1,0]
            rows.append({'k':k0,'matched_cross':mp.nstr(matched,70),
                         'k_times_matched_cross':mp.nstr(k*matched,70),
                         'mismatched_cross':mp.nstr(mismatched,70),
                         'k_times_mismatched_cross':mp.nstr(k*mismatched,70)})
        if not (abs(mp.mpf(rows[0]['matched_cross'])/mp.mpf(rows[0]['k'])+mp.mpf(7)/40)<mp.mpf('1e-7')):
            raise AssertionError('check failed: code/revision_checks.py:120')
        if not (abs(mp.mpf(rows[-1]['k_times_matched_cross']))<mp.mpf('1e-45')):
            raise AssertionError('check failed: code/revision_checks.py:121')
        if not (abs(mp.mpf(rows[-1]['k_times_mismatched_cross'])+mp.mpf(1)/12)<mp.mpf('1e-45')):
            raise AssertionError('check failed: code/revision_checks.py:122')
        if not (abs(mp.mpf(rows[3]['matched_cross']))>mp.mpf('0.001')):
            raise AssertionError('check failed: code/revision_checks.py:123')
    return rows


def numerical_energies(root_records: list[dict]) -> list[dict]:
    out=[]
    for rec in root_records:
        law,reg=rec['law'],rec['regularisation']
        load=mp.mpf(rec['load'])
        roots=[mp.mpf(x) for x in rec['roots']]
        passes=[]
        for precision in (40,70):
            with mp.workdps(precision):
                def state(r):
                    return -mp.log(mp.lambertw(mp.exp(r-3))) if law=='ageing' else (3-r*mp.log(r))/(1+r)
                def gfun(r): return 1+mp.log(r)+state(r)
                def mufun(r):
                    if r==0: return mp.mpf(0)
                    g=gfun(r)
                    return max(g,0) if reg=='truncated' else mp.asinh(mp.exp(g)/2)
                cutoff=(mp.lambertw((mp.e-1)*mp.exp(-4))/(mp.e-1)
                        if law=='ageing' else mp.lambertw(mp.exp(-4)))
                if not (abs(gfun(cutoff))<mp.mpf(10)**(-precision+5)):
                    raise AssertionError('check failed: code/revision_checks.py:145')
                vals=[]
                for r in roots:
                    lower=cutoff if reg=='truncated' else mp.mpf(0)
                    # Split at a fixed inner point and at the exact activation
                    # speed; the integral never crosses a truncated kink.
                    points=[lower]
                    for point in [mp.mpf('0.01'),mp.mpf('0.1'),mp.mpf(1),mp.mpf(3)]:
                        if lower < point < r: points.append(point)
                    points.append(r)
                    energy=r*r/16-load*r+mp.quad(mufun,points)
                    if law=='ageing':
                        t=mp.exp(-state(r)); gp=1/r-1/(1+t)
                    else:
                        gp=(1/r+1-3-mp.log(r))/(1+r)**2
                    sigma=mp.mpf(1) if reg=='truncated' else (mp.exp(gfun(r))/2)/mp.sqrt(1+(mp.exp(gfun(r))/2)**2)
                    curvature=mp.mpf(1)/8+sigma*gp
                    vals.append((+energy,+curvature))
                passes.append(vals)
        err=max(abs(passes[0][i][j]-passes[1][i][j]) for i in range(3) for j in range(2))
        if not (err<mp.mpf('1e-30')):
            raise AssertionError('check failed: code/revision_checks.py:165')
        if not ([mp.sign(v[1]) for v in passes[1]]==[1,-1,1]):
            raise AssertionError('check failed: code/revision_checks.py:166')
        winner=min(range(3),key=lambda i:passes[1][i][0])+1
        if not (winner==(1 if law=='ageing' else 3)):
            raise AssertionError('check failed: code/revision_checks.py:168')
        out.append({'law':law,'regularisation':reg,
                    'roots':[mp.nstr(x,25) for x in roots],
                    'energies':[mp.nstr(x[0],30) for x in passes[1]],
                    'scalar_curvatures':[mp.nstr(x[1],30) for x in passes[1]],
                    'lowest_among_displayed_roots':winner,
                    'precision_comparison_max_difference':mp.nstr(err,6),
                    'assurance':'Numerical quadrature at 40/70 digits, not interval certification or global-minimizer classification.'})
    return out


def run(root_records: list[dict], out_dir: Path) -> dict:
    shift=exact_coordinate_shift()
    elastic=exact_elastic_checks()
    spectrum=spectral_samples()
    energies=numerical_energies(root_records)
    with (out_dir/'matched_layer_spectrum.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(spectrum[0])); writer.writeheader(); writer.writerows(spectrum)
    report={'status':'PASS','groups':[
        {'name':'GKP logarithmic-state coordinate identity','kind':'exact symbolic'},
        {'name':'two-layer modes, tractions, reflection and small-k coefficient','kind':'exact symbolic'},
        {'name':'matched and mismatched finite-layer compliance','kind':'100-digit numerical illustration'},
        {'name':'stationary-root energies and scalar curvatures','kind':'40/70-digit numerical illustration'}],
        'coordinate_shift':shift,'elastic':elastic,'finite_layer_spectrum':spectrum,'energies':energies,
        'analytic_scope':'No computation proves the pseudomonotone VI theorem or the curved-boundary compactness argument.'}
    (out_dir/'revision_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    (out_dir/'revision_verification.log').write_text('\n'.join('PASS ['+x['kind']+'] '+x['name'] for x in report['groups'])+'\n')
    return report
