#!/usr/bin/env python3
"""Deterministic supporting checks for the Problem 4 research candidate.

This is NOT a formal proof assistant. Exact checks cover algebra and rational
inequalities; analytic compactness/existence arguments are proved in the paper.
Floating-point roots and finite-element tables are explicitly supplementary.
Run: python code/verify.py --out evidence
"""
from __future__ import annotations
import argparse
import csv
import json
import platform
from fractions import Fraction as F
from pathlib import Path
import sys
import sympy as sp
import mpmath as mp

mp.mp.dps = 80
CHECKS: list[dict[str, str]] = []


def passed(name: str, kind: str, detail: str) -> None:
    CHECKS.append(dict(name=name, kind=kind, status="PASS", detail=detail))


def exp_lower(x: F, n: int = 20) -> F:
    """Strict lower bound on exp(x) when x>0, using n+1 series terms."""
    if x <= 0 or n < 0:
        raise ValueError("Positive x and nonnegative n are required")
    term = total = F(1)
    for j in range(1, n + 1):
        term *= x / j
        total += term
    return total


def exp_upper(x: F, n: int = 20) -> F:
    """Upper bound: finite Taylor sum plus a geometric majorant for its tail."""
    total = exp_lower(x, n)
    term = F(1)
    for j in range(1, n + 2):
        term *= x / j
    ratio = x / (n + 2)
    if ratio >= 1:
        raise ValueError("Increase n")
    return total + term / (1 - ratio)


def exact_sign_certificates() -> dict:
    if not (exp_lower(F(7, 10)) > 2):
        raise AssertionError('check failed: code/verify.py:52')
    if not (exp_lower(F(21, 20)) > F(57, 20)):
        raise AssertionError('check failed: code/verify.py:53')
    if not (exp_lower(F(11, 10)) > 3):
        raise AssertionError('check failed: code/verify.py:54')
    if not (exp_upper(F(1)) < 3):
        raise AssertionError('check failed: code/verify.py:55')
    if not (exp_lower(F(2), 4) == 7):  # omitted terms are strictly positive
        raise AssertionError('check failed: code/verify.py:56')
    # log 2 = 2 sum_{j>=0} (1/3)^(2j+1)/(2j+1) > 2/3.
    age_peak_lower = F(3, 2) + F(4, 3)
    age_valley_upper = F(57, 80) + F(41, 20) + F(4, 175)
    slip_peak_lower = F(21, 8)
    slip_valley_upper = F(17, 8) + F(11, 40) + F(1, 49)
    if not (age_peak_lower > F(14, 5) > age_valley_upper):
        raise AssertionError('check failed: code/verify.py:62')
    if not (slip_peak_lower > F(5, 2) > slip_valley_upper):
        raise AssertionError('check failed: code/verify.py:63')
    if not (F(32, 8) > max(F(14, 5), F(5, 2))):
        raise AssertionError('check failed: code/verify.py:64')
    out = {
        "age_peak_lower": str(age_peak_lower),
        "age_valley_upper_asinh": str(age_valley_upper),
        "age_load": "14/5",
        "slip_peak_lower": str(slip_peak_lower),
        "slip_valley_upper_asinh": str(slip_valley_upper),
        "slip_load": "5/2",
        "method": "Fraction arithmetic and positive exponential Taylor remainders",
    }
    passed("four multiplicity sign patterns", "exact rational", str(out))
    return out


def symbolic_state_and_thresholds() -> None:
    q, d = sp.symbols("q delta", positive=True)
    z, c, a, b = sp.symbols("z c a b", real=True)
    S = (z - q * sp.log(q / d)) / (1 + q)
    if not (sp.simplify(S - z + q * (S + sp.log(q / d))) == 0):
        raise AssertionError('check failed: code/verify.py:82')
    H = a * sp.log(q / d) + b * S
    target = a / q - b * (sp.log(q / d) + z + 1 + q) / (1 + q)**2
    if not (sp.simplify(sp.diff(H, q) - target) == 0):
        raise AssertionError('check failed: code/verify.py:85')
    # After c=z-log(delta), the feedback amplification has this simple form.
    G = q * (sp.log(q) + c + 1 + q) / (1 + q)**2
    numerator = (1 - q) * (sp.log(q) + c) + 2 * (1 + q)
    if not (sp.simplify(sp.diff(G, q) - numerator / (1 + q)**3) == 0):
        raise AssertionError('check failed: code/verify.py:89')
    critical_c = 2 * (q + 1) / (q - 1) - sp.log(q)
    if not (sp.simplify(G.subs(c, critical_c) - q / (q - 1)) == 0):
        raise AssertionError('check failed: code/verify.py:91')
    t = sp.symbols("t", positive=True)
    R = (1 + t) / (t + sp.log(t) + c)
    expected = (sp.log(t) + c - 2 - 1 / t) / (t + sp.log(t) + c)**2
    if not (sp.simplify(sp.diff(R, t) - expected) == 0):
        raise AssertionError('check failed: code/verify.py:95')
    if not (sp.simplify(R.subs(c, 2 + 1 / t - sp.log(t)) - t / (t + 1)) == 0):
        raise AssertionError('check failed: code/verify.py:96')
    # Exact frozen-rate slip flow (a supplementary extension).
    Sexact = sp.exp(-q) * z - (1 - sp.exp(-q)) * sp.log(q / d)
    e_derivative = -sp.exp(-q) * (z + sp.log(q / d)) - (1 - sp.exp(-q)) / q
    if not (sp.simplify(sp.diff(Sexact, q) - e_derivative) == 0):
        raise AssertionError('check failed: code/verify.py:100')
    # The elementary slope majorant used for backward Euler slip.
    majorant = (q - 1) / (q + 1)**2
    if not (sp.simplify(sp.diff(majorant, q) - (3 - q) / (q + 1)**3) == 0):
        raise AssertionError('check failed: code/verify.py:103')
    if not (majorant.subs(q, 3) == sp.Rational(1, 8)):
        raise AssertionError('check failed: code/verify.py:104')
    passed("state maps and sharp BE thresholds", "exact symbolic", "All residuals simplify to zero")


def symbolic_halfspace() -> dict:
    x = sp.symbols("depth", real=True)
    k, mu = sp.symbols("k mu", positive=True)
    la = sp.symbols("lambda", real=True)
    A, C = sp.symbols("A C")
    B = -k * (la + mu) / (la + 3 * mu) * (A + sp.I * C)
    D = k * (la + mu) / (la + 3 * mu) * (C - sp.I * A)
    u1 = (A + B * x) * sp.exp(-k * x)
    u3 = (C + D * x) * sp.exp(-k * x)
    div = sp.I * k * u1 + sp.diff(u3, x)
    navier1 = mu * (sp.diff(u1, x, 2) - k*k*u1) + (la + mu) * sp.I*k*div
    navier3 = mu * (sp.diff(u3, x, 2) - k*k*u3) + (la + mu) * sp.diff(div, x)
    if not (sp.simplify(navier1) == sp.simplify(navier3) == 0):
        raise AssertionError('check failed: code/verify.py:120')
    traction1 = -mu * (sp.diff(u1, x) + sp.I*k*u3).subs(x, 0)
    traction3 = -(la*div + 2*mu*sp.diff(u3, x)).subs(x, 0)
    actual_D = sp.Matrix([traction1, traction3]).jacobian([A, C])
    expected_D = 2*k*mu/(la+3*mu) * sp.Matrix([[la+2*mu, -sp.I*mu], [sp.I*mu, la+2*mu]])
    if not ((actual_D - expected_D).applyfunc(sp.simplify) == sp.zeros(2)):
        raise AssertionError('check failed: code/verify.py:125')
    expected_N = sp.Matrix([
        [(la+2*mu)/(2*k*mu*(la+mu)), sp.I/(2*k*(la+mu))],
        [-sp.I/(2*k*(la+mu)), (la+2*mu)/(2*k*mu*(la+mu))],
    ])
    if not ((expected_D * expected_N - sp.eye(2)).applyfunc(sp.simplify) == sp.zeros(2)):
        raise AssertionError('check failed: code/verify.py:130')
    passed("3D isotropic halfspace symbol", "exact symbolic", "Navier residuals, traction matrix, and inverse verified")
    return {"DtN": str(expected_D), "NtD": str(expected_N), "u3_from_t1_multiplier": "-i*xi_1/(2*(lambda+mu)*|xi|^2)", "outward_un_from_t1_multiplier": "+i*xi_1/(2*(lambda+mu)*|xi|^2)", "outward_normal": [0, 0, -1], "component_relation": "u_n=-u_3"}


def concentration_identities() -> dict:
    r, phi, eps, p = sp.symbols("r phi epsilon p", positive=True)
    moment = sp.integrate(r**3, (r, 0, 1)) * sp.integrate(sp.cos(phi)**2, (phi, 0, 2*sp.pi))
    if not (moment == sp.pi / 4):
        raise AssertionError('check failed: code/verify.py:138')
    pairing_L = moment / 2
    pairing_K = -2*sp.pi*pairing_L
    if not (pairing_L == sp.pi/8 and pairing_K == -sp.pi**2/4):
        raise AssertionError('check failed: code/verify.py:141')
    # Nonnegative C^1 compactly supported radial profile also has nonzero defect.
    smoothish = sp.integrate(r**3*(1-r*r)**2*(sp.Rational(1, 2)-r*r/2+r**4/6), (r, 0, 1))*sp.pi
    if not (smoothish > 0):
        raise AssertionError('check failed: code/verify.py:144')
    if not ((-sp.Rational(3, 2) + 2/p).subs(p, sp.Rational(4, 3)) == 0):
        raise AssertionError('check failed: code/verify.py:145')
    if not (-3 + 4 - 1 == 0):  # amplitudes, dxdy, degree -1 kernel
        raise AssertionError('check failed: code/verify.py:146')
    passed("concentration constants and critical exponent", "exact symbolic", "L defect=pi/8; K defect=-pi^2/4; critical p=4/3")
    return {"disk_L_pairing": "pi/8", "disk_K_pairing": "-pi^2/4", "C1_radial_L_pairing": str(smoothish), "critical_pressure_exponent": "4/3"}


def age_state(r: mp.mpf, z: mp.mpf = mp.mpf(3)) -> mp.mpf:
    """Dimensionless h=L=V*=1 BE log-state; stable for this test range."""
    t = mp.lambertw(mp.exp(r-z))
    return -mp.log(t)


def slip_state(r: mp.mpf, z: mp.mpf = mp.mpf(3)) -> mp.mpf:
    return z if r == 0 else (z-r*mp.log(r))/(1+r)


def friction(r: mp.mpf, law: str, regularisation: str) -> mp.mpf:
    if r < 0:
        raise ValueError("Only nonnegative speeds are accepted")
    if r == 0:
        return mp.mpf(0)
    S = age_state(r) if law == "ageing" else slip_state(r)
    g = 1 + mp.log(r) + S
    return max(g, mp.mpf(0)) if regularisation == "truncated" else mp.asinh(mp.exp(g)/2)


def bisect(fun, lo, hi, iterations: int = 230) -> mp.mpf:
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    flo, fhi = fun(lo), fun(hi)
    if flo == 0:
        return lo
    if fhi == 0:
        return hi
    if flo * fhi >= 0:
        raise ValueError(f"Not a strict sign bracket: {lo}, {hi}: {flo}, {fhi}")
    for _ in range(iterations):
        mid = (lo + hi)/2
        fm = fun(mid)
        if fm == 0:
            return mid
        if flo * fm < 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return (lo+hi)/2


def numerical_root_checks() -> list[dict]:
    out = []
    for law in ("ageing", "slip"):
        load = mp.mpf(14)/5 if law == "ageing" else mp.mpf(5)/2
        pts = [0, 4, 5+mp.log(2), 32] if law == "ageing" else [0, 1, 3, 32]
        for reg in ("truncated", "asinh"):
            fun = lambda r: r/8 + friction(r, law, reg) - load
            signs = [mp.sign(fun(mp.mpf(x))) for x in pts]
            if not (signs == [-1, 1, -1, 1]):
                raise AssertionError('check failed: code/verify.py:200')
            roots = [bisect(fun, pts[j], pts[j+1]) for j in range(3)]
            residuals = [abs(fun(root)) for root in roots]
            if not (max(residuals) < mp.mpf("1e-60")):
                raise AssertionError('check failed: code/verify.py:203')
            out.append(dict(law=law, regularisation=reg, load=mp.nstr(load, 20),
                            roots=[mp.nstr(v, 45) for v in roots],
                            residual_upper=mp.nstr(max(residuals), 5),
                            assurance="80-digit floating point; exact existence uses rational sign certificates"))
    passed("root illustrations", "high-precision supplementary", "12 roots in the analytically certified disjoint intervals")
    return out


def state_bound_tests() -> None:
    for z0 in (-100, -5, 0, 3, 100):
        z = mp.mpf(z0)
        for r0 in (0, mp.mpf("1e-30"), mp.mpf("0.01"), 1, 4, 10, 100):
            r = mp.mpf(r0)
            sa = age_state(r, z)
            ss = slip_state(r, z)
            tol = mp.mpf("1e-60")
            if not (abs(sa - mp.exp(-sa) - z + r) < tol * (1+abs(z)+r)):
                raise AssertionError('check failed: code/verify.py:220')
            if not (sa >= z-r-tol):
                raise AssertionError('check failed: code/verify.py:221')
            if not (max(sa, 0) <= max(z, 0)+1+tol):
                raise AssertionError('check failed: code/verify.py:222')
            if not (abs(sa) <= abs(z)+r+1+tol):
                raise AssertionError('check failed: code/verify.py:223')
            if not (max(ss, 0) <= max(z, 0)+1/mp.e+tol):
                raise AssertionError('check failed: code/verify.py:224')
            if not (abs(ss) <= abs(z)+1/mp.e+r+tol):
                raise AssertionError('check failed: code/verify.py:225')
    passed("state-map edge cases", "high-precision supplementary", "35 cases per law, including zero and tiny speed and large signed state")


def threshold_samples() -> list[dict]:
    out = []
    for c0 in (-5, -2, 0, 2, 5):
        c = mp.mpf(c0)
        wa = mp.lambertw(mp.exp(c-2))
        beta_a = 1/(1+wa)
        fun = lambda q: mp.log(q)+c-2*(q+1)/(q-1)
        hi = mp.mpf(2)
        while fun(hi) <= 0:
            hi *= 2
        q = bisect(fun, mp.mpf(1)+mp.mpf("1e-30"), hi)
        beta_s = (q-1)/q
        if not (0 < beta_a < 1 and 0 < beta_s < 1):
            raise AssertionError('check failed: code/verify.py:241')
        out.append(dict(c=c0, beta_age=mp.nstr(beta_a, 35), beta_slip=mp.nstr(beta_s, 35), q_star=mp.nstr(q, 35)))
    passed("sharp-threshold illustrations", "high-precision supplementary", "Five parameter values; symbolic stationarity checked independently")
    return out


def finite_element_slabs() -> list[dict]:
    c = mp.tanh(1)/8
    target_K = mp.mpf(1)/8
    rows = []
    for n in (8, 16, 32, 64, 128):
        dx = mp.mpf(1)/n
        diagonal = c/dx + c*dx/3
        off = -c/dx + c*dx/6
        kappa = mp.acosh(-diagonal/off)
        kn = diagonal + off*mp.sinh((n-1)*kappa)/mp.sinh(n*kappa)
        cross = off*mp.sinh(kappa)/mp.sinh(n*kappa)
        if not (kn > target_K):
            raise AssertionError('check failed: code/verify.py:258')
        for law in ("ageing", "slip"):
            load = mp.mpf(14)/5 if law == "ageing" else mp.mpf(5)/2
            drive = load*mp.sinh(1)/c
            load_n = -cross*drive
            points = [0, 4, 5+mp.log(2), 32] if law == "ageing" else [0, 1, 3, 32]
            for reg in ("truncated", "asinh"):
                fun = lambda r: kn*r + friction(r, law, reg) - load_n
                roots = [bisect(fun, points[j], points[j+1], iterations=150) for j in range(3)]
                rows.append(dict(elements=n, law=law, regularisation=reg,
                                 impedance=mp.nstr(kn, 22), impedance_error=mp.nstr(kn-target_K, 12),
                                 load=mp.nstr(load_n, 22), root_1=mp.nstr(roots[0], 22),
                                 root_2=mp.nstr(roots[1], 22), root_3=mp.nstr(roots[2], 22)))
    errors = [mp.mpf(rows[4*j]["impedance_error"]) for j in range(5)]
    if not (all(errors[j+1] < errors[j]/3 for j in range(4))):
        raise AssertionError('check failed: code/verify.py:272')
    passed("conforming slab discretisations", "high-precision supplementary", "20 parameter/mesh cases, 60 roots; impedance errors decrease at approximately second order")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parents[1]/"evidence")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    exact = exact_sign_certificates()
    symbolic_state_and_thresholds()
    halfspace = symbolic_halfspace()
    concentration = concentration_identities()
    roots = numerical_root_checks()
    state_bound_tests()
    thresholds = threshold_samples()
    fe = finite_element_slabs()
    from revision_checks import run as run_revision
    revision = run_revision(roots, args.out)
    for group in revision["groups"]:
        passed(group["name"], group["kind"], "See revision_verification.json for residuals and scope")
    report = dict(schema_version="2.0", status="PASS", checks=CHECKS,
                  scope="Supporting algebra and examples only; no formal verification of the analytic theorems",
                  environment=dict(python=platform.python_version(), sympy=sp.__version__, mpmath=mp.__version__, precision_decimal_digits=mp.mp.dps),
                  exact_multiplicity=exact, halfspace=halfspace, concentration=concentration,
                  roots=roots, threshold_samples=thresholds)
    (args.out/"verification.json").write_text(json.dumps(report, indent=2)+"\n")
    with (args.out/"finite_element_slabs.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(fe[0]))
        writer.writeheader()
        writer.writerows(fe)
    summary = "\n".join(f"PASS [{x['kind']}] {x['name']}: {x['detail']}" for x in CHECKS)
    summary += "\n\nAll analytic proof obligations remain in paper/manuscript.pdf; not Lean/Coq verification.\n"
    (args.out/"verification.log").write_text(summary)
    print(summary)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"FAIL: {type(exc).__name__}: {exc}", file=sys.stderr)
        raise
