# Proof and implementation audit — v2.0

## Corpus and status
The starting objects were the exact v1.0 package, its PDF/source and eight-group verifier, the supplied 430-line review, and the seven-file referee audit archive. The original 17-file manifest and referee seven-file manifest both matched. The unchanged referee script was replayed; its JSON was byte-identical to the supplied result. This is regression evidence for v1.0, not an independent audit of new v2.0 proofs.

## New contact proof: load-bearing checks

1. **Dual pairings:** negative H^{-1/2} distributions are measures. Absolute values of H^{1/2} vector fields remain in H^{1/2}, continuously and with a norm bound. All pressure products are positive dual pairings, interpreted by appropriate Sobolev representatives, not arbitrary point values.
2. **Frozen minimiser:** symmetric coercive A and a finite continuous convex nonnegative friction term give a unique minimiser. Comparing with zero bounds it on bounded reaction/data sets, including nonzero old slip.
3. **Traction reconstruction:** the Euler residual annihilates the tangential-trace kernel. Joint trace surjectivity yields the tangential reaction; testing arbitrary traces proves domination and the exact friction-work equality. It is not postulated independently.
4. **Normal strong monotonicity:** opposite frozen inequalities control displacement difference. A lift with prescribed normal trace and zero tangent trace controls the normal reaction difference. This is the crucial quantitative estimate, independent of the frozen friction weight.
5. **Compactness:** for fixed normal load tau, the variable frozen coefficient enters the normal displacement only through C_nt times a bounded tangential reaction. Compactness is assumed at the correct endpoint H^{-1/2} to H^{1/2}; ordinary trace compactness into L2 would not suffice.
6. **Type-(S)+:** strong precompactness of S(t_j,t), paired with weak-null t_j-t, plus the normal strong-monotonicity estimate, forces norm convergence. Strong continuity then gives pseudomonotonicity.
7. **Coercivity:** pure-normal lifts bound ||t|| by ||u||+||ell||. Friction work bounds <F(t),t> below by elastic energy minus load and an old-slip term linear in ||t||. Young's inequality supplies positive quadratic reaction growth.
8. **Contact existence:** the established Brézis theorem, as stated in Ballard–Iurlano Appendix A, applies on the nonpositive cone. Tests 0, 2t and t minus a positive test function recover complementarity and nonpenetration. The tangential Euler inequality supplies actual Coulomb friction. No pressure regularisation was introduced.
9. **Weighted variant:** domination permits s=f sigma with |sigma|≤-t. The same H^{-1/2} bound applies to sigma even where f vanishes, after setting it to zero there. Compact C_nt M_f suffices.
10. **Data stability:** variation of old-slip friction functionals is uniformly bounded by <-ft,|w_j-w|>. Strong convexity controls frozen minimiser changes. Testing the reaction VI by the weak limit yields the type-(S)+ criterion. No uniqueness or full-sequence convergence is silently inserted.
11. **Reaction approximation:** finite-dimensional projected VIs exist by coercivity and Brouwer. Strong cone-recovery tests yield the same type-(S)+ argument. This is not an implemented finite-element contact solver.

## Smooth bounded bodies and localised obstruction
The relative trace is onto because the shared component and clamps are separated. Korn coercivity holds on each anchored body. The contact Dirichlet-to-Neumann map is the contact-to-contact block of the full smooth-boundary map with zero remote Dirichlet data. Its inverse exists by coercive mixed elasticity and is an order-minus-one pseudodifferential operator. The principal symbol is checked against Tan–Liu Theorem 1.1 and the exact half-space calculation.

Opposite normals and opposite relative forces produce the DIFFERENCE of chi=1/[2(lambda+mu)]. Equality of lambda+mu cancels the entire principal normal–tangential block. The remainder is order at most minus two, hence maps H^{-1/2} to H^{3/2}, compactly into H^{1/2}. Curvature, smooth coefficient derivatives and remote clamping may retain lower-order coupling.

For mismatch, concentrating smooth radial profiles are localised in a geodesic chart with density/Jacobian factors shared between pressure and traction, preserving domination. The leading kernel frozen at the point yields c(x0)pi/8. The variable leading-kernel error is bounded by omega(C epsilon)/|x-y| and therefore by the bounded positive Riesz energy times o(1). The order-minus-two remainder is compact. No claim that these fields solve the nonlinear contact equations is made.

## Concrete finite-layer check
Both growing and decaying nilpotent longitudinal modes satisfy the Navier equations symbolically. Their traction matrices and the top-clamp reflection identity are verified. The exact compliance inverse gives high-frequency chi/k plus exponential corrections. The low-frequency boundary ODE independently yields (mu-lambda)H^2/[2mu(lambda+2mu)] times k. For (2,1) versus (1,2) this is -1/8-1/20=-7/40, so cancellation is not exact decoupling. The periodic-layer result is a second, elementary realisation of the abstract hypothesis, not a proof of curved-boundary symbol calculus by simulation.

## Referee-mandated interpretive checks
The GKP logarithmic state shift is exact for both laws. Its integrator is not new. The thesis uses an exact flow; the cited earlier Theorem 5.14 concerns ageing, not slip. The revised comparison is limited to the inspected statements. Numerically ranked energies describe only the displayed roots; the middle restricted branch is a negative-curvature stationary point. Positive restricted curvature cannot establish full-space local minimality. Global-minimiser Galerkin convergence and dynamic stability remain distinct.

The characteristic theorem now fixes material speeds and states its input and output topologies explicitly. The JSON names Cartesian u3 and outward un=-u3 separately.

## Check failure found and corrected during this revision
An initial spectral assertion failed because it compared a 40-significant-digit serialised approximation to -1/12 with an absolute tolerance of 10^{-45}. The calculation itself ran at 100 digits. The stored spectral values were increased to 70 significant digits, retaining the 100-digit calculation and the original tolerance. All exact identities and the subsequent twelve-group replay passed. No theorem or tolerance was weakened to suppress that representation error.

A notation collision was also removed: H_plus/H_minus denoted traction matrices and had briefly been reused for layer thicknesses. The final example says “unit layer thicknesses”.

## Remaining audit gates
A specialist must still inspect the analytic proof, the precise pseudodifferential boundary calculus application, and originality against older contact literature. The new theorems have not been independently reviewed. Formal verification is absent. The stronger theorem does not remove the general continuous-time boundary-trace gap or solve mismatched/rough/contact-edge configurations.
