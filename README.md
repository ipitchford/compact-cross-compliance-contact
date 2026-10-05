# Cross-compliance compactness in three-dimensional frictional contact
## Revised Evidence Press candidate v2.1 — 5 October 2026

DOI: https://doi.org/10.5281/zenodo.23171696 . See [REVIEW_RESPONSE_V21.md](REVIEW_RESPONSE_V21.md) for the latest review dispositions and [AI_INDEX.md](AI_INDEX.md) for the agent entry point. This publication supersedes the supplied local v2.0 draft, not an earlier public DOI version. Supporting checks also run under Python -O; six corrupted-check invocations are rejected by `python code/test_negative.py`.

**Manuscript:** `paper/manuscript.pdf` (32 pages), editable source `paper/manuscript.tex`.
**Current review response:** `REVIEW_RESPONSE_V21.md`. The earlier v2.0 response is retained in `response/response.pdf` and `response/RESPONSE_TO_REVIEW.md` as revision history, not a substitute for the current dispositions.
**Release status:** revised research candidate, not a publication acceptance, independently certified theorem, or awarded REF rating.

### What changed
The revision addresses all four major and five minor referee comments. The new mathematical core is an actual three-dimensional incremental Signorini–Coulomb existence theorem for smooth anchored isotropic bodies with a closed matched interface. Here the normal pressure is an unknown reaction. Matching means equality of lambda+mu across the interface, not identity of the materials. The proof uses compactness of the normal–tangential elastic cross-compliance block, type-(S)+ coercivity, and the established pseudomonotone variational-inequality theorem.

New consequences include strong solution stability under data perturbations, convergent reaction-space approximations, and an if-and-only-if classification of the unrestricted dominated weak-product property in this smooth two-body class. A finite-layer calculation demonstrates nonzero lower-frequency coupling despite high-frequency cancellation.

The logarithmic-Euler prior-art comparison now includes Gräser, Kornhuber and Podlesny (2024). The state integrator is not claimed new. The prescribed-pressure existence theorem is presented as a self-contained proof with specified extensions. The new energy table distinguishes stationary roots, restricted local minimality, global minimisers and dynamic stability.

### Reproduce
Python 3.13.5 was used. Install the pinned requirements in an isolated environment if needed. The two main modules require SymPy and mpmath; the unchanged reviewer script also requires NumPy and SciPy. From this directory:

```sh
python code/verify.py --out /tmp/friction_v2_replay
python code/check_manifest.py
python provenance/reviewer_independent_checks.py --out /tmp/reviewer_replay.json
```

The main command creates its output directory and runs 12 check groups. Compare the generated `verification.json`, `verification.log`, `revision_verification.json`, `revision_verification.log`, `finite_element_slabs.csv`, and `matched_layer_spectrum.csv` against `evidence/`. No network access or proprietary solver is required. The reviewer script requires the parent output directory to exist.

To build the manuscript, with a TeX installation containing newtx, amsmath, amsthm, mathtools, microtype, booktabs, enumitem, fancyhdr, hyperref and longtable:

```sh
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=paper paper/manuscript.tex
```

Numerical JSON outputs are deterministic in the recorded environment. PDF hashes may change on recompilation because of TeX timestamps; the manifest checks the distributed files, not byte-identical recompilation.

### Evidence boundaries
Exact algebra and rational sign certificates are separate from numerical root, energy and spectral evaluations. None formalises the infinite-dimensional proofs. The supplied reviewer audit was AI-assisted and reviewed v1.0; replaying it does not review the new v2.0 contact theorems. No human specialist endorsement or exhaustive novelty search is claimed.

The existence result does not cover arbitrary mismatched interfaces, contact edges, rough coefficients, evolving geometry, or the multidimensional continuous-time rate-and-state limit. The weak-product counterexample does not prove contact nonexistence. A failure of the matching condition does not establish failure of contact existence. General Problem 4 is not claimed solved.

### Package guide
`notes/` contains the claim ledger, proof audit, source comparison, falsifier ledger, and historical research gates. `provenance/` preserves supplied audit code/results and input hashes; private review prose is excluded. `evidence/` contains supporting checks and replay comparisons. `STATUS.json` and `REVIEW_RESPONSE_V21.md` describe the current revision; `response/REVIEW_ACTIONS.json` records the earlier revision. Root `RESEARCH_GATES.json` is the current publication-intake record.
