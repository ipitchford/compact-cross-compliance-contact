# Agent entry point — cross-compliance contact

Anonymous; v2.1.0-candidate; unrefereed. DOI: 10.5281/zenodo.23171696.

## Claim and scope
[paper/manuscript.tex](paper/manuscript.tex) and [PDF](paper/manuscript.pdf) contain numbered written proofs. Theorem 4.1: compact cross compliance yields unknown-pressure incremental Coulomb existence without small friction. Theorem 5.2 realises this for smooth anchored bodies with a closed interface and equal lambda+mu across it. Theorem 5.3 classifies a dominated weak-product property, NOT contact existence. Sections 7–13 concern separate prescribed-pressure rate-and-state models.

## Evidence and replay
[code/verify.py](code/verify.py) runs twelve supporting groups using the pinned [requirements.txt](requirements.txt). Run `python code/verify.py --out /tmp/friction-replay` in a fresh copy. Repeat with `python -O`; explicit exceptions preserve checks. [code/test_negative.py](code/test_negative.py) rejects damaged rational signs, material matching and a low-frequency coefficient in both modes. [code/check_manifest.py](code/check_manifest.py) checks distributed hashes. [evidence/verification.json](evidence/verification.json) and [evidence/revision_verification.json](evidence/revision_verification.json) are retained input evidence; new replay does not formalise the continuum proof.

## Dependencies and limitations
[notes/PROOF_AUDIT.md](notes/PROOF_AUDIT.md) maps analytic dependencies; [notes/CLAIM_LEDGER.md](notes/CLAIM_LEDGER.md) separates claims. [ASSURANCE.md](ASSURANCE.md) fixes the assurance ceiling. [notes/LITERATURE_AUDIT.md](notes/LITERATURE_AUDIT.md) records the earlier bounded source audit; [REVIEW_RESPONSE_V21.md](REVIEW_RESPONSE_V21.md) records current corrections. No exhaustive novelty conclusion is available. Unknown-pressure and prescribed-pressure models must not be merged.

## Provenance and reuse
[PROVENANCE.md](PROVENANCE.md), [LICENSES.md](LICENSES.md), and [MANIFEST.sha256](MANIFEST.sha256) describe inputs, rights and exact bytes. Original prose/data CC0; original code MIT. Private review text is excluded from the publication archive. Existing imported reviewer code retains NOASSERTION rather than being relabelled as original. New work should target independent analytic scrutiny, mismatched interfaces or explicitly consistent discretisation, not infer these from the finite replay.
