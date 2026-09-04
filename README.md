# Failure Pattern Operational Intelligence (FPOI)

FPOI is a versioned application ontology and operational knowledge-graph framework for representing radiation-oncology failure-pattern evidence, provenance, QA state, uncertainty, and governed interpretation.

## Repository

This is the public development repository for FPOI: `https://github.com/linfanluntan/fpoi-ontology`.

The ontology/semantic release is v1.0.1. The manuscript submission release remains **pre-submission** until the remaining release gates below are closed.

## Verified ontology evaluation

FPOI v1.0.1 contains 126 ontology triples, 22 OWL classes, 19 object properties, 15 datatype properties, two explicit disjointness axioms, two SHACL node shapes, and eight `sh:property` constraints (7 top-level + 1 nested qualified constraint). The compact validation graph has 46 triples. All seven competency queries return the expected state. The unmodified graph conforms under pySHACL 0.40.1 and four malformed controls are rejected.

The corrected SHACL semantics intentionally permit a non-passing `TreatmentTargetCorrespondenceAssessment` to remain valid evidence. A normalized `FailurePattern` instead requires at least one supporting target-correspondence assessment whose verdict is `PASS`.

## Quick validation

```bash
python -m pip install -r requirements.txt
python validation/check_artifact_stats.py
python validation/run_competency_queries.py
python validation/run_pyshacl_validation.py
```

## Pre-submission release gates

- C4: complete and commit the published-rule comparison for the OpenKBP classification experiment.
- C8: execute and commit systematic DICOM-RT object-ablation negative controls.
- C9: merge the authoritative PROTEAS gate configuration into `main` through this release-preparation change set and cite the resulting commit SHA in the manuscript/release provenance.
- Archive a tagged release in Zenodo and add the DOI to `CITATION.cff` and the manuscript.
- Register `https://w3id.org/fpoi/` and verify dereferenceability before claiming the persistent namespace is active.

## Data

Patient imaging is not redistributed. Obtain PROTEAS and OpenKBP from their original sources and terms. The five-package DICOM-RT test set is not redistributed because its provenance and redistribution status were not independently established.

## Licensing

See `LICENSE.md`.
