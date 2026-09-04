# Failure Pattern Operational Intelligence (FPOI)

FPOI is a versioned application ontology and operational knowledge-graph framework for representing radiation-oncology failure-pattern evidence, provenance, QA state, uncertainty, and governed interpretation.

## Release state
This repository package is **GitHub-ready but not yet the public release**. Before manuscript submission:
- publish this repository and replace `__GITHUB_OWNER__`;
- commit/tag the authoritative PROTEAS gate configuration (closes C9);
- complete C4 published-rule comparison;
- execute and commit C8 DICOM-RT ablation results;
- create a Zenodo release/DOI;
- register `https://w3id.org/fpoi/` and verify dereferenceability.

## Verified ontology evaluation
FPOI v1.0.1 currently contains 126 ontology triples, 22 OWL classes, 19 object properties, 15 datatype properties, two explicit disjointness axioms, two SHACL node shapes, and eight `sh:property` constraints (7 top-level + 1 nested qualified constraint). The compact validation graph has 46 triples. All seven competency queries return the expected state. The unmodified graph conforms under pySHACL 0.40.1 and four malformed controls are rejected.

## Quick validation
```bash
python -m pip install -r requirements.txt
python validation/check_artifact_stats.py
python validation/run_competency_queries.py
python validation/run_pyshacl_validation.py
```

## Data
Patient imaging is not redistributed. Obtain PROTEAS and OpenKBP from their original sources and terms.

## Licensing
See `LICENSE.md`.
