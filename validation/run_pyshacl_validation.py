#!/usr/bin/env python3
"""Validate FPOI positive and negative SHACL controls with pySHACL."""
from pathlib import Path
import csv, copy
import pyshacl
from pyshacl import validate
from rdflib import Graph, Namespace, RDF, Literal

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "fpoi_validation_graph_v1.0.ttl"
SHAPES = ROOT.parent / "shapes" / "fpoi_shapes_v1.0.1.ttl"
ONT = ROOT.parent / "ontology" / "fpoi_ontology_v1.0.1.ttl"
FPOI = Namespace("https://w3id.org/fpoi/")

def run(g):
    conforms, report_graph, report_text = validate(
        data_graph=g,
        shacl_graph=str(SHAPES),
        ont_graph=str(ONT),
        inference="none",
        abort_on_first=False,
        allow_infos=False,
        allow_warnings=False,
        meta_shacl=True,
        advanced=False,
        debug=False,
    )
    SH = Namespace("http://www.w3.org/ns/shacl#")
    n_results = len(set(report_graph.subjects(RDF.type, SH.ValidationResult)))
    return bool(conforms), int(n_results)

base = Graph().parse(DATA, format="turtle")
rows=[]

c,n=run(base)
rows.append(["positive_unmodified","CONFORM","CONFORM" if c else "REJECT",n])

# 1 missing target correspondence
g=Graph()
for t in base: g.add(t)
g.remove((FPOI.fp_P23b_FU1,FPOI.supportedByTargetCorrespondence,None))
c,n=run(g)
rows.append(["missing_target_correspondence","REJECT","CONFORM" if c else "REJECT",n])

# 2 nonpassing target correspondence
g=Graph()
for t in base: g.add(t)
g.remove((FPOI.tca_P23b,FPOI.verdict,None))
g.add((FPOI.tca_P23b,FPOI.verdict,Literal("FAIL")))
c,n=run(g)
rows.append(["nonpassing_target_correspondence","REJECT","CONFORM" if c else "REJECT",n])

# 3 missing classification scheme
g=Graph()
for t in base: g.add(t)
g.remove((FPOI.fp_P23b_FU1,FPOI.classifiedUnder,None))
c,n=run(g)
rows.append(["missing_classification_scheme","REJECT","CONFORM" if c else "REJECT",n])

# 4 missing computation run
g=Graph()
for t in base: g.add(t)
g.remove((FPOI.fp_P23b_FU1,FPOI.generatedBy,None))
c,n=run(g)
rows.append(["missing_computation_run","REJECT","CONFORM" if c else "REJECT",n])

with open(ROOT/"shacl_validation_results.csv","w",newline="",encoding="utf-8") as f:
    w=csv.writer(f)
    w.writerow(["test","expected","observed","validation_results"])
    w.writerows(rows)

(ROOT/"shacl_engine_version.txt").write_text(
    f"pySHACL {pyshacl.__version__}\n", encoding="utf-8"
)
print("pySHACL", pyshacl.__version__)
for r in rows: print(*r)
