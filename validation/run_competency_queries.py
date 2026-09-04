#!/usr/bin/env python3
from pathlib import Path
import csv
from rdflib import Graph
ROOT=Path(__file__).resolve().parent.parent
g=Graph()
g.parse(ROOT/"ontology/fpoi_ontology_v1.0.1.ttl",format="turtle")
g.parse(ROOT/"validation/fpoi_validation_graph_v1.0.ttl",format="turtle")
rows=[]
for q in sorted((ROOT/"queries").glob("CQ*.rq")):
    result=list(g.query(q.read_text()))
    rows.append((q.stem,len(result)," | ".join(str(v) for v in result[0]) if result else ""))
with open(ROOT/"validation/competency_query_results_reproduced.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["query","rows","values"]); w.writerows(rows)
for r in rows: print(*r,sep="\t")
