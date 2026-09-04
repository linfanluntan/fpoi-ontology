#!/usr/bin/env python3
from pathlib import Path
from rdflib import Graph,Namespace,RDF,OWL
import json
ROOT=Path(__file__).resolve().parent.parent
SH=Namespace("http://www.w3.org/ns/shacl#")
g=Graph().parse(ROOT/"ontology/fpoi_ontology_v1.0.1.ttl",format="turtle")
s=Graph().parse(ROOT/"shapes/fpoi_shapes_v1.0.1.ttl",format="turtle")
v=Graph().parse(ROOT/"validation/fpoi_validation_graph_v1.0.ttl",format="turtle")
print(json.dumps({
"ontology_triples":len(g),
"classes":len(set(g.subjects(RDF.type,OWL.Class))),
"object_properties":len(set(g.subjects(RDF.type,OWL.ObjectProperty))),
"datatype_properties":len(set(g.subjects(RDF.type,OWL.DatatypeProperty))),
"disjointwith_axioms":len(list(g.triples((None,OWL.disjointWith,None)))),
"shacl_node_shapes":len(set(s.subjects(RDF.type,SH.NodeShape))),
"shacl_property_constraints":len(list(s.triples((None,SH.property,None)))),
"validation_graph_triples":len(v)},indent=2))
