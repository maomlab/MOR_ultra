
import requests
from pathlib import Path
import polars as pl
from rdflib import Graph

# Download Guide To Pharmacology data
url = "https://www.guidetopharmacology.org/DATA/rdf/2025.2/gtp-rdf.n3"
iuphar_rdf_path = Path("data/IUPHAR/gtp-rdf.n3")
iuphar_rdf_path.parent.mkdir(parents=True, exist_ok=True)

output_path = Path("intermediate/IUPHAR")
output_path.parent.mkdir(parents=True, exist_ok=True)

# Download data
iuphar_rdf_path.parent.mkdir(parents=True, exist_ok=True)
response = requests.get(url)
response.raise_for_status()  # Raise error for bad status codes
with open(iuphar_rdf_path, "wb") as f:
    f.write(response.content)
print(f"Guide to Pharmacology RDF data downloaded to {iuphar_rdf_path}")



def sparql_results_to_polars(results):
    """
    Convert rdflib SPARQL SELECT query results to a Polars DataFrame.
    """
    # Extract variable names (column headers) from the result
    columns = results.vars
    # Build list of rows
    rows = []
    for row in results:
        row_dict = {str(var): str(row[var]) if row[var] is not None else None for var in columns}
        rows.append(row_dict)
    # Convert to Polars DataFrame
    return pl.DataFrame(rows)


# Load RDF graph
g = Graph()
g.parse(iuphar_rdf_path, format="n3")


# Use Protege to inspect schema and build queries
############ Targets ###########

q = """
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX gtpo: <https://rdf.guidetopharmacology.org/ns/gtpo#>

SELECT ?targetLabel
WHERE {
  ?target gtpo:hasTargetFamily ?targetFamily .
  ?targetFamily gtpo:targetFamilyName \"Opioid receptors\" .

  ?target rdfs:label ?targetLabel
}
"""

results = g.query(q)
df = sparql_results_to_polars(results)
df.write_csv(f"{output_path}/opioid_receptors.tsv", separator='\t')

##### Ligands ##########
q = """
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX gtpo: <https://rdf.guidetopharmacology.org/ns/gtpo#>

SELECT DISTINCT ?ligandLabel ?approved ?canonicalSMILES ?inChI ?xref
WHERE {
	?interaction gtpo:hasTarget ?target .
	?interaction gtpo:hasLigand ?ligand .
	?interaction gtpo:hasAffinity ?affinity .

	?target gtpo:hasTargetFamily ?targetFamily .
	?targetFamily gtpo:targetFamilyName "Opioid receptors" .

	?ligand rdfs:label ?ligandLabel .
	OPTIONAL { ?ligand gtpo:approved ?approved . }
	OPTIONAL { ?ligand gtpo:canonicalSMILES ?canonicalSMILES . }
	OPTIONAL { ?ligand gtpo:inChI ?inChI . }
	OPTIONAL { ?ligand gtpo:xref ?xref }
} 
"""
results = g.query(q)
df = sparql_results_to_polars(results)
df.write_csv(f"{output_path}/opioid_ligands.tsv", separator='\t')


######## Interactions

q = """
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX gtpo: <https://rdf.guidetopharmacology.org/ns/gtpo#>

SELECT ?targetLabel ?targetSpecies ?ligandLabel ?reference ?action ?medianAffinity ?affinityUnits
WHERE {
	?interaction gtpo:hasTarget ?target .
	?interaction gtpo:hasLigand ?ligand .
	?interaction gtpo:hasAffinity ?affinity .

	?target gtpo:hasTargetFamily ?targetFamily .
	?targetFamily gtpo:targetFamilyName "Opioid receptors" .

	?target rdfs:label ?targetLabel .
        OPTIONAL { ?interaction gtpo:hasTaxonomy ?targetSpecies . }
	?ligand rdfs:label ?ligandLabel .

	OPTIONAL { ?interaction gtpo:hasReference ?reference . }
	OPTIONAL { ?interaction gtpo:hasAction ?action . }

	OPTIONAL { ?affinity gtpo:hasMedianValue ?medianAffinity . }
	OPTIONAL { ?affinity gtpo:hasUnits ?affinityUnits }
}
"""
results = g.query(q)
df = sparql_results_to_polars(results)
df.write_csv(f"{output_path}/opioid_interactions.tsv", separator='\t')



g = Graph()
g.parse("/tmp/bao_complete.owl", format="xml")

# URI of the instance
instance_uri = URIRef("http://www.bioassayontology.org/bao/bao_complete.owl#BAO_0190004")

# Print all triples about this term
for s, p, o in g.triples((instance_uri, None, None)):
    print(f"{p} --> {o}")
