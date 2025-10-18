


import requests
import pandas as pd
import time
import yaml

with open("parameters.yaml", 'r') as file:
    parameters = yaml.safe_load(file)

KeBOB_referances = pd.read_csv(
    "intermediate/01.1_gather_KeBOB/KeBOB_references.tsv",
    sep = "\t")

papers = []

for index, reference in KeBOB_references.iterrows():
    print(f"Downloading paper for reference '{reference["reference"]}'")
    
    r = requests.get(
        f"https://api.unpaywall.org/v2/{reference["DOI"]}",
        params={"email": parameters["web_tokens"]["user_agent"]})
    
    if r.ok:
        data = r.json()
        oa = data.get("best_oa_location", {})
        papers.append({
            "doi": doi,
            "title": data.get("title"),
            "pdf_url": oa.get("url_for_pdf")
        })
    time.sleep(0.1)  # be polite

df = pd.DataFrame(papers)
print(df)
