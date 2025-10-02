"""SEARCH HYPATIA"""

#%% 
"""For searching elements in the Hypatia Catalog"""
import requests
params = {"name": ['hd 11695','hd 18884','hd 18907','hd 22049','hd 23249',
    'hd 128621','hd 10700','hd 100407'], "element": ["Eu_II"], 
          "solarnorm": ["grevesse07"]}
entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/star", params=params)
print(entry.json())

#%%
"""For finding the solar normilisation values"""
import requests
entry = requests.get("https://hypatiacatalog.com/hypatia/api/v2/solarnorm")
print(entry.json())
# %%
"""From GPT-4: For finding all the elements available in the Hypatia Catalog"""
#!/usr/bin/env python3
"""
Query Hypatia Catalog API for Ba II and Eu II abundances and list source catalog citations.

Saves results to hypatia_hyper_sources.csv
"""
import requests
import pandas as pd
from time import sleep

API_BASE = "https://hypatiacatalog.com/hypatia/api/v2"

# stars you gave
stars = [
    'hd 11695','hd 18884','hd 18907','hd 22049','hd 23249',
    'hd 128621','hd 10700','hd 100407'
]

# elements: use Hypatia element naming for ions (Ba II, Eu II)
elements = ["Ba_II", "Eu_II"]

# 1) fetch catalog list and build mapping id -> author (year) + DOI
cat_resp = requests.get(f"{API_BASE}/catalog")
cat_resp.raise_for_status()
catalog_list = cat_resp.json()
catalog_map = {
    c["id"]: {
        "citation": f"{c.get('author','')} ({c.get('year','')})".strip(),
        "doi": c.get("doi")
    }
    for c in catalog_list
}

# container for rows
rows = []

# Helper to call composition endpoint for one star/element
def query_composition(star_name, element, solarnorm="grevesse07"):
    url = f"{API_BASE}/composition"
    params = {"name": [star_name], "element": [element]}
    if solarnorm:
        params["solarnorm"] = [solarnorm]
    resp = requests.get(url, params=params)
    resp.raise_for_status()
    data = resp.json()
    if not data:
        return None
    return data[0]  # first (and only) group

# iterate stars/elements
for star in stars:
    for el in elements:
        try:
            comp = query_composition(star, el)
            sleep(0.1)
        except Exception as e:
            print(f"Failed for {star} {el}: {e}")
            continue

        if comp is None:
            rows.append({
                "requested_name": star,
                "match_name": None,
                "element": el,
                "median_value": None,
                "median_type": "no_data",
                "catalog_id": None,
                "catalog_citation": None,
                "catalog_doi": None,
                "reported_value": None,
                "note": "no composition returned"
            })
            continue

        median_value = comp.get("median_value", None)
        mean_value = comp.get("mean", None)

        # summary row
        rows.append({
            "requested_name": comp.get("requested_name"),
            "match_name": comp.get("match_name"),
            "element": comp.get("element"),
            "median_value": median_value,
            "mean_value": mean_value,
            "reported_value": None,
            "catalog_id": None,
            "catalog_citation": None,
            "catalog_doi": None,
            "note": "SUMMARY (median/mean)"
        })

        # individual values
        all_vals = comp.get("all_values", [])
        for av in all_vals:
            value = av.get("value")
            catalog_obj = av.get("catalog", {}) or {}
            cat_id = catalog_obj.get("id") if isinstance(catalog_obj, dict) else catalog_obj
            if cat_id is None:
                cat_id = av.get("catalog_id") or av.get("source")

            cat_info = catalog_map.get(cat_id, {})
            cat_cite = cat_info.get("citation", "")
            cat_doi = cat_info.get("doi")

            spread = comp.get("spread")
            std_dev = comp.get("std_dev")

            if spread is not None:
                median_uncertainty = spread / 2
            elif std_dev is not None:
                median_uncertainty = std_dev
            elif len(all_vals) == 1:
                median_uncertainty = all_vals[0].get("err")
            else:
                median_uncertainty = None

            rows.append({
                "requested_name": comp.get("requested_name"),
                "match_name": comp.get("match_name"),
                "element": comp.get("element"),
                "solarnorm": comp.get("solarnorm"),
                "median_value": None,
                "median_uncertainty": median_uncertainty,
                "mean_value": None,
                "std_dev": comp.get("std_dev"),
                "reported_value": value,
                "catalog_id": cat_id,
                "catalog_citation": cat_cite,
                "catalog_doi": cat_doi,
                "note": "individual measurement"
            })

import pandas as pd

# Create the DataFrame (your code)
df = pd.DataFrame(rows)

# --- NEW PART: collect unique bibcodes ---
# Some rows may have None, so drop them
bibcodes = df["catalog_citation"].dropna().unique()

# Convert to a string you can paste into ADS
bibcode_query = " OR ".join([f"bibcode:{b}" for b in bibcodes])

print("Paste this into ADS search bar:\n")
print(bibcode_query)


# collect unique citations and DOIs per star/element
citation_summary = (
    df[df["catalog_citation"].notna()]
    .groupby(["match_name", "element"])
    .agg({
        "catalog_citation": lambda x: "; ".join(sorted(set(filter(None, x)))),
        "catalog_doi": lambda x: "; ".join(sorted(set(filter(None, x))))
    })
    .reset_index()
    .rename(columns={
        "catalog_citation": "unique_citations",
        "catalog_doi": "unique_dois"
    })
)

# merge back
df = df.merge(citation_summary, on=["match_name", "element"], how="left")

summary_df = df[df["note"] == "SUMMARY (median/mean)"].copy()
summary_df = summary_df[[
    "requested_name", "match_name", "element",
    "median_value", "mean_value", "unique_citations", "unique_dois"
]]

summary_df.to_csv("/home/users/qai11/Documents/Isotope-Pipeline/hypatia_ba_eu_summary.csv", index=False)
# print(f"Wrote summary with citations and DOIs to hypatia_ba_eu_summary.csv")

out_csv = "/home/users/qai11/Documents/Isotope-Pipeline/hypatia_ba_eu_sources.csv"
df.to_csv(out_csv, index=False)
# print(f"Wrote {len(df)} rows to {out_csv}")
# print(df.head(40).to_string(index=False))



# %%
