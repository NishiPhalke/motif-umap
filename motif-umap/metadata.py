#!/usr/bin/env python3

import sys
import requests

MARKS = [ "H3K4me1", "H3K4me3", "H3K27ac", "CTCF", "WGBS", "H3K9me3", "H3K27me3", "DNase-seq" ]

def reference_epigenome(e):
    j = requests.get("https://www.encodeproject.org/experiments/%s/?format=json" % e).json()["related_datasets"]
    if e == "ENCSR518BPP": j.append(requests.get("https://www.encodeproject.org/experiments/ENCSR017BUL/?format=json").json())
    if e == "ENCSR193SZM": j.append(requests.get("https://www.encodeproject.org/experiments/ENCSR892UWR/?format=json").json())
    return { k: [ x for x in j if x["assay_title"] == k or ("description" in x and k in x["description"]) or ("target" in x and "label" in x["target"] and k == x["target"]["label"]) ][0] for k in MARKS }

def reference_epigenome_signal_files(e):
    e = reference_epigenome(e)
    return { k: (x["accession"], best_signal_file(x["accession"], x["assay_title"] == "WGBS")) for k, x in e.items() if k != "DNase-seq" }

def best_signal_file(e, wgbs = False):
    print(e, file = sys.stderr)
    j = requests.get("https://www.encodeproject.org/experiments/%s/?format=json" % e).json()
    a = [ x for x in reversed(sorted(j["analyses"], key = lambda x: x["date_created"])) ]
    bw = sorted([ x for x in j["files"] if "fold" in x["output_type"] or (wgbs and "signal" in x["output_type"]) ], key = lambda x: -len(x["biological_replicates"]))
    m = max([ len(x["biological_replicates"]) for x in bw ])
    best_bws = { x["@id"] for x in bw if len(x["biological_replicates"]) == m }
    for x in a:
        for f in x["files"]:
            if f in best_bws: return f.split('/')[-2]
