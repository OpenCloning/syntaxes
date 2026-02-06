"""
Add plasmids from the kits to the plasmids directory if they are not already present.
"""

import glob
from opencloning.dna_functions import request_from_addgene
import asyncio
import time

kits = glob.glob("kits/*")
existing_plasmids = glob.glob("addgene_plasmids/*.gb")
existing_addgene_ids = [p.split("/")[-1][:-3] for p in existing_plasmids]

for kit in kits:
    plasmids = glob.glob(f"{kit}/plasmids.tsv")
    if len(plasmids) != 1:
        raise Exception(f"Expected 1 plasmid file in {kit}, but found {len(plasmids)}")
    plasmids_file = plasmids[0]
    with open(plasmids_file, "r") as f:
        plasmids_lines = f.readlines()
    addgene_ids = [line.split("\t")[2] for line in plasmids_lines[1:]]

    requests_made = 0
    for addgene_id in addgene_ids:
        if addgene_id in existing_addgene_ids:
            continue
        print(f"Requesting {addgene_id}")
        seq = asyncio.run(request_from_addgene(addgene_id))
        requests_made += 1
        with open(f"addgene_plasmids/{addgene_id}.gb", "w") as f:
            f.write(seq.format("genbank"))

        if requests_made % 10 == 0:
            time.sleep(2)
