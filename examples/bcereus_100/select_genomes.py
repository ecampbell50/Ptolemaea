#!/usr/bin/env python3
"""
select_genomes.py
Build the fixed example genome list (accessions.tsv) used by run_example.sh.

You do not need to run this - accessions.tsv is already in the repo so the example
is reproducible. It is kept to document exactly how the 100 genomes were chosen.

Method:
  1. Query NCBI Datasets for every complete, non-atypical RefSeq genome in the
     Bacillus cereus group (NCBI taxon 86661).
  2. Group by species and sort each species' genomes by accession.
  3. Take genomes round-robin across species (one from each species, then a second
     from each, ...) until N are selected, so rare species are represented and the
     set is not dominated by B. anthracis...

Standard library only.

Usage:
    python3 select_genomes.py [--n 100] [--output accessions.tsv]
"""

import argparse
import json
import urllib.request
from collections import defaultdict
from datetime import date

TAXON = 86661  # Bacillus cereus group
API = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/{taxon}/dataset_report"
       "?filters.assembly_source=refseq"
       "&filters.assembly_level=complete_genome"
       "&filters.exclude_atypical=true"
       "&page_size=1000")


def fetch_reports(taxon):
    """Return every genome report for the taxon, following page tokens."""
    reports, token = [], None
    while True:
        url = API.format(taxon=taxon) + (f"&page_token={token}" if token else "")
        with urllib.request.urlopen(url) as resp:
            data = json.load(resp)
        reports.extend(data.get("reports", []))
        token = data.get("next_page_token")
        if not token:
            return reports


def species_of(report):
    """'Bacillus cereus ATCC 14579' -> 'Bacillus cereus' (brackets stripped)."""
    words = report["organism"]["organism_name"].replace("[", "").replace("]", "").split()
    return " ".join(words[:2])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--n", type=int, default=100, help="number of genomes (default 100)")
    parser.add_argument("--output", default="accessions.tsv", help="output TSV")
    args = parser.parse_args()

    reports = fetch_reports(TAXON)
    print(f"{len(reports)} candidate genomes from NCBI")

    by_species = defaultdict(list)
    for r in reports:
        by_species[species_of(r)].append(r)
    for genomes in by_species.values():
        genomes.sort(key=lambda r: r["accession"])

    # Round-robin across species (alphabetical order for determinism)
    selected, depth = [], 0
    while len(selected) < args.n:
        added = False
        for sp in sorted(by_species):
            if depth < len(by_species[sp]) and len(selected) < args.n:
                selected.append(by_species[sp][depth])
                added = True
        if not added:
            break
        depth += 1

    with open(args.output, "w") as out:
        out.write(f"# {len(selected)} complete RefSeq B. cereus group genomes, "
                  f"selected with select_genomes.py on {date.today().isoformat()}\n")
        out.write("accession\tspecies\torganism_name\n")
        for r in sorted(selected, key=lambda r: (species_of(r), r["accession"])):
            out.write(f"{r['accession']}\t{species_of(r)}\t{r['organism']['organism_name']}\n")

    print(f"Wrote {len(selected)} genomes from {len({species_of(r) for r in selected})} "
          f"species to {args.output}")


if __name__ == "__main__":
    main()
