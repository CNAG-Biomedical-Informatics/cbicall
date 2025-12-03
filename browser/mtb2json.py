#!/usr/bin/env python3
"""
mtb2json.py - Parse MToolBox prioritized variants output

Usage:
  mtb2json.py -i mit_prioritized_variants.txt -f json4html > mtDNA.json
"""

import argparse
import csv
import json
import sys
from pprint import pprint


KEYS2REPORT = [
    "Variant_Allele",
    "Sample",
    "Locus",
    "Nt_Variability",
    "Codon_Position",
    "Aa_Change",
    "Aa_Variability",
    "REF",
    "ALT",
    "GT",
    "DP",
    "HF",
    "tRNA_Annotation",
    "Disease_Score",
    "RNA_predictions",
    "Mitomap_Associated_Disease(s)",
    "Mitomap_Homoplasmy",
    "Mitomap_Heteroplasmy",
    "Somatic_Mutations",
    "SM_Homoplasmy",
    "SM_Heteroplasmy",
    "ClinVar",
    "OMIM_link",
    "dbSNP_ID",
    "Mamit-tRNA_link",
    "PhastCons20Way",
    "PhyloP20Way",
    "AC/AN_1000_Genomes",
    "1000_Genomes_Homoplasmy",
    "1000_Genomes_Heteroplasmy",
]

KEYS4HTML = [
    "Sample",
    "Locus",
    "Variant_Allele",
    "REF",
    "ALT",
    "Aa_Change",
    "GT",
    "DP",
    "HF",
    "tRNA_Annotation",
    "Disease_Score",
    "RNA_predictions",
    "Mitomap_Associated_Disease(s)",
    "Mitomap_Homoplasmy",
    "Mitomap_Heteroplasmy",
    "ClinVar",
    "OMIM_link",
    "dbSNP_ID",
    "Mamit-tRNA_link",
    "AC/AN_1000_Genomes",
    "1000_Genomes_Homoplasmy",
    "1000_Genomes_Heteroplasmy",
]


def parse_args():
    p = argparse.ArgumentParser(
        description="Parse MToolBox prioritized variants and output hash/json/json4html/tsv"
    )
    p.add_argument("-i", "--input", required=True, help="MToolBox prioritized_variants.txt")
    p.add_argument(
        "-f",
        "--format",
        default="hash",
        choices=["hash", "json", "json4html", "tsv"],
        help="Output format",
    )
    p.add_argument(
        "--HF",
        type=float,
        default=0.30,  # default HF cutoff
        help="Heteroplasmic Fraction cutoff for 01P [0.30]",
    )
    p.add_argument(
        "--MAF",
        type=float,
        default=0.01,  # default MAF cutoff
        help="MAF filter [0.01]; keep variant if 1000G MAF < MAF",
    )
    p.add_argument(
        "--denovo",
        action="store_true",
        help="Denovo mutations (reserved, not implemented)",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Print debugging info to stderr",
    )
    return p.parse_args()


def normalize_header(header_fields):
    """
    Convert header names like 'Variant Allele' to 'Variant_Allele'.
    """
    normalized = []
    for h in header_fields:
        h = h.strip().replace(" ", "_")
        normalized.append(h)
    return normalized


def sort_sample_alphabetically(sample_str):
    parts = [s for s in sample_str.split(",") if s]
    if not parts:
        return sample_str
    parts.sort()
    return ",".join(parts)


def passes_maf_filter(row, col_idx, maf_cutoff):
    if maf_cutoff <= 0:
        return True
    idx = col_idx.get("AC/AN_1000_Genomes")
    if idx is None:
        return True
    val = row[idx].strip()
    if not val:
        return True
    try:
        if float(val) >= maf_cutoff:
            return False
    except ValueError:
        return True
    return True


def max_hf_for_proband(hf_str):
    """
    HF patterns:
      single mode => "0.02,0.08"
      cohort mode => "01P:0.02,0.8|02M:N/A|03F:1.0"
    Only consider the first block (01P), turn N/A to 0, take max.
    """
    if not hf_str:
        return 0.0

    first_block = hf_str.split("|", 1)[0]
    first_block = first_block.replace("01P:", "")
    first_block = first_block.replace("N/A", "0")

    parts = [x.strip() for x in first_block.split(",") if x.strip()]
    if not parts:
        return 0.0

    vals = []
    for x in parts:
        try:
            vals.append(float(x))
        except ValueError:
            continue

    return max(vals) if vals else 0.0


def passes_hf_filter(row, col_idx, hf_cutoff):
    if hf_cutoff <= 0:
        return True
    idx = col_idx.get("HF")
    if idx is None:
        return True
    hf_str = row[idx].strip()
    m = max_hf_for_proband(hf_str)
    return m > hf_cutoff


def build_hash_out(rows, header, hf_cutoff, maf_cutoff, debug=False):
    col_idx = {name: i for i, name in enumerate(header)}
    hash_out = {}

    for line_no, row in enumerate(rows, start=2):
        if len(row) < len(header):
            if debug:
                print(f"Skipping short line {line_no}", file=sys.stderr)
            continue

        # Discard synonymous variants
        aa_idx = col_idx.get("Aa_Change")
        if aa_idx is not None and row[aa_idx].strip() == "syn":
            continue

        if not passes_maf_filter(row, col_idx, maf_cutoff):
            continue

        if not passes_hf_filter(row, col_idx, hf_cutoff):
            continue

        # Normalize sample
        sample_idx = col_idx.get("Sample")
        if sample_idx is not None:
            row[sample_idx] = sort_sample_alphabetically(row[sample_idx])

        variant_idx = col_idx.get("Variant_Allele")
        locus_idx = col_idx.get("Locus")
        if variant_idx is None or locus_idx is None or sample_idx is None:
            continue

        index_key = "{}_{}_{}".format(
            row[sample_idx].strip(),
            row[locus_idx].strip(),
            row[variant_idx].strip(),
        )

        record = {}
        for key in KEYS2REPORT:
            idx = col_idx.get(key)
            if idx is None or idx >= len(row):
                record[key] = ""
            else:
                record[key] = row[idx].strip()

        hash_out[index_key] = record

    return hash_out


def hash2array(hash_out):
    """
    Build AoA for DataTables, with HTML links for Locus, dbSNP_ID, OMIM_link,
    but NO color logic (that lives in DataTables JS).
    """
    rows = []

    for key1 in sorted(hash_out.keys()):
        record = hash_out[key1]
        row = []

        for key2 in KEYS4HTML:
            tmp = record.get(key2, "")

            # LINKS (restore previous behavior)
            if key2 == "Locus" and tmp:
                tmp = (
                    '<a target="_blank" href="https://ghr.nlm.nih.gov/gene/'
                    + tmp
                    + '#conditions">'
                    + tmp
                    + "</a>"
                )

            if key2 == "dbSNP_ID" and tmp:
                tmp = (
                    '<a target="_blank" href="http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs='
                    + tmp
                    + '">'
                    + tmp
                    + "</a>"
                )

            if key2 == "OMIM_link" and tmp:
                tmp = '<a target="_blank" href="' + tmp + '">' + tmp + "</a>"

            # Optional: keep the old + -> _plus_ safeguard in Mitomap
            if key2 == "Mitomap_Associated_Disease(s)" and tmp:
                tmp = tmp.replace("+", "_plus_")

            # Trim sample suffix
            if key2 == "Sample" and tmp:
                tmp = tmp.replace("-DNA_MIT", "")

            # GT/DP/HF: '|' -> ',<br />'
            if key2 in ("GT", "DP", "HF") and tmp:
                tmp = tmp.replace("|", ",<br />")

            row.append(tmp)

        rows.append(row)

    return rows


def main():
    args = parse_args()

    with open(args.input, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        try:
            header_raw = next(reader)
        except StopIteration:
            sys.exit("Empty input file")

        header = normalize_header(header_raw)
        data_rows = list(reader)

    hash_out = build_hash_out(
        data_rows,
        header,
        hf_cutoff=args.HF,
        maf_cutoff=args.MAF,
        debug=args.debug,
    )

    fmt = args.format

    if fmt == "hash":
        pprint(hash_out)
    elif fmt == "json":
        print(json.dumps(hash_out, sort_keys=True))
    elif fmt == "json4html":
        array = hash2array(hash_out)
        print(json.dumps({"data": array}))
    elif fmt == "tsv":
        writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
        writer.writerow(KEYS2REPORT)
        for key in sorted(hash_out.keys()):
            record = hash_out[key]
            row = [record.get(k, "") for k in KEYS2REPORT]
            writer.writerow(row)
    else:
        sys.exit("Unknown format: {}".format(fmt))


if __name__ == "__main__":
    main()

