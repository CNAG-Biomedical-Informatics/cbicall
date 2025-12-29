#!/usr/bin/env python3
"""
mtb2json.py - Parse MToolBox prioritized variants output

Usage:
  mtb2json.py -i mit_prioritized_variants.txt -f json4html > mtDNA.json
  mtb2json.py -i mit_prioritized_variants.txt -f json > mit.raw.json

Notes:
  - By default, HF filtering (if enabled) is applied using the maximum HF found
    in ANY sample present in the HF field. This works for single-sample and
    multi-sample/cohort inputs.
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
        default=0.30,  # existing default behavior
        help="Heteroplasmic Fraction cutoff [0.30]. Set to 0 to disable.",
    )
    p.add_argument(
        "--MAF",
        type=float,
        default=0.01,  # existing default behavior
        help="MAF filter [0.01]; keep variant if 1000G MAF < MAF. Set to 0 to disable.",
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
    """Convert header names like 'Variant Allele' to 'Variant_Allele'."""
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


def max_hf_any_sample(hf_str):
    """
    HF patterns:
      single mode => "0.02,0.08"
      cohort mode => "01F:0.02|01M:N/A|01P:1.0"  (order may vary)
      cohort multiallelic => "01P:0.161,0.387|01M:0.161,0.387|01F:0.161,0.387"
    Consider ANY sample. Turn N/A to 0. Take maximum numeric HF found anywhere.
    """
    if not hf_str:
        return 0.0

    s = hf_str.strip().replace("N/A", "0")

    # Single-sample case (no per-sample blocks)
    if "|" not in s and ":" not in s:
        vals = []
        for x in [p.strip() for p in s.split(",") if p.strip()]:
            try:
                vals.append(float(x))
            except ValueError:
                continue
        return max(vals) if vals else 0.0

    # Cohort case: parse each block after "SAMPLE:"
    max_val = 0.0
    for block in [b.strip() for b in s.split("|") if b.strip()]:
        rhs = block.split(":", 1)[1] if ":" in block else block
        for x in [p.strip() for p in rhs.split(",") if p.strip()]:
            try:
                v = float(x)
            except ValueError:
                continue
            if v > max_val:
                max_val = v
    return max_val


def passes_hf_filter(row, col_idx, hf_cutoff):
    if hf_cutoff <= 0:
        return True
    idx = col_idx.get("HF")
    if idx is None:
        return True
    hf_str = row[idx].strip()
    m = max_hf_any_sample(hf_str)
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

        # Normalize sample (cosmetic ordering)
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

            # LINKS
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

