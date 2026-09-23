#!/usr/bin/env python3
"""HOMD reference acquisition and release-specific parsing.

Requires numpy and scikit-learn. Run: python homd_data.py
Downloads four public HOMD files to data/ if absent.
"""
from __future__ import annotations

import hashlib
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import normalize


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)
BASE = "https://www.homd.org/ftp/16S_rRNA_refseq/HOMD_16S_rRNA_RefSeq/"
FILES = {
    "new_fasta": "V16.03/HOMD_16S_rRNA_RefSeq_V16.03_full.fasta",
    "new_tax": "V16.03/HOMD_16S_rRNA_RefSeq_V16.03.qiime.taxonomy",
    "old_fasta": "V15.23/HOMD_16S_rRNA_RefSeq_V15.23.p9.fasta",
    "old_tax": "V15.23/HOMD_16S_rRNA_RefSeq_V15.23.qiime.taxonomy",
}
SEED = 20260923


def acquire():
    for name, relative in FILES.items():
        target = DATA / Path(relative).name
        if not target.exists():
            with urlopen(Request(BASE + relative, headers={"User-Agent": "HOMD-release-benchmark/1.0"}), timeout=60) as response:
                target.write_bytes(response.read())
        print(name, target.stat().st_size, "bytes")


def read_release(fasta: Path, taxonomy: Path, new: bool):
    tax = {}
    for line in taxonomy.read_text().splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            tax[key] = value
    rows = []
    header, seq = None, []

    def add():
        if header is None:
            return
        ident = header.split()[0] if new else header.split(" | ")[0]
        match = re.search(r"HMT-\d+", header)
        if not match or ident not in tax:
            raise ValueError("Unmatched HOMD sequence or taxonomy: " + header)
        hierarchy = dict(p.split("__", 1) for p in tax[ident].split(";") if "__" in p)
        rows.append({"id": ident, "hmt": match.group(), "seq": "".join(seq).upper(),
                     "genus": hierarchy.get("g", ""), "species": hierarchy.get("s", "")})

    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            add()
            header, seq = line[1:], []
        else:
            seq.append(line.strip())
    add()
    return rows



if __name__ == "__main__":
    acquire()
