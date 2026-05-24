#!/usr/bin/env python3
"""
Build framework-only BLAST databases using ANARCII numbering with GPU.

Rebuilds FR-only BLAST databases for T20 scoring. Uses ANARCII on GPU
for proper IMGT framework extraction (CDRs excluded).

Usage:
    cd PanIg
    python3 -u scripts/build_fr_blastdbs.py 2>&1 | tee fr_build.log
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

# Force unbuffered output
sys.stdout.reconfigure(line_buffering=True)

# Setup
os.environ["PYTHONUNBUFFERED"] = "1"
PANIG_DIR = Path(__file__).parent.parent
MAKEBLASTDB = shutil.which("makeblastdb") or "makeblastdb"

# IMGT FR positions (1-indexed, inclusive)
FR_POSITIONS = set(range(1, 27)) | set(range(39, 56)) | set(range(66, 105)) | set(range(118, 129))

SPECIES_VH = [
    ("dog", "VH"),
    ("cat", "VH"),
    ("horse", "VH"),
    ("cattle", "VH"),
    ("pig", "VH"),
    ("sheep", "VH"),
    ("goat", "VH"),
    ("rabbit", "VH"),
    ("hamster", "VH"),
    ("human", "VH"),
    ("human", "VL"),  # no human VL BLAST DB needed, skip if not exists
]

# VHH BLAST DB uses alpaca + other camelid sequences
SPECIES_VHH = [
    ("vhh", "VH"),  # vhh_VH_protein.fasta contains camelid VHH sequences
]


def parse_fasta(fp):
    records, header, parts = [], None, []
    with open(fp) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header:
                    records.append((header, "".join(parts)))
                header, parts = line[1:], []
            else:
                parts.append(line)
    if header:
        records.append((header, "".join(parts)))
    return records


def extract_fr_anarcii(anarcii, header, seq):
    """Extract framework from one sequence using ANARCII numbering."""
    try:
        results = anarcii.number([seq])
        entry = results[f"Sequence 1"]
        if entry.get("error"):
            return None
        fr_seq = ""
        for (pos, ins), res in entry["numbering"]:
            if pos in FR_POSITIONS and res != "-":
                fr_seq += res
        return fr_seq if len(fr_seq) >= 70 else None
    except Exception:
        return None


def build_species_fr(anarcii, species, chain_type):
    """Build FR-only BLAST database for one species/chain combo."""
    suffix = chain_type  # VH or VL
    prot_path = PANIG_DIR / "data" / species / f"{species}_{suffix}_protein.fasta"

    if not prot_path.exists():
        print(f"  [SKIP] {prot_path} not found")
        return 0

    records = parse_fasta(prot_path)
    print(f"  {len(records)} sequences in {prot_path.name}")

    fr_path = PANIG_DIR / "data" / species / f"{species}_{suffix}_protein_FR.fasta"
    fr_records = []
    skipped = 0
    start = time.time()

    for i, (header, seq) in enumerate(records):
        fr_seq = extract_fr_anarcii(anarcii, header, seq)
        if fr_seq:
            fr_records.append((header, fr_seq))
        else:
            skipped += 1

        if (i + 1) % 10 == 0:
            elapsed = time.time() - start
            rate = (i + 1) / elapsed
            eta = (len(records) - i - 1) / rate if rate > 0 else 0
            print(
                f"    {i + 1}/{len(records)}  "
                f"ok={len(fr_records)} skip={skipped}  "
                f"{rate:.1f} seq/s  ETA={eta:.0f}s"
            )

    elapsed = time.time() - start
    print(f"  Extracted: {len(fr_records)} FR, skipped: {skipped}, time: {elapsed:.1f}s")

    # Write FR FASTA
    with open(fr_path, "w") as f:
        for header, seq in fr_records:
            f.write(f">{header}\n")
            for j in range(0, len(seq), 60):
                f.write(seq[j : j + 60] + "\n")

    # Build BLAST database
    db_dir = PANIG_DIR / "blastdb" / f"{species}_{suffix}_blastdb_FR"
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / f"{species}_{suffix}_FR"

    cmd = [MAKEBLASTDB, "-in", str(fr_path), "-dbtype", "prot", "-out", str(db_prefix)]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0:
        db_size = sum(f.stat().st_size for f in db_dir.iterdir() if f.is_file()) / 1024
        print(f"  Built {db_prefix.name} ({db_size:.0f} KB)")
    else:
        print(f"  BLAST DB FAILED: {result.stderr[:200]}")

    return len(fr_records)


def main():
    print("=" * 60)
    print("Building FR-only BLAST databases (ANARCII + GPU)")
    print("=" * 60)

    # Initialize ANARCII with GPU
    print("Loading ANARCII on GPU...", flush=True)
    from anarcii import Anarcii

    anarcii = Anarcii(seq_type="antibody", mode="accuracy", cpu=False)
    print("ANARCII loaded on GPU.", flush=True)

    total_start = time.time()
    summary = {}

    for species, chain_type in SPECIES_VH + SPECIES_VHH:
        label = f"{species}_{chain_type}"
        print(f"\n--- {label} ---")
        nseq = build_species_fr(anarcii, species, chain_type)
        summary[label] = nseq

    total_elapsed = time.time() - total_start

    # Summary
    print("\n" + "=" * 60)
    print(f"DONE in {total_elapsed / 60:.1f} min")
    print("=" * 60)
    for label, nseq in summary.items():
        print(f"  {label:<20} {nseq} seqs")


if __name__ == "__main__":
    main()
