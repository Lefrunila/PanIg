#!/usr/bin/env python3
"""
Build BLAST protein databases with framework-only sequences for T20 scoring.

Extracts FR1+FR2+FR3+FR4 from each sequence using ANARCII numbering,
then builds BLAST databases from the framework-only sequences.

This gives more accurate T20 scores for xenotypization evaluation,
since CDR regions are preserved during xenotypization and shouldn't
affect the species-likeness measurement.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

PANIG_DIR = Path(__file__).parent.parent
MAKEBLASTDB = shutil.which("makeblastdb") or "makeblastdb"

SPECIES = ["dog", "cat", "horse", "cattle", "pig", "sheep", "goat", "rabbit", "hamster", "human"]


def extract_framework(sequence: str, chain_type: str = "heavy") -> str:
    """Extract framework regions from an antibody sequence using ANARCII."""
    from panig.numbering import Numberer

    numberer = Numberer(scheme="imgt")
    try:
        numbered = numberer.number_sequence(sequence, "temp", chain_type)
        return numbered.get_framework_sequence() or ""
    except Exception:
        return ""


def parse_fasta(filepath):
    """Parse a FASTA file, return list of (header, sequence) tuples."""
    records = []
    header = None
    seq_parts = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, ''.join(seq_parts)))
    return records


def write_fasta(records, filepath):
    """Write list of (header, sequence) tuples as FASTA."""
    with open(filepath, 'w') as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')


def build_fr_database(species: str, chain_type: str = "heavy"):
    """Build framework-only BLAST database for a species."""
    suffix = "VH" if chain_type == "heavy" else "VL"
    prot_path = PANIG_DIR / "data" / species / f"{species}_{suffix}_protein.fasta"
    fr_path = PANIG_DIR / "data" / species / f"{species}_{suffix}_protein_FR.fasta"
    db_dir = PANIG_DIR / "blastdb" / f"{species}_{suffix}_blastdb_FR"
    db_prefix = db_dir / f"{species}_{suffix}_FR"

    if not prot_path.exists():
        print(f"  [SKIP] {prot_path} not found")
        return 0

    # Parse protein sequences
    records = parse_fasta(prot_path)
    print(f"  Found {len(records)} sequences in {prot_path.name}")

    # Extract framework regions
    fr_records = []
    skipped = 0
    for header, seq in records:
        fr_seq = extract_framework(seq, chain_type)
        if fr_seq and len(fr_seq) >= 50:
            fr_records.append((header, fr_seq))
        else:
            skipped += 1

    print(f"  Extracted framework from {len(fr_records)} sequences (skipped {skipped})")

    if not fr_records:
        print(f"  [ERROR] No framework sequences extracted")
        return 0

    # Write framework FASTA
    write_fasta(fr_records, fr_path)

    # Build BLAST database
    db_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        MAKEBLASTDB,
        "-in", str(fr_path),
        "-dbtype", "prot",
        "-out", str(db_prefix),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] makeblastdb failed: {result.stderr[:200]}")
        return 0

    total_size = sum(f.stat().st_size for f in db_dir.iterdir() if f.is_file())
    size_kb = total_size / 1024
    print(f"  Built {db_prefix} ({len(fr_records)} sequences, {size_kb:.1f} KB)")
    return len(fr_records)


def main():
    print("=" * 60)
    print("Building framework-only BLAST databases for PanIg T20 scoring")
    print("=" * 60)

    summary = {}
    for species in SPECIES:
        print(f"\n{species} VH:")
        nseq = build_fr_database(species, "heavy")
        summary[f"{species}_VH"] = nseq

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n{'Database':<25} {'Sequences':<12}")
    print("-" * 37)
    for db_name, nseq in summary.items():
        print(f"{db_name:<25} {nseq:<12}")


if __name__ == "__main__":
    main()
