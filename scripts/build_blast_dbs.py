#!/usr/bin/env python3
"""
Build BLAST protein databases for T20 scoring in PanIg.

For each species:
  1. Translate nucleotide VH sequences to protein (best of 3 reading frames)
  2. Create a BLAST protein database with makeblastdb

Also creates a BLAST database from PLAbDab-nano VHH sequences.
"""

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

PANIG_DIR = Path(__file__).parent.parent
PYTHON = sys.executable
MAKEBLASTDB = shutil.which("makeblastdb") or "makeblastdb"

SPECIES = ["dog", "cat", "horse", "cattle", "pig", "sheep", "rabbit"]

# Standard codon table
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def translate_frame(seq, frame):
    """Translate nucleotide sequence in a given reading frame (0, 1, 2)."""
    seq = seq[frame:].upper()
    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        aa = CODON_TABLE.get(codon, 'X')
        protein.append(aa)
    return ''.join(protein)


def best_translation(nuc_seq):
    """
    Try all 3 reading frames, pick the one with fewest stop codons.
    Return (protein_seq, frame) or None if best has stops or < 50 aa.
    """
    best = None
    best_stops = float('inf')
    for frame in range(3):
        protein = translate_frame(nuc_seq, frame)
        stops = protein.count('*')
        if stops < best_stops:
            best_stops = stops
            best = (protein, frame)
    protein, frame = best
    # Remove trailing stop if present
    if protein.endswith('*'):
        protein = protein[:-1]
    if best_stops > 0 or len(protein) < 50:
        return None
    return protein


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
                # Handle malformed headers where sequence is appended
                # Look for a long stretch of nucleotide characters (ATCG) at the end
                content = line[1:]
                # Check if the last token looks like a nucleotide sequence
                # (at least 10 chars of only ATCG)
                import re
                match = re.search(r'\s([ATCG]{10,})\s*$', content)
                if match:
                    # Split: header part + nucleotide part in the header
                    header = content[:match.start()]
                    seq_parts = [match.group(1)]
                else:
                    header = content
                    seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        records.append((header, ''.join(seq_parts)))
    return records


def deduplicate_records(records):
    """Deduplicate records by adding a suffix to duplicate IDs."""
    seen = {}
    deduped = []
    for header, seq in records:
        # Extract the ID (first token)
        seq_id = header.split()[0]
        if seq_id in seen:
            seen[seq_id] += 1
            new_id = f"{seq_id}_{seen[seq_id]}"
            # Replace just the first token
            rest = header[len(seq_id):]
            header = new_id + rest
        else:
            seen[seq_id] = 0
        deduped.append((header, seq))
    return deduped


def write_fasta(records, filepath):
    """Write list of (header, sequence) tuples as FASTA."""
    with open(filepath, 'w') as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')


def translate_ncbi_fasta(species):
    """Translate nucleotide FASTA to protein for a species."""
    nuc_path = PANIG_DIR / "data" / species / f"{species}_VH_ncbi.fasta"
    prot_path = PANIG_DIR / "data" / species / f"{species}_VH_protein.fasta"

    if not nuc_path.exists():
        print(f"  [SKIP] {nuc_path} not found")
        return 0, 0

    records = parse_fasta(nuc_path)
    total = len(records)
    translated = []
    skipped = 0

    for header, seq in records:
        protein = best_translation(seq)
        if protein is None:
            skipped += 1
        else:
            translated.append((header, protein))

    # Deduplicate before writing
    translated = deduplicate_records(translated)
    write_fasta(translated, prot_path)
    print(f"  Translated {len(translated)}/{total} sequences (skipped {skipped}: stops or too short)")
    return len(translated), skipped


def build_blast_db(species, parse_seqids=True):
    """Build BLAST protein database for a species."""
    prot_path = PANIG_DIR / "data" / species / f"{species}_VH_protein.fasta"
    db_dir = PANIG_DIR / "blastdb" / f"{species}_VH_blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / f"{species}_VH"

    if not prot_path.exists():
        print(f"  [SKIP] {prot_path} not found")
        return 0

    # Count sequences
    nseq = sum(1 for line in open(prot_path) if line.startswith('>'))

    cmd = [
        MAKEBLASTDB,
        "-in", str(prot_path),
        "-dbtype", "prot",
        "-out", str(db_prefix),
    ]
    if parse_seqids:
        cmd.append("-parse_seqids")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] makeblastdb failed for {species}: {result.stderr}")
        return 0

    # Get database file sizes
    total_size = sum(f.stat().st_size for f in db_dir.iterdir() if f.is_file())
    size_mb = total_size / (1024 * 1024)
    print(f"  Built {db_prefix} ({nseq} sequences, {size_mb:.1f} MB)")
    return nseq


def build_vhh_blast_db():
    """Build BLAST database from PLAbDab-nano VHH sequences."""
    csv_path = PANIG_DIR / "data" / "vhh" / "vhh_sequences.csv"
    prot_path = PANIG_DIR / "data" / "vhh" / "vhh_VH_protein.fasta"
    db_dir = PANIG_DIR / "blastdb" / "vhh_blastdb"
    db_dir.mkdir(parents=True, exist_ok=True)
    db_prefix = db_dir / "vhh"

    if not csv_path.exists():
        print(f"  [SKIP] {csv_path} not found")
        return 0

    # Parse VHH CSV - extract sequence and ID columns
    records = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            seq_id = row.get('ID', '').strip()
            sequence = row.get('sequence', '').strip()
            if seq_id and sequence:
                header = f"{seq_id}"
                # Add definition if available
                definition = row.get('definition', '').strip()
                if definition:
                    header += f" {definition}"
                records.append((header, sequence))

    # Deduplicate before writing
    records = deduplicate_records(records)
    write_fasta(records, prot_path)
    nseq = len(records)
    print(f"  Extracted {nseq} VHH sequences from CSV")

    # Don't use -parse_seqids for VHH since IDs contain underscores
    # that BLAST interprets as special formats (e.g. 7STG_B)
    cmd = [
        MAKEBLASTDB,
        "-in", str(prot_path),
        "-dbtype", "prot",
        "-out", str(db_prefix),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] makeblastdb failed for VHH: {result.stderr}")
        return 0

    total_size = sum(f.stat().st_size for f in db_dir.iterdir() if f.is_file())
    size_mb = total_size / (1024 * 1024)
    print(f"  Built {db_prefix} ({nseq} sequences, {size_mb:.1f} MB)")
    return nseq


def main():
    print("=" * 60)
    print("Building BLAST databases for PanIg T20 scoring")
    print("=" * 60)

    # Step 1: Translate nucleotide to protein for each species
    print("\n--- Step 1: Translate nucleotide to protein ---")
    translation_summary = {}
    for species in SPECIES:
        prot_path = PANIG_DIR / "data" / species / f"{species}_VH_protein.fasta"
        nuc_path = PANIG_DIR / "data" / species / f"{species}_VH_ncbi.fasta"

        if prot_path.exists() and not nuc_path.exists():
            # Already has protein, no nucleotide to translate
            nseq = sum(1 for line in open(prot_path) if line.startswith('>'))
            print(f"{species}: protein file already exists ({nseq} seqs), skipping translation")
            translation_summary[species] = (nseq, 0)
        elif prot_path.exists() and nuc_path.exists():
            # Cat case: protein exists but might be from a prior run
            # Check if ncbi has more sequences than protein
            nuc_count = sum(1 for line in open(nuc_path) if line.startswith('>'))
            prot_count = sum(1 for line in open(prot_path) if line.startswith('>'))
            if prot_count > 0:
                print(f"{species}: protein file already exists ({prot_count} seqs from {nuc_count} ncbi), keeping as-is")
                translation_summary[species] = (prot_count, 0)
            else:
                print(f"{species}: protein file empty, re-translating...")
                n, s = translate_ncbi_fasta(species)
                translation_summary[species] = (n, s)
        elif nuc_path.exists():
            print(f"{species}: translating from nucleotide...")
            n, s = translate_ncbi_fasta(species)
            translation_summary[species] = (n, s)
        else:
            print(f"{species}: no nucleotide file found, skipping")

    # Step 2: Build BLAST databases for each species
    print("\n--- Step 2: Build BLAST protein databases ---")
    db_summary = {}
    for species in SPECIES:
        print(f"\n{species}:")
        nseq = build_blast_db(species)
        db_summary[species] = nseq

    # Step 3: Build VHH BLAST database
    print("\n--- Step 3: Build VHH BLAST database from PLAbDab-nano ---")
    vhh_nseq = build_vhh_blast_db()

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n{'Species':<12} {'Translated':<12} {'Skipped':<10} {'DB Seqs':<10}")
    print("-" * 44)
    for species in SPECIES:
        trans = translation_summary.get(species, (0, 0))
        db_n = db_summary.get(species, 0)
        print(f"{species:<12} {trans[0]:<12} {trans[1]:<10} {db_n:<10}")
    print(f"{'VHH (PLAbDab)':<12} {'-':<12} {'-':<10} {vhh_nseq:<10}")

    # Report database file sizes
    print(f"\n{'Database':<30} {'Size':<12}")
    print("-" * 42)
    blastdb_dir = PANIG_DIR / "blastdb"
    for db_path in sorted(blastdb_dir.iterdir()):
        if db_path.is_dir() and list(db_path.glob("*.psq")):
            total_size = sum(f.stat().st_size for f in db_path.iterdir() if f.is_file())
            size_mb = total_size / (1024 * 1024)
            print(f"{db_path.name:<30} {size_mb:>8.1f} MB")


if __name__ == "__main__":
    main()
