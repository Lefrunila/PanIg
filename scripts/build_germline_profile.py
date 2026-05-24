#!/usr/bin/env python3
"""
Build species-specific frequency profiles from IMGT germline sequences.

This uses the downloaded IMGT V-gene germline sequences to build
frequency profiles that represent the "ideal" antibody framework
for each species.

Usage:
    python build_germline_profile.py --species dog --output profiles/dog_VH.json
"""

import argparse
import json
import logging
import sys
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_fasta(fasta_path: str) -> dict:
    """Parse a FASTA file into {name: sequence} dictionary."""
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(current_seq)
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        sequences[current_name] = "".join(current_seq)

    return sequences


def is_nucleotide(seq: str) -> bool:
    """Check if a sequence is nucleotide (DNA) rather than protein."""
    # If >80% of characters are A, C, G, T, it's likely DNA
    dna_chars = set('ACGTacgt')
    dna_count = sum(1 for c in seq if c in dna_chars)
    return dna_count / len(seq) > 0.8 if seq else False


def translate_dna_to_protein(seq: str) -> str:
    """Translate a DNA sequence to protein."""
    # Standard genetic code
    codon_table = {
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

    protein = []
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3].upper()
        aa = codon_table.get(codon, 'X')
        if aa == '*':
            break  # Stop codon
        protein.append(aa)

    return ''.join(protein)


def parse_germline_fasta(fasta_path: str) -> dict:
    """
    Parse IMGT germline FASTA file.

    IMGT germline files may contain nucleotide or amino acid sequences.
    If nucleotide, we translate to protein first.
    Then we strip gaps and extract the sequence.

    Returns dict of {gene_name: sequence}
    """
    sequences = {}
    current_name = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    # Remove gaps and convert to uppercase
                    seq = "".join(current_seq).replace(".", "").replace("-", "").upper()
                    if seq:
                        # Check if it's nucleotide and translate if needed
                        if is_nucleotide(seq):
                            seq = translate_dna_to_protein(seq)
                        if seq and len(seq) >= 50:  # Minimum viable protein length
                            sequences[current_name] = seq
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        seq = "".join(current_seq).replace(".", "").replace("-", "").upper()
        if seq:
            if is_nucleotide(seq):
                seq = translate_dna_to_protein(seq)
            if seq and len(seq) >= 50:
                sequences[current_name] = seq

    return sequences


def build_germline_frequency_profile(
    germline_sequences: dict,
    scheme: str = "imgt",
) -> dict:
    """
    Build a frequency profile from germline sequences.

    Since germline sequences are already numbered according to IMGT scheme,
    we can directly count amino acids at each position.

    Args:
        germline_sequences: Dict of {gene_name: sequence}
        scheme: Numbering scheme

    Returns:
        Dict of {position: {amino_acid: frequency}}
    """
    # IMGT positions are 1-128 for VH
    # The germline sequences are already aligned to IMGT numbering
    position_counts = defaultdict(lambda: defaultdict(int))

    for gene_name, seq in germline_sequences.items():
        for i, aa in enumerate(seq):
            position = i + 1  # IMGT positions start at 1
            if aa != "." and aa != "-" and aa != "X":
                position_counts[position][aa] += 1

    # Convert counts to frequencies
    profile = {}
    for position, counts in position_counts.items():
        total = sum(counts.values())
        if total > 0:
            profile[position] = {
                aa: count / total
                for aa, count in counts.items()
            }

    return profile


def merge_profiles(profiles: list) -> dict:
    """
    Merge multiple frequency profiles by averaging.

    Args:
        profiles: List of profile dicts

    Returns:
        Merged profile dict
    """
    merged = defaultdict(lambda: defaultdict(list))

    for profile in profiles:
        for position, aa_freq in profile.items():
            for aa, freq in aa_freq.items():
                merged[position][aa].append(freq)

    result = {}
    for position, aa_freqs in merged.items():
        result[position] = {
            aa: sum(freqs) / len(freqs)
            for aa, freqs in aa_freqs.items()
        }

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Build species-specific frequency profiles from IMGT germlines"
    )
    parser.add_argument(
        "--species", "-s",
        required=True,
        help="Species name (e.g., dog, cat, horse)",
    )
    parser.add_argument(
        "--germline-dir",
        default="data/germlines/imgt",
        help="Directory with IMGT germline FASTA files",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output JSON file path",
    )
    parser.add_argument(
        "--chain-type",
        default="VH",
        choices=["VH", "VL", "VHH"],
        help="Chain type (default: VH)",
    )

    args = parser.parse_args()

    # Map species names to IMGT directory names
    species_map = {
        "dog": "Canis_lupus_familiaris",
        "cat": "Felis_catus",
        "horse": "Equus_caballus",
        "cattle": "Bos_taurus",
        "pig": "Sus_scrofa",
        "sheep": "Ovis_aries",
        "goat": "Capra_hircus",
        "rabbit": "Oryctolagus_cuniculus",
        "alpaca": "Vicugna_pacos",
        "llama": "Vicugna_pacos",
        "dromedary": "Camelus_dromedarius",
    }

    species_dir = species_map.get(args.species.lower())
    if species_dir is None:
        logger.error(f"Unknown species: {args.species}")
        logger.info(f"Available: {', '.join(species_map.keys())}")
        sys.exit(1)

    germline_dir = Path(args.germline_dir) / species_dir

    if not germline_dir.exists():
        logger.error(f"Germline directory not found: {germline_dir}")
        sys.exit(1)

    # Load germline sequences
    if args.chain_type == "VH":
        v_file = germline_dir / "IGHV.fasta"
        if not v_file.exists() or v_file.stat().st_size == 0:
            logger.error(f"Germline file not found or empty: {v_file}")
            sys.exit(1)
        logger.info(f"Loading germline sequences from {v_file}")
        germline_sequences = parse_germline_fasta(str(v_file))
        logger.info(f"Loaded {len(germline_sequences)} germline genes")
    elif args.chain_type == "VL":
        # VL includes both kappa (IGKV) and lambda (IGLV) genes
        germline_sequences = {}
        kappa_file = germline_dir / "IGKV.fasta"
        lambda_file = germline_dir / "IGLV.fasta"
        if kappa_file.exists() and kappa_file.stat().st_size > 0:
            logger.info(f"Loading kappa germline sequences from {kappa_file}")
            kappa_seqs = parse_germline_fasta(str(kappa_file))
            logger.info(f"Loaded {len(kappa_seqs)} kappa germline genes")
            germline_sequences.update(kappa_seqs)
        if lambda_file.exists() and lambda_file.stat().st_size > 0:
            logger.info(f"Loading lambda germline sequences from {lambda_file}")
            lambda_seqs = parse_germline_fasta(str(lambda_file))
            logger.info(f"Loaded {len(lambda_seqs)} lambda germline genes")
            germline_sequences.update(lambda_seqs)
        if not germline_sequences:
            logger.error("No germline sequences found (neither IGKV nor IGLV)")
            sys.exit(1)
        logger.info(f"Total: {len(germline_sequences)} germline genes")
    elif args.chain_type == "VHH":
        # VHH is part of IGHV in camelids
        v_file = germline_dir / "IGHV.fasta"
        if not v_file.exists() or v_file.stat().st_size == 0:
            logger.error(f"Germline file not found or empty: {v_file}")
            sys.exit(1)
        logger.info(f"Loading germline sequences from {v_file}")
        germline_sequences = parse_germline_fasta(str(v_file))
        logger.info(f"Loaded {len(germline_sequences)} germline genes")
    else:
        logger.error(f"Unsupported chain type: {args.chain_type}")
        sys.exit(1)

    if not germline_sequences:
        logger.error("No germline sequences found")
        sys.exit(1)

    # Build frequency profile
    logger.info(f"Building {args.species} {args.chain_type} profile...")
    profile = build_germline_frequency_profile(germline_sequences)

    # Save profile
    data = {
        "species": args.species,
        "chain_type": args.chain_type,
        "profile": {str(k): v for k, v in profile.items()},
        "source": "IMGT germline reference",
        "scheme": "IMGT",
        "num_germline_genes": len(germline_sequences),
    }

    with open(args.output, "w") as f:
        json.dump(data, f, indent=2)

    logger.info(f"Saved profile to {args.output}")
    logger.info(f"Profile covers {len(profile)} positions from {len(germline_sequences)} germline genes")


if __name__ == "__main__":
    main()
