#!/usr/bin/env python3
"""
Build camelid VHH (nanobody) frequency profile from Llamanade Nb database.

This script processes the 50K nanobody sequences from Llamanade to build
a frequency profile for camelid nanobodies.

Usage:
    python build_camelid_vhh_profile.py \
        --input /path/to/Llamanade/Nb\ database/FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4_Martin_annotated_Nb_50K.fasta \
        --output profiles/camelid_VHH.json
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


def number_sequence_anarcii(sequence: str, scheme: str = "imgt") -> list:
    """
    Number a sequence using ANARCII.

    Returns list of (position, residue) tuples.
    """
    try:
        from anarcii import Anarcii

        # Initialize ANARCII with the correct API
        a = Anarcii(seq_type='antibody', mode='accuracy', cpu=True)

        # Run numbering
        results = a.number([sequence])

        if results and len(results) > 0:
            # Get the first result
            result = list(results.values())[0]
            numbering = result.get("numbering", [])

            if numbering is None:
                return []

            positions = []
            for pos_data in numbering:
                # ANARCII v2 format: ((position, insertion), residue)
                if isinstance(pos_data, (list, tuple)) and len(pos_data) >= 2:
                    pos_info, residue = pos_data[0], pos_data[1]
                    if isinstance(pos_info, (list, tuple)) and len(pos_info) >= 1:
                        pos = pos_info[0]
                    else:
                        pos = pos_info

                    if residue != "-" and residue != " ":
                        positions.append((pos, residue))

            return positions

    except ImportError:
        logger.error("ANARCII not installed. Install with: pip install anarcii")
        return []
    except Exception as e:
        logger.warning(f"ANARCII failed: {e}")
        return []


def build_frequency_profile(
    sequences: dict,
    scheme: str = "imgt",
    max_sequences: int = 10000,
) -> dict:
    """
    Build amino acid frequency profile from sequences.

    Args:
        sequences: Dictionary of {name: sequence}
        scheme: Numbering scheme
        max_sequences: Maximum sequences to process (for speed)

    Returns:
        Dictionary of {position: {amino_acid: frequency}}
    """
    # Count amino acids at each position
    position_counts = defaultdict(lambda: defaultdict(int))

    total = min(len(sequences), max_sequences)
    processed = 0
    failed = 0

    for i, (name, seq) in enumerate(sequences.items(), 1):
        if i > max_sequences:
            break

        if i % 100 == 0:
            logger.info(f"Processing {i}/{total}...")

        positions = number_sequence_anarcii(seq, scheme)
        if positions:
            processed += 1
            for pos, residue in positions:
                position_counts[pos][residue] += 1
        else:
            failed += 1

    logger.info(f"Processed {processed} sequences, failed {failed}")

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


def save_profile(profile: dict, species: str, chain_type: str, output_path: str):
    """Save frequency profile to JSON."""
    data = {
        "species": species,
        "chain_type": chain_type,
        "profile": {str(k): v for k, v in profile.items()},
        "source": "Llamanade Nb database (50K camelid nanobodies)",
        "scheme": "IMGT",
    }

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)

    logger.info(f"Saved profile to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build camelid VHH frequency profile"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file (Llamanade Nb database)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output JSON file path",
    )
    parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "kabat", "chothia", "martin", "aho"],
        help="Numbering scheme (default: imgt)",
    )
    parser.add_argument(
        "--max-sequences",
        type=int,
        default=10000,
        help="Maximum sequences to process (default: 10000)",
    )

    args = parser.parse_args()

    # Load sequences
    logger.info(f"Loading sequences from {args.input}")
    sequences = parse_fasta(args.input)
    logger.info(f"Loaded {len(sequences)} sequences")

    if not sequences:
        logger.error("No sequences found in input file")
        sys.exit(1)

    # Build profile
    logger.info("Building camelid VHH profile...")
    profile = build_frequency_profile(sequences, args.scheme, args.max_sequences)

    # Save profile
    save_profile(profile, "camelid", "VHH", args.output)

    logger.info("Done!")


if __name__ == "__main__":
    main()
