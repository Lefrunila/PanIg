#!/usr/bin/env python3
"""
Build species-specific antibody frequency profiles.

This script processes antibody sequences to build amino acid frequency
profiles for each position in the antibody framework regions.

Usage:
    python build_species_profile.py \
        --input data/dog_VH.fasta \
        --species dog \
        --chain-type VH \
        --output profiles/dog_VH.json
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
    max_sequences: int = 0,
) -> dict:
    """
    Build amino acid frequency profile from sequences.

    Args:
        sequences: Dictionary of {name: sequence}
        scheme: Numbering scheme
        max_sequences: Maximum sequences to process (0 = all)

    Returns:
        Dictionary of {position: {amino_acid: frequency}}
    """
    from anarcii import Anarcii

    # Initialize ANARCII once for all sequences
    a = Anarcii(seq_type='antibody', mode='accuracy', cpu=True)

    # Count amino acids at each position
    position_counts = defaultdict(lambda: defaultdict(int))

    total = len(sequences)
    if max_sequences > 0:
        total = min(total, max_sequences)

    # Process sequences in batches
    seq_list = list(sequences.items())
    if max_sequences > 0:
        seq_list = seq_list[:max_sequences]

    batch_size = 50
    for batch_start in range(0, len(seq_list), batch_size):
        batch = seq_list[batch_start:batch_start + batch_size]
        batch_seqs = [seq for _, seq in batch]

        logger.info(f"Processing batch {batch_start//batch_size + 1} ({batch_start + 1}-{min(batch_start + batch_size, len(seq_list))}/{len(seq_list)})...")

        try:
            results = a.number(batch_seqs)
            if results:
                for (name, _), result in zip(batch, results.values()):
                    numbering = result.get("numbering", [])
                    if numbering:
                        for pos_data in numbering:
                            if isinstance(pos_data, (list, tuple)) and len(pos_data) >= 2:
                                pos_info, residue = pos_data[0], pos_data[1]
                                if isinstance(pos_info, (list, tuple)) and len(pos_info) >= 1:
                                    pos = pos_info[0]
                                else:
                                    pos = pos_info
                                if residue != "-" and residue != " ":
                                    position_counts[pos][residue] += 1
        except Exception as e:
            logger.warning(f"Batch processing failed: {e}")
            # Fall back to individual processing
            for name, seq in batch:
                positions = number_sequence_anarcii(seq, scheme)
                for pos, residue in positions:
                    position_counts[pos][residue] += 1

    # Convert counts to frequencies
    profile = {}
    for position, counts in position_counts.items():
        total_count = sum(counts.values())
        if total_count > 0:
            profile[position] = {
                aa: count / total_count
                for aa, count in counts.items()
            }

    return profile


def save_profile(profile: dict, species: str, chain_type: str, output_path: str):
    """Save frequency profile to JSON."""
    data = {
        "species": species,
        "chain_type": chain_type,
        "profile": {str(k): v for k, v in profile.items()},
    }

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)

    logger.info(f"Saved profile to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build species-specific antibody frequency profiles"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with antibody sequences",
    )
    parser.add_argument(
        "--species", "-s",
        required=True,
        help="Species name (e.g., dog, cat)",
    )
    parser.add_argument(
        "--chain-type",
        default="VH",
        choices=["VH", "VL", "VHH"],
        help="Chain type (default: VH)",
    )
    parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "kabat", "chothia", "martin", "aho"],
        help="Numbering scheme (default: imgt)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output JSON file path",
    )
    parser.add_argument(
        "--max-sequences",
        type=int,
        default=0,
        help="Maximum sequences to process (0 = all)",
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
    logger.info(f"Building {args.species} {args.chain_type} profile...")
    profile = build_frequency_profile(sequences, args.scheme, args.max_sequences)

    # Save profile
    save_profile(profile, args.species, args.chain_type, args.output)

    logger.info("Done!")


if __name__ == "__main__":
    main()
