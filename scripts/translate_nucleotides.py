#!/usr/bin/env python3
"""
Translate nucleotide sequences to protein sequences.

This script translates nucleotide sequences from FASTA files
to protein sequences using BioPython's translate function.

Usage:
    python translate_nucleotides.py --input data/cat/cat_VH_ncbi.fasta --output data/cat/cat_VH_protein.fasta
"""

import argparse
import logging
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def translate_nucleotides(input_file: str, output_file: str, frame: int = 0):
    """
    Translate nucleotide sequences to protein sequences.

    Args:
        input_file: Input FASTA file with nucleotide sequences
        output_file: Output FASTA file with protein sequences
        frame: Reading frame (0, 1, or 2). If -1, try all frames and pick best.
    """
    try:
        from Bio import SeqIO
        from Bio.Seq import Seq
    except ImportError:
        logger.error("BioPython not installed. Install with: pip install biopython")
        return

    # Read input sequences
    sequences = []
    with open(input_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            sequences.append(record)

    logger.info(f"Loaded {len(sequences)} nucleotide sequences")

    # Translate sequences
    translated = []
    for record in sequences:
        # Get the sequence
        seq = record.seq

        if frame == -1:
            # Try all three reading frames and pick the best one
            best_protein = None
            best_frame = 0
            best_score = -1

            for f in range(3):
                frame_seq = seq[f:]
                protein = frame_seq.translate()

                # Score: prefer sequences without stop codons
                stop_count = protein.count('*')
                length = len(protein.rstrip('*'))

                # Perfect score: no stops and long sequence
                if stop_count == 0:
                    score = length
                else:
                    # Penalize stop codons heavily
                    score = length - (stop_count * 1000)

                if score > best_score:
                    best_score = score
                    best_protein = protein.rstrip('*')
                    best_frame = f

            protein_seq = best_protein
            if best_frame > 0:
                logger.debug(f"Using frame {best_frame} for {record.id}")
        else:
            # Use specified frame
            if frame > 0:
                seq = seq[frame:]
            protein_seq = seq.translate().rstrip('*')

        # Check if translation is valid
        if protein_seq is None:
            logger.warning(f"Skipping {record.id}: no valid translation found")
            continue

        if len(protein_seq) < 50:
            logger.warning(f"Skipping {record.id}: translated sequence too short ({len(protein_seq)} aa)")
            continue

        # Check for too many stop codons (indicates wrong frame)
        if '*' in protein_seq:
            logger.warning(f"Skipping {record.id}: contains stop codons")
            continue

        # Create new record
        new_record = record[:]
        new_record.seq = protein_seq
        new_record.description = f"{record.description} [translated]"
        translated.append(new_record)

    logger.info(f"Translated {len(translated)} sequences")

    # Write output
    with open(output_file, 'w') as f:
        SeqIO.write(translated, f, 'fasta')

    logger.info(f"Saved translated sequences to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Translate nucleotide sequences to protein sequences"
    )
    parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with nucleotide sequences",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output FASTA file with protein sequences",
    )
    parser.add_argument(
        "--frame",
        type=int,
        default=-1,
        help="Reading frame (0, 1, or 2, or -1 for auto-detect)",
    )

    args = parser.parse_args()

    translate_nucleotides(args.input, args.output, args.frame)


if __name__ == "__main__":
    main()
