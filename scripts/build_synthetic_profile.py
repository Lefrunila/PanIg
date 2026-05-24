#!/usr/bin/env python3
"""
Build frequency profiles from synthetic repertoire FASTA files using ANARCII numbering.

Uses ANARCII batch API with correct output format parsing.

Usage:
    python build_synthetic_profile.py --input data/synthetic/cat_VL_synthetic.fasta --output profiles/cat_VL_synthetic.json --species cat
"""

import argparse
import json
import logging
import sys
import time
from collections import defaultdict
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def parse_fasta(fasta_path: str) -> list:
    """Parse a FASTA file into [(name, sequence)] list."""
    sequences = []
    current_name = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    sequences.append((current_name, "".join(current_seq)))
                current_name = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

    if current_name is not None:
        sequences.append((current_name, "".join(current_seq)))

    return sequences


def number_sequences_batch(sequences: list, batch_size: int = 50) -> list:
    """
    Number sequences using ANARCII batch API.
    
    ANARCII.number() returns a dict: {name: {"numbering": [[[pos, ins], aa], ...], ...}}
    
    Returns:
        List of dicts {imgt_position: amino_acid} for each successfully numbered sequence.
    """
    from anarcii import Anarcii
    
    start_time = time.time()
    
    # Create single ANARCII instance
    anarcii = Anarcii(seq_type='antibody', mode='accuracy', cpu=True)
    
    numbered_sequences = []
    total = len(sequences)
    failed = 0
    
    # Process in batches
    for batch_start in range(0, total, batch_size):
        batch_end = min(batch_start + batch_size, total)
        batch = sequences[batch_start:batch_end]
        
        # Prepare input for ANARCII
        seq_pairs = [(name, seq) for name, seq in batch]
        
        # Run ANARCII batch numbering
        try:
            result = anarcii.number(seq_pairs)
            
            # Parse results - result is a dict: {name: {"numbering": [...]}}
            for name, data in result.items():
                try:
                    numbering = data.get('numbering', [])
                    imgt_positions = {}
                    
                    for entry in numbering:
                        if len(entry) >= 2:
                            pos_info = entry[0]
                            residue = entry[1]
                            
                            if len(pos_info) >= 1:
                                imgt_pos = pos_info[0]
                                # Only include valid amino acids
                                if residue and residue not in ('-', '.', 'X', ' '):
                                    imgt_positions[imgt_pos] = residue
                    
                    if imgt_positions:
                        numbered_sequences.append(imgt_positions)
                    else:
                        failed += 1
                except Exception as e:
                    failed += 1
                    if failed <= 5:
                        logger.warning(f"Failed to parse result for {name}: {e}")
        except Exception as e:
            logger.error(f"ANARCII batch failed: {e}")
            failed += len(batch)
        
        # Progress update
        elapsed = time.time() - start_time
        processed = batch_end
        rate = processed / elapsed if elapsed > 0 else 0
        eta = (total - processed) / rate if rate > 0 else 0
        logger.info(f"Progress: {processed}/{total} ({failed} failed) [{rate:.1f} seq/s, ETA: {eta:.0f}s]")
    
    elapsed = time.time() - start_time
    logger.info(f"Successfully numbered {len(numbered_sequences)}/{total} sequences in {elapsed:.1f}s ({failed} failed)")
    
    return numbered_sequences


def build_frequency_profile(numbered_sequences: list) -> dict:
    """
    Build an amino acid frequency profile from IMGT-numbered sequences.
    
    Args:
        numbered_sequences: List of dicts {imgt_position: amino_acid}
    
    Returns:
        Dict of {position: {amino_acid: frequency}}
    """
    position_counts = defaultdict(lambda: defaultdict(int))
    
    for numbered in numbered_sequences:
        for pos, aa in numbered.items():
            position_counts[pos][aa] += 1
    
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


def main():
    parser = argparse.ArgumentParser(
        description="Build frequency profiles from synthetic repertoire FASTA files using ANARCII"
    )
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input FASTA file with synthetic sequences',
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output JSON profile file',
    )
    parser.add_argument(
        '--species', '-s',
        required=True,
        help='Species name (e.g., cat, goat)',
    )
    parser.add_argument(
        '--chain-type',
        default='VL',
        choices=['VL', 'VH'],
        help='Chain type (default: VL)',
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=50,
        help='Batch size for ANARCII numbering (default: 50)',
    )
    
    args = parser.parse_args()
    
    # Create output directory
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    
    # Load sequences
    sequences = parse_fasta(args.input)
    logger.info(f"Loaded {len(sequences)} sequences from {args.input}")
    
    if not sequences:
        logger.error("No sequences found!")
        sys.exit(1)
    
    # Number sequences with ANARCII
    numbered_sequences = number_sequences_batch(sequences, batch_size=args.batch_size)
    
    if not numbered_sequences:
        logger.error("No sequences were successfully numbered!")
        sys.exit(1)
    
    # Build frequency profile
    profile = build_frequency_profile(numbered_sequences)
    logger.info(f"Built profile with {len(profile)} positions")
    
    # Save profile
    data = {
        "species": args.species,
        "chain_type": args.chain_type,
        "profile": {str(k): v for k, v in profile.items()},
        "source": "synthetic_repertoire",
        "scheme": "IMGT",
        "num_synthetic_sequences": len(numbered_sequences),
    }
    
    with open(args.output, 'w') as f:
        json.dump(data, f, indent=2)
    
    logger.info(f"Saved profile to {args.output}")
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Synthetic Profile Summary")
    print(f"{'='*60}")
    print(f"Species: {args.species}")
    print(f"Chain type: {args.chain_type}")
    print(f"Sequences numbered: {len(numbered_sequences)}")
    print(f"Profile positions: {len(profile)}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
