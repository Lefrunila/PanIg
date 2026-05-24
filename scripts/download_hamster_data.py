#!/usr/bin/env python3
"""
Download hamster antibody sequences from NCBI.

Hamster antibody data is limited, so we search broadly across
multiple hamster species and antibody-related terms.

Usage:
    python download_hamster_data.py --output data/hamster --max-sequences 5000
"""

import argparse
import logging
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def download_hamster_sequences(output_dir: str, max_sequences: int = 5000, email: str = None):
    """
    Download hamster antibody sequences from NCBI.

    Args:
        output_dir: Output directory
        max_sequences: Maximum sequences to download
        email: Email for NCBI Entrez
    """
    try:
        from Bio import Entrez, SeqIO
    except ImportError:
        logger.error("BioPython not installed. Install with: pip install biopython")
        return

    if email:
        Entrez.email = email

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Broad search queries for hamster antibody sequences
    queries = [
        # Golden/Syrian hamster (Mesocricetus auratus)
        'immunoglobulin heavy chain[Title] AND Mesocricetus auratus[Organism]',
        'immunoglobulin[Title] AND Mesocricetus auratus[Organism]',
        'antibody[Title] AND Mesocricetus auratus[Organism]',
        'VH[Title] AND Mesocricetus auratus[Organism]',
        'immunoglobulin variable[Title] AND Mesocricetus auratus[Organism]',
        'BCR[Title] AND Mesocricetus auratus[Organism]',
        'immunoglobulin[Gene] AND Mesocricetus auratus[Organism]',
        'IGH[Gene] AND Mesocricetus auratus[Organism]',
        # Chinese hamster (Cricetulus griseus)
        'immunoglobulin heavy chain[Title] AND Cricetulus griseus[Organism]',
        'immunoglobulin[Title] AND Cricetulus griseus[Organism]',
        'antibody[Title] AND Cricetulus griseus[Organism]',
        # Armenian hamster
        'immunoglobulin[Title] AND Cricetulus migratorius[Organism]',
        'antibody[Title] AND Cricetulus migratorius[Organism]',
        # Generic hamster searches
        'immunoglobulin heavy chain[Title] AND hamster[Organism]',
        'immunoglobulin[Title] AND hamster[Organism]',
        'antibody[Title] AND hamster[Organism]',
    ]

    all_ids = set()
    query_results = {}

    for query in queries:
        try:
            handle = Entrez.esearch(
                db='nucleotide',
                term=query,
                retmax=max_sequences,
                idtype='acc'
            )
            record = Entrez.read(handle)
            handle.close()

            ids = record['IdList']
            if ids:
                logger.info(f"Query: {query[:70]}... -> {len(ids)} sequences")
                query_results[query] = len(ids)
                all_ids.update(ids)

            time.sleep(0.5)  # Be nice to NCBI

        except Exception as e:
            logger.warning(f"Query failed: {query[:70]}... -> {e}")
            continue

    logger.info(f"\nTotal unique sequences found: {len(all_ids)}")

    if not all_ids:
        logger.warning("No hamster antibody sequences found")
        return

    # Download sequences
    output_file = output_path / "hamster_VH_ncbi.fasta"
    downloaded = 0
    failed = 0

    logger.info(f"Downloading sequences to {output_file}...")

    with open(output_file, "w") as f:
        for i, seq_id in enumerate(list(all_ids)[:max_sequences], 1):
            if i % 10 == 0:
                logger.info(f"Downloading {i}/{min(len(all_ids), max_sequences)}...")

            try:
                handle = Entrez.efetch(
                    db='nucleotide',
                    id=seq_id,
                    rettype='fasta',
                    retmode='text'
                )
                record = SeqIO.read(handle, 'fasta')
                handle.close()

                # Check if it's actually an antibody sequence
                description = record.description.lower()
                is_antibody = any(term in description for term in [
                    'immunoglobulin', 'antibody', 'vh', 'vl', 'ig ', 'igh', 'igk', 'igl',
                    'variable', 'heavy chain', 'light chain', 'fab', 'scfv', 'nanobody'
                ])

                if is_antibody:
                    SeqIO.write(record, f, 'fasta')
                    downloaded += 1
                else:
                    logger.debug(f"Skipping non-antibody: {record.description[:50]}")

                time.sleep(0.5)

            except Exception as e:
                failed += 1
                logger.debug(f"Failed {seq_id}: {e}")
                continue

    logger.info(f"\nDownload complete!")
    logger.info(f"Downloaded: {downloaded} antibody sequences")
    logger.info(f"Failed: {failed}")
    logger.info(f"Output: {output_file}")

    # Also search for hamster protein sequences
    logger.info("\nSearching for hamster protein sequences...")
    protein_queries = [
        'immunoglobulin heavy chain[Protein] AND Mesocricetus auratus[Organism]',
        'immunoglobulin[Protein] AND Mesocricetus auratus[Organism]',
        'antibody[Protein] AND Mesocricetus auratus[Organism]',
    ]

    protein_ids = set()
    for query in protein_queries:
        try:
            handle = Entrez.esearch(
                db='protein',
                term=query,
                retmax=max_sequences,
                idtype='acc'
            )
            record = Entrez.read(handle)
            handle.close()

            ids = record['IdList']
            if ids:
                logger.info(f"Protein query: {query[:70]}... -> {len(ids)} sequences")
                protein_ids.update(ids)

            time.sleep(0.5)

        except Exception as e:
            logger.warning(f"Protein query failed: {e}")
            continue

    if protein_ids:
        protein_file = output_path / "hamster_protein_ncbi.fasta"
        logger.info(f"Downloading {len(protein_ids)} protein sequences...")

        with open(protein_file, "w") as f:
            for i, seq_id in enumerate(list(protein_ids)[:max_sequences], 1):
                if i % 10 == 0:
                    logger.info(f"Downloading protein {i}/{min(len(protein_ids), max_sequences)}...")

                try:
                    handle = Entrez.efetch(
                        db='protein',
                        id=seq_id,
                        rettype='fasta',
                        retmode='text'
                    )
                    record = SeqIO.read(handle, 'fasta')
                    handle.close()

                    SeqIO.write(record, f, 'fasta')
                    time.sleep(0.5)

                except Exception as e:
                    logger.debug(f"Failed {seq_id}: {e}")
                    continue

        logger.info(f"Protein sequences saved to {protein_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Download hamster antibody sequences from NCBI"
    )
    parser.add_argument(
        "--output", "-o",
        default="data/hamster",
        help="Output directory (default: data/hamster)",
    )
    parser.add_argument(
        "--max-sequences",
        type=int,
        default=5000,
        help="Maximum sequences to download (default: 5000)",
    )
    parser.add_argument(
        "--email",
        help="Email for NCBI Entrez (recommended)",
    )

    args = parser.parse_args()

    download_hamster_sequences(
        output_dir=args.output,
        max_sequences=args.max_sequences,
        email=args.email,
    )


if __name__ == "__main__":
    main()
