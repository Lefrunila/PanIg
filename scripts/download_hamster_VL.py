#!/usr/bin/env python3
"""
Download hamster immunoglobulin light chain (VL) sequences from NCBI.

Searches for both kappa and lambda chains across hamster species.

Usage:
    python download_hamster_VL.py --output data/hamster --max-sequences 500
"""

import argparse
import logging
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def download_hamster_VL_sequences(output_dir: str, max_sequences: int = 500, email: str = None):
    """
    Download hamster immunoglobulin light chain sequences from NCBI.
    """
    from Bio import Entrez, SeqIO

    if email:
        Entrez.email = email

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Search queries for hamster light chain sequences
    queries = [
        # Kappa light chain searches
        'immunoglobulin kappa light chain[Title] AND Mesocricetus auratus[Organism]',
        'immunoglobulin kappa[Title] AND Mesocricetus auratus[Organism]',
        'kappa light chain[Title] AND Mesocricetus auratus[Organism]',
        'IGK[Gene] AND Mesocricetus auratus[Organism]',
        # Lambda light chain searches
        'immunoglobulin lambda light chain[Title] AND Mesocricetus auratus[Organism]',
        'immunoglobulin lambda[Title] AND Mesocricetus auratus[Organism]',
        'lambda light chain[Title] AND Mesocricetus auratus[Organism]',
        'IGL[Gene] AND Mesocricetus auratus[Organism]',
        # Generic light chain searches
        'immunoglobulin light chain[Title] AND Mesocricetus auratus[Organism]',
        'light chain[Title] AND Mesocricetus auratus[Organism]',
        'immunoglobulin light chain[Title] AND Cricetulus griseus[Organism]',
        'light chain[Title] AND Cricetulus griseus[Organism]',
        'immunoglobulin light chain[Title] AND Cricetulus migratorius[Organism]',
        # Broader hamster searches
        'immunoglobulin light chain[Title] AND hamster[Organism]',
        'light chain[Title] AND hamster[Organism]',
        'antibody light chain[Title] AND hamster[Organism]',
        # Protein-level searches
        'immunoglobulin light chain[Protein] AND Mesocricetus auratus[Organism]',
        'immunoglobulin kappa[Protein] AND Mesocricetus auratus[Organism]',
        'immunoglobulin lambda[Protein] AND Mesocricetus auratus[Organism]',
        'immunoglobulin light chain[Protein] AND Cricetulus griseus[Organism]',
    ]

    all_nuc_ids = set()
    all_prot_ids = set()

    for query in queries:
        is_protein = '[Protein]' in query
        db = 'protein' if is_protein else 'nucleotide'

        try:
            handle = Entrez.esearch(
                db=db,
                term=query,
                retmax=max_sequences,
                idtype='acc'
            )
            record = Entrez.read(handle)
            handle.close()

            ids = record['IdList']
            if ids:
                logger.info(f"Query: {query[:70]}... -> {len(ids)} sequences")
                if is_protein:
                    all_prot_ids.update(ids)
                else:
                    all_nuc_ids.update(ids)

            time.sleep(0.4)

        except Exception as e:
            logger.warning(f"Query failed: {query[:70]}... -> {e}")
            continue

    logger.info(f"\nTotal unique nucleotide IDs found: {len(all_nuc_ids)}")
    logger.info(f"Total unique protein IDs found: {len(all_prot_ids)}")

    # Download nucleotide sequences
    output_file = output_path / "hamster_VL_ncbi.fasta"
    downloaded = 0
    failed = 0

    if all_nuc_ids:
        logger.info(f"Downloading nucleotide sequences to {output_file}...")

        with open(output_file, "w") as f:
            for i, seq_id in enumerate(list(all_nuc_ids)[:max_sequences], 1):
                if i % 10 == 0:
                    logger.info(f"Downloading nuc {i}/{min(len(all_nuc_ids), max_sequences)}...")

                try:
                    handle = Entrez.efetch(
                        db='nucleotide',
                        id=seq_id,
                        rettype='fasta',
                        retmode='text'
                    )
                    record = SeqIO.read(handle, 'fasta')
                    handle.close()

                    SeqIO.write(record, f, 'fasta')
                    downloaded += 1
                    time.sleep(0.4)

                except Exception as e:
                    failed += 1
                    logger.debug(f"Failed {seq_id}: {e}")
                    continue

    # Download protein sequences
    if all_prot_ids:
        prot_file = output_path / "hamster_VL_protein_ncbi.fasta"
        logger.info(f"Downloading protein sequences to {prot_file}...")
        prot_downloaded = 0

        with open(prot_file, "w") as f:
            for i, seq_id in enumerate(list(all_prot_ids)[:max_sequences], 1):
                if i % 10 == 0:
                    logger.info(f"Downloading protein {i}/{min(len(all_prot_ids), max_sequences)}...")

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
                    prot_downloaded += 1
                    time.sleep(0.4)

                except Exception as e:
                    failed += 1
                    logger.debug(f"Failed {seq_id}: {e}")
                    continue

        logger.info(f"Downloaded {prot_downloaded} protein sequences")

    logger.info(f"\nDownload complete!")
    logger.info(f"Nucleotide sequences downloaded: {downloaded}")
    logger.info(f"Failed: {failed}")


def main():
    parser = argparse.ArgumentParser(
        description="Download hamster immunoglobulin light chain sequences from NCBI"
    )
    parser.add_argument("--output", "-o", default="data/hamster", help="Output directory")
    parser.add_argument("--max-sequences", type=int, default=500, help="Maximum sequences to download")
    parser.add_argument("--email", help="Email for NCBI Entrez (recommended)")

    args = parser.parse_args()

    download_hamster_VL_sequences(
        output_dir=args.output,
        max_sequences=args.max_sequences,
        email=args.email,
    )


if __name__ == "__main__":
    main()
