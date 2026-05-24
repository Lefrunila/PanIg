#!/usr/bin/env python3
"""
Download hamster immunoglobulin light chain (VL) sequences from NCBI v2.
Searches protein database first, falls back to CDS nucleotide sequences.
Saves protein sequences to hamster_VL_protein.fasta.
Saves nucleotide CDS sequences to hamster_VL_ncbi.fasta.
"""

import argparse
import logging
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def download_sequences(output_dir: str, max_sequences: int = 500, email: str = None):
    from Bio import Entrez, SeqIO

    if email:
        Entrez.email = email

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Protein-focused queries for hamster light chains
    queries = [
        ('protein', 'immunoglobulin kappa light chain[Title] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin lambda light chain[Title] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin light chain[Title] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin kappa[Title] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin lambda[Title] AND Mesocricetus auratus[Organism]'),
        ('protein', 'light chain variable[Title] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin kappa chain[Protein Name] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin lambda chain[Protein Name] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin light chain[Protein Name] AND Mesocricetus auratus[Organism]'),
        ('protein', 'IGLV[Gene] AND Mesocricetus auratus[Organism]'),
        ('protein', 'IGKV[Gene] AND Mesocricetus auratus[Organism]'),
        ('protein', 'immunoglobulin light chain[Title] AND Cricetulus griseus[Organism]'),
        ('protein', 'immunoglobulin kappa[Title] AND Cricetulus griseus[Organism]'),
        ('protein', 'immunoglobulin lambda[Title] AND Cricetulus griseus[Organism]'),
        ('protein', 'immunoglobulin light chain[Protein Name] AND Cricetulus griseus[Organism]'),
        ('protein', 'immunoglobulin light chain[Title] AND Cricetulus migratorius[Organism]'),
        ('protein', 'immunoglobulin light chain[Title] AND hamster[Organism]'),
        ('protein', 'immunoglobulin light chain[Protein Name] AND hamster[Organism]'),
        ('protein', 'antibody light chain[Title] AND hamster[Organism]'),
        # CDS nucleotide queries
        ('nucleotide', 'immunoglobulin light chain[Title] AND Mesocricetus auratus[Organism] AND cds[Feature]'),
        ('nucleotide', 'immunoglobulin kappa[Title] AND Mesocricetus auratus[Organism] AND cds[Feature]'),
        ('nucleotide', 'immunoglobulin lambda[Title] AND Mesocricetus auratus[Organism] AND cds[Feature]'),
        ('nucleotide', 'immunoglobulin light chain[Title] AND Cricetulus griseus[Organism] AND cds[Feature]'),
        ('nucleotide', 'immunoglobulin light chain[Title] AND Cricetulus migratorius[Organism] AND cds[Feature]'),
        ('nucleotide', 'immunoglobulin light chain[Title] AND hamster[Organism] AND cds[Feature]'),
        ('nucleotide', 'kappa light chain[Title] AND Mesocricetus auratus[Organism]'),
        ('nucleotide', 'lambda light chain[Title] AND Mesocricetus auratus[Organism]'),
        ('nucleotide', 'light chain[Title] AND Mesocricetus auratus[Organism]'),
        ('nucleotide', 'immunoglobulin light chain[Title] AND hamster[Organism]'),
    ]

    all_ids = {}  # db -> set of ids
    for db, query in queries:
        try:
            handle = Entrez.esearch(db=db, term=query, retmax=max_sequences, idtype='acc')
            record = Entrez.read(handle)
            handle.close()
            ids = record['IdList']
            if ids:
                logger.info(f"{db}: {query[:70]}... -> {len(ids)} sequences")
                if db not in all_ids:
                    all_ids[db] = set()
                all_ids[db].update(ids)
            time.sleep(0.4)
        except Exception as e:
            logger.warning(f"Failed: {query[:70]}... -> {e}")
            continue

    prot_ids = all_ids.get('protein', set())
    nuc_ids = all_ids.get('nucleotide', set())
    logger.info(f"\nTotal unique protein IDs: {len(prot_ids)}")
    logger.info(f"Total unique nucleotide IDs: {len(nuc_ids)}")

    # Download protein sequences
    prot_count = 0
    prot_file = output_path / "hamster_VL_protein.fasta"
    if prot_ids:
        logger.info(f"Downloading protein sequences to {prot_file}...")
        with open(prot_file, "w") as f:
            for i, seq_id in enumerate(list(prot_ids)[:max_sequences], 1):
                if i % 20 == 0:
                    logger.info(f"Downloading protein {i}/{min(len(prot_ids), max_sequences)}...")
                try:
                    handle = Entrez.efetch(db='protein', id=seq_id, rettype='fasta', retmode='text')
                    record = SeqIO.read(handle, 'fasta')
                    handle.close()
                    SeqIO.write(record, f, 'fasta')
                    prot_count += 1
                    time.sleep(0.4)
                except Exception as e:
                    logger.debug(f"Failed {seq_id}: {e}")
                    continue
        logger.info(f"Downloaded {prot_count} protein sequences")

    # Download nucleotide CDS sequences
    nuc_count = 0
    nuc_file = output_path / "hamster_VL_ncbi.fasta"
    if nuc_ids:
        logger.info(f"Downloading nucleotide sequences to {nuc_file}...")
        with open(nuc_file, "w") as f:
            for i, seq_id in enumerate(list(nuc_ids)[:max_sequences], 1):
                if i % 20 == 0:
                    logger.info(f"Downloading nuc {i}/{min(len(nuc_ids), max_sequences)}...")
                try:
                    handle = Entrez.efetch(db='nucleotide', id=seq_id, rettype='fasta', retmode='text')
                    record = SeqIO.read(handle, 'fasta')
                    handle.close()
                    SeqIO.write(record, f, 'fasta')
                    nuc_count += 1
                    time.sleep(0.4)
                except Exception as e:
                    logger.debug(f"Failed {seq_id}: {e}")
                    continue
        logger.info(f"Downloaded {nuc_count} nucleotide sequences")

    logger.info(f"\nSummary:")
    logger.info(f"  Protein sequences: {prot_count}")
    logger.info(f"  Nucleotide sequences: {nuc_count}")


def main():
    parser = argparse.ArgumentParser(description="Download hamster VL sequences from NCBI (v2)")
    parser.add_argument("--output", "-o", default="data/hamster")
    parser.add_argument("--max-sequences", type=int, default=500)
    parser.add_argument("--email")
    args = parser.parse_args()
    download_sequences(args.output, args.max_sequences, args.email)


if __name__ == "__main__":
    main()
