#!/usr/bin/env python3
"""
Download veterinary antibody sequences from NCBI.

This script uses Entrez to download antibody heavy chain sequences
for veterinary species (dog, cat, horse, cattle, pig).

Usage:
    python download_veterinary_antibodies.py --species dog --output data/dog/
    python download_veterinary_antibodies.py --all --output data/
"""

import argparse
import logging
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def download_ncbi_sequences(
    species: str,
    chain_type: str,
    output_dir: str,
    max_sequences: int = 1000,
    email: str = None,
):
    """
    Download antibody sequences from NCBI using Entrez.

    Args:
        species: Species name
        chain_type: Chain type ('VH' or 'VL')
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

    # Species-specific search terms
    species_terms = {
        "dog": "Canis lupus familiaris",
        "cat": "Felis catus",
        "horse": "Equus caballus",
        "cattle": "Bos taurus",
        "pig": "Sus scrofa",
        "sheep": "Ovis aries",
        "goat": "Capra hircus",
        "rabbit": "Oryctolagus cuniculus",
        "llama": "Lama glama",
        "alpaca": "Vicugna pacos",
        "camel": "Camelus dromedarius",
    }

    chain_terms = {
        "VH": "immunoglobulin heavy chain",
        "VL": "immunoglobulin light chain",
        "VHH": "single domain antibody",
    }

    if species.lower() not in species_terms:
        logger.error(f"Unknown species: {species}")
        logger.info(f"Available: {', '.join(species_terms.keys())}")
        return

    species_latin = species_terms[species.lower()]
    chain_query = chain_terms.get(chain_type, "immunoglobulin heavy chain")

    # Build search query
    query = f'{chain_query}[Title] AND "{species_latin}"[Organism]'

    logger.info(f"Searching NCBI for {species} {chain_type} sequences...")
    logger.info(f"Query: {query}")

    try:
        # Search for sequences
        handle = Entrez.esearch(
            db="nucleotide",
            term=query,
            retmax=max_sequences,
            idtype="acc",
        )
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]
        logger.info(f"Found {len(id_list)} sequences")

        if not id_list:
            logger.warning(f"No sequences found for {species} {chain_type}")
            return

        # Download sequences
        output_file = output_path / f"{species}_{chain_type}_ncbi.fasta"

        logger.info(f"Downloading sequences to {output_file}...")

        with open(output_file, "w") as f:
            for i, seq_id in enumerate(id_list, 1):
                if i % 10 == 0:
                    logger.info(f"Downloading {i}/{len(id_list)}...")

                try:
                    handle = Entrez.efetch(
                        db="nucleotide",
                        id=seq_id,
                        rettype="fasta",
                        retmode="text",
                    )
                    record = SeqIO.read(handle, "fasta")
                    handle.close()

                    # Write to FASTA
                    SeqIO.write(record, f, "fasta")

                    time.sleep(0.5)  # Be nice to NCBI

                except Exception as e:
                    logger.warning(f"Failed to download {seq_id}: {e}")
                    continue

        logger.info(f"Downloaded {len(id_list)} sequences to {output_file}")

    except Exception as e:
        logger.error(f"Failed to download from NCBI: {e}")
        return


def download_all_species(output_dir: str, max_sequences: int = 1000, email: str = None):
    """
    Download antibody sequences for all veterinary species.

    Args:
        output_dir: Output directory
        max_sequences: Maximum sequences per species
        email: Email for NCBI Entrez
    """
    # Priority species for veterinary use
    priority_species = [
        "dog",
        "cat",
        "horse",
        "cattle",
        "pig",
        "sheep",
        "goat",
        "rabbit",
    ]

    for species in priority_species:
        species_dir = Path(output_dir) / species
        species_dir.mkdir(parents=True, exist_ok=True)

        try:
            download_ncbi_sequences(
                species=species,
                chain_type="VH",
                output_dir=str(species_dir),
                max_sequences=max_sequences,
                email=email,
            )
        except Exception as e:
            logger.warning(f"Failed to download {species}: {e}")
            continue

        time.sleep(2)  # Be nice to NCBI


def main():
    parser = argparse.ArgumentParser(
        description="Download veterinary antibody sequences from NCBI"
    )
    parser.add_argument(
        "--species", "-s",
        help="Species name (e.g., dog, cat, horse)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download all veterinary species",
    )
    parser.add_argument(
        "--chain-type",
        default="VH",
        choices=["VH", "VL", "VHH"],
        help="Chain type (default: VH)",
    )
    parser.add_argument(
        "--output", "-o",
        default="data",
        help="Output directory (default: data/)",
    )
    parser.add_argument(
        "--max-sequences",
        type=int,
        default=1000,
        help="Maximum sequences to download (default: 1000)",
    )
    parser.add_argument(
        "--email",
        help="Email for NCBI Entrez (recommended)",
    )

    args = parser.parse_args()

    if args.all:
        download_all_species(args.output, args.max_sequences, args.email)
    elif args.species:
        download_ncbi_sequences(
            species=args.species,
            chain_type=args.chain_type,
            output_dir=args.output,
            max_sequences=args.max_sequences,
            email=args.email,
        )
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
