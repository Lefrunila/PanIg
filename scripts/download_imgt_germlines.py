#!/usr/bin/env python3
"""
Download germline sequences from IMGT for veterinary species.

This script downloads V-gene germline sequences from IMGT for building
species-specific antibody frequency profiles and BLAST databases.

Usage:
    python download_imgt_germlines.py --species dog --output data/dog/
    python download_imgt_germlines.py --all --output data/
"""

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# IMGT species codes (Latin name -> IMGT species name)
IMGT_SPECIES = {
    # Veterinary species
    "dog": "Canis lupus familiaris",
    "cat": "Felis catus",
    "horse": "Equus caballus",
    "cattle": "Bos taurus",
    "pig": "Sus scrofa",
    "sheep": "Ovis aries",
    "goat": "Capra hircus",
    "rabbit": "Oryctolagus cuniculus",
    "ferret": "Mustela putorius furo",
    "chicken": "Gallus gallus",
    # Camelids (for nanobodies)
    "llama": "Lama glama",
    "alpaca": "Vicugna pacos",
    "camel": "Camelus dromedarius",
    # Common lab species
    "human": "Homo sapiens",
    "mouse": "Mus musculus",
    "rat": "Rattus norvegicus",
}

# IMGT germline gene types
GERMLINE_TYPES = {
    "IGHV": "Heavy chain V-genes",
    "IGHD": "Heavy chain D-genes",
    "IGHJ": "Heavy chain J-genes",
    "IGLV": "Lambda light chain V-genes",
    "IGLJ": "Lambda light chain J-genes",
    "IGKV": "Kappa light chain V-genes",
    "IGKJ": "Kappa light chain J-genes",
}


def download_imgt_germlines(
    species: str,
    gene_type: str,
    output_dir: str,
):
    """
    Download germline sequences from IMGT.

    Args:
        species: Species name (e.g., 'dog', 'cat')
        gene_type: Gene type (e.g., 'IGHV', 'IGHJ')
        output_dir: Output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if species.lower() not in IMGT_SPECIES:
        logger.error(f"Unknown species: {species}")
        logger.info(f"Available: {', '.join(IMGT_SPECIES.keys())}")
        return

    if gene_type not in GERMLINE_TYPES:
        logger.error(f"Unknown gene type: {gene_type}")
        logger.info(f"Available: {', '.join(GERMLINE_TYPES.keys())}")
        return

    species_latin = IMGT_SPECIES[species.lower()]

    # IMGT API URL for germline sequences
    # Note: IMGT requires specific URL construction
    base_url = "http://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi"

    output_file = output_path / f"{species}_{gene_type}_germlines.fasta"

    logger.info(f"Downloading {species} {gene_type} germlines from IMGT...")
    logger.info(f"Species: {species_latin}")
    logger.info(f"Gene type: {gene_type}")

    # For now, we'll use a simplified approach
    # In production, this would use the IMGT API or scrape the website
    logger.info(
        "Note: IMGT requires web browser interaction for bulk downloads. "
        "Please visit the IMGT website to download germline sequences manually."
    )
    logger.info(f"URL: {base_url}")
    logger.info(f"Select species: {species_latin}")
    logger.info(f"Select gene type: {gene_type}")

    return output_file


def create_imgt_germline_fasta(
    species: str,
    gene_type: str,
    output_file: str,
    sequences: dict = None,
):
    """
    Create a FASTA file with IMGT germline sequences.

    Args:
        species: Species name
        gene_type: Gene type
        output_file: Output FASTA file path
        sequences: Dictionary of {gene_name: sequence}
    """
    if sequences is None:
        sequences = {}

    with open(output_file, "w") as f:
        for gene_name, seq in sequences.items():
            f.write(f">{species}_{gene_type}_{gene_name}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i + 60] + "\n")

    logger.info(f"Created {output_file} with {len(sequences)} sequences")


def download_all_species(output_dir: str, gene_types: list = None):
    """
    Download germlines for all veterinary species.

    Args:
        output_dir: Output directory
        gene_types: List of gene types to download (default: all)
    """
    if gene_types is None:
        gene_types = ["IGHV", "IGHJ"]  # Most important for animalization

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

        for gene_type in gene_types:
            try:
                download_imgt_germlines(species, gene_type, str(species_dir))
            except Exception as e:
                logger.warning(f"Failed to download {species} {gene_type}: {e}")
                continue

        time.sleep(1)  # Be nice to IMGT servers


def main():
    parser = argparse.ArgumentParser(
        description="Download germline sequences from IMGT"
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
        "--gene-type",
        default="IGHV",
        choices=list(GERMLINE_TYPES.keys()),
        help="Gene type to download (default: IGHV)",
    )
    parser.add_argument(
        "--output", "-o",
        default="data",
        help="Output directory (default: data/)",
    )

    args = parser.parse_args()

    if args.all:
        download_all_species(args.output)
    elif args.species:
        download_imgt_germlines(args.species, args.gene_type, args.output)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
