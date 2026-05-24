#!/usr/bin/env python3
"""
Download species-specific antibody data from OAS (Observed Antibody Space).

This script downloads antibody sequences from OAS for building
species-specific frequency profiles.

Usage:
    python download_oas_data.py --species dog --output data/dog/
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# OAS data URLs (public datasets)
OAS_URLS = {
    "human": {
        "VH": "https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/?species=Human&chain=Heavy",
        "VL": "https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/?species=Human&chain=Light",
    },
    "mouse": {
        "VH": "https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/?species=Mouse&chain=Heavy",
        "VL": "https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/?species=Mouse&chain=Light",
    },
    # Note: OAS primarily has human and mouse data
    # For other species, we need alternative data sources
}


def download_oas_data(species: str, chain_type: str, output_dir: str):
    """
    Download antibody data from OAS.

    Args:
        species: Species name
        chain_type: Chain type ('VH' or 'VL')
        output_dir: Output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    species_lower = species.lower()

    if species_lower not in OAS_URLS:
        logger.warning(
            f"OAS data not directly available for {species}. "
            f"Available: {', '.join(OAS_URLS.keys())}"
        )
        logger.info(
            "For other species, consider using IMGT or literature data. "
            "See scripts/build_species_profile.py for custom data."
        )
        return

    if chain_type not in OAS_URLS[species_lower]:
        logger.warning(f"Chain type {chain_type} not available for {species}")
        return

    url = OAS_URLS[species_lower][chain_type]
    output_file = output_path / f"{species}_{chain_type}_oas.fasta"

    logger.info(f"Downloading {species} {chain_type} data from OAS...")
    logger.info(f"URL: {url}")
    logger.info(f"Output: {output_file}")

    # Note: OAS requires web browser interaction for bulk downloads
    # For automated downloads, consider using their API or pre-downloaded datasets
    logger.info(
        "Note: OAS requires manual download for bulk data. "
        "Please visit the URL above and download the data manually, "
        "then place it in the output directory."
    )


def download_from_imgt(species: str, chain_type: str, output_dir: str):
    """
    Download germline sequences from IMGT.

    IMGT provides germline V-gene sequences for multiple species.

    Args:
        species: Species name
        chain_type: Chain type ('VH' or 'VL')
        output_dir: Output directory
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # IMGT species codes
    imgt_species = {
        "human": "Homo sapiens",
        "mouse": "Mus musculus",
        "dog": "Canis lupus familiaris",
        "cat": "Felis catus",
        "horse": "Equus caballus",
        "cattle": "Bos taurus",
        "pig": "Sus scrofa",
        "sheep": "Ovis aries",
        "goat": "Capra hircus",
        "rabbit": "Oryctolagus cuniculus",
    }

    if species.lower() not in imgt_species:
        logger.warning(f"Species {species} not in IMGT database")
        return

    logger.info(
        f"For IMGT germline sequences for {species}, visit:\n"
        f"http://www.imgt.org/3Dstructure-DB/cgi/DomainGapAlign.cgi\n"
        f"Select species: {imgt_species[species.lower()]}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Download antibody sequence data"
    )
    parser.add_argument(
        "--species", "-s",
        required=True,
        help="Species name (e.g., human, mouse, dog)",
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
        "--source",
        default="oas",
        choices=["oas", "imgt"],
        help="Data source (default: oas)",
    )

    args = parser.parse_args()

    if args.source == "oas":
        download_oas_data(args.species, args.chain_type, args.output)
    elif args.source == "imgt":
        download_from_imgt(args.species, args.chain_type, args.output)


if __name__ == "__main__":
    main()
