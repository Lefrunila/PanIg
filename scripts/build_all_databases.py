#!/usr/bin/env python3
"""
Master script to build all PanIg databases.

This script orchestrates the entire database building process:
1. Download sequences from NCBI
2. Build frequency profiles
3. Create BLAST databases
4. Upload to Google Drive

Usage:
    python build_all_databases.py --email your@email.com --output data/
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Species to process
SPECIES = [
    "dog",
    "cat",
    "horse",
    "cattle",
    "pig",
    "sheep",
    "goat",
    "rabbit",
]


def run_command(cmd: list, description: str) -> bool:
    """Run a command and log the result."""
    logger.info(f"Running: {description}")
    logger.info(f"Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True,
        )
        logger.info(f"Success: {description}")
        return True
    except subprocess.CalledProcessError as e:
        logger.error(f"Failed: {description}")
        logger.error(f"Error: {e.stderr}")
        return False
    except FileNotFoundError as e:
        logger.error(f"Command not found: {e}")
        return False


def download_sequences(species: str, output_dir: str, email: str = None, max_sequences: int = 1000):
    """Download sequences from NCBI."""
    cmd = [
        "python", "scripts/download_veterinary_antibodies.py",
        "--species", species,
        "--output", output_dir,
        "--max-sequences", str(max_sequences),
    ]
    if email:
        cmd.extend(["--email", email])

    return run_command(cmd, f"Download {species} sequences")


def build_profile(species: str, input_file: str, output_file: str):
    """Build frequency profile from sequences."""
    cmd = [
        "python", "scripts/build_species_profile.py",
        "--input", input_file,
        "--species", species,
        "--chain-type", "VH",
        "--output", output_file,
    ]

    return run_command(cmd, f"Build {species} profile")


def create_blastdb(species: str, input_file: str, output_dir: str):
    """Create BLAST database from sequences."""
    db_name = f"{species}_VH_blastdb"
    db_path = Path(output_dir) / db_name

    cmd = [
        "makeblastdb",
        "-in", input_file,
        "-dbtype", "prot",
        "-out", str(db_path),
        "-title", f"{species} VH antibody database",
    ]

    return run_command(cmd, f"Create {species} BLAST database")


def upload_to_gdrive(local_dir: str, remote_dir: str):
    """Upload files to Google Drive using rclone."""
    cmd = [
        "rclone", "copy",
        local_dir,
        f"gdrive:{remote_dir}",
        "--progress",
    ]

    return run_command(cmd, f"Upload to Google Drive: {remote_dir}")


def build_camelid_profile(nb_database: str, output_file: str):
    """Build camelid VHH profile from Llamanade Nb database."""
    cmd = [
        "python", "scripts/build_camelid_vhh_profile.py",
        "--input", nb_database,
        "--output", output_file,
        "--max-sequences", "10000",
    ]

    return run_command(cmd, "Build camelid VHH profile")


def main():
    parser = argparse.ArgumentParser(
        description="Build all PanIg databases"
    )
    parser.add_argument(
        "--email",
        help="Email for NCBI Entrez (recommended)",
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
        help="Maximum sequences per species (default: 1000)",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Specific species to process (default: all)",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Skip sequence download (use existing data)",
    )
    parser.add_argument(
        "--skip-upload",
        action="store_true",
        help="Skip upload to Google Drive",
    )
    parser.add_argument(
        "--nb-database",
        help="Path to Llamanade Nb database for camelid VHH profile",
    )

    args = parser.parse_args()

    # Determine species to process
    species_list = args.species if args.species else SPECIES

    # Create output directories
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    profiles_dir = output_dir / "profiles"
    profiles_dir.mkdir(exist_ok=True)

    blastdb_dir = output_dir / "blastdb"
    blastdb_dir.mkdir(exist_ok=True)

    sequences_dir = output_dir / "sequences"
    sequences_dir.mkdir(exist_ok=True)

    # Track success/failure
    results = {}

    # Process each species
    for species in species_list:
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Processing {species}")
        logger.info(f"{'=' * 60}")

        species_seq_dir = sequences_dir / species
        species_seq_dir.mkdir(exist_ok=True)

        # Step 1: Download sequences
        if not args.skip_download:
            success = download_sequences(
                species=species,
                output_dir=str(species_seq_dir),
                email=args.email,
                max_sequences=args.max_sequences,
            )
            if not success:
                logger.warning(f"Failed to download {species} sequences")
                results[species] = "download_failed"
                continue

        # Step 2: Build profile
        fasta_file = species_seq_dir / f"{species}_VH_ncbi.fasta"
        if not fasta_file.exists():
            logger.warning(f"FASTA file not found: {fasta_file}")
            results[species] = "fasta_missing"
            continue

        profile_file = profiles_dir / f"{species}_VH.json"
        success = build_profile(
            species=species,
            input_file=str(fasta_file),
            output_file=str(profile_file),
        )
        if not success:
            logger.warning(f"Failed to build {species} profile")
            results[species] = "profile_failed"
            continue

        # Step 3: Create BLAST database
        success = create_blastdb(
            species=species,
            input_file=str(fasta_file),
            output_dir=str(blastdb_dir),
        )
        if not success:
            logger.warning(f"Failed to create {species} BLAST database")
            results[species] = "blastdb_failed"
            continue

        results[species] = "success"
        logger.info(f"Successfully processed {species}")

    # Build camelid VHH profile if Nb database provided
    if args.nb_database:
        logger.info(f"\n{'=' * 60}")
        logger.info("Building camelid VHH profile")
        logger.info(f"{'=' * 60}")

        camelid_profile = profiles_dir / "camelid_VHH.json"
        success = build_camelid_profile(
            nb_database=args.nb_database,
            output_file=str(camelid_profile),
        )
        if success:
            results["camelid"] = "success"
        else:
            results["camelid"] = "failed"

    # Upload to Google Drive
    if not args.skip_upload:
        logger.info(f"\n{'=' * 60}")
        logger.info("Uploading to Google Drive")
        logger.info(f"{'=' * 60}")

        # Upload profiles
        upload_to_gdrive(
            str(profiles_dir),
            "PanIg_databases/profiles/",
        )

        # Upload BLAST databases
        upload_to_gdrive(
            str(blastdb_dir),
            "PanIg_databases/blastdb/",
        )

    # Print summary
    logger.info(f"\n{'=' * 60}")
    logger.info("Summary")
    logger.info(f"{'=' * 60}")

    for species, status in results.items():
        logger.info(f"  {species}: {status}")

    success_count = sum(1 for s in results.values() if s == "success")
    total_count = len(results)
    logger.info(f"\nCompleted: {success_count}/{total_count} species")


if __name__ == "__main__":
    main()
