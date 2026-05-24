#!/usr/bin/env python3
"""
Extract immunoglobulin germline V/D/J genes from genome assemblies using BLAST.

Uses mouse germline sequences from IMGT as queries to find orthologous
genes in target species genomes. This is the proper way to extract
germline genes for species without IMGT annotations.

Usage:
    python extract_germlines_blast.py --species hamster --output data/germlines/hamster/
"""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# Genome assembly accessions for target species
GENOME_ASSEMBLIES = {
    "hamster": {
        "name": "BCM_Maur_2.0",
        "accession": "GCF_000349665.1",
        "organism": "Mesocricetus auratus",
        "taxon": 10036,
    },
    "dog": {
        "name": "ROS_Cfam_1.0",
        "accession": "GCF_000002285.5",
        "organism": "Canis lupus familiaris",
        "taxon": 9615,
    },
    "cat": {
        "name": "Felis_catus_9.0",
        "accession": "GCF_000181335.3",
        "organism": "Felis catus",
        "taxon": 9685,
    },
    "horse": {
        "name": "EquCab3.0",
        "accession": "GCF_002863925.1",
        "organism": "Equus caballus",
        "taxon": 9796,
    },
    "cattle": {
        "name": "ARS-UCD1.3",
        "accession": "GCF_002263795.3",
        "organism": "Bos taurus",
        "taxon": 9913,
    },
    "pig": {
        "name": "Sscrofa11.1",
        "accession": "GCF_000003025.6",
        "organism": "Sus scrofa",
        "taxon": 9823,
    },
    "sheep": {
        "name": "Oar_rambouillet_v1.0",
        "accession": "GCF_002742125.1",
        "organism": "Ovis aries",
        "taxon": 9940,
    },
    "goat": {
        "name": "ARS1",
        "accession": "GCF_001704415.1",
        "organism": "Capra hircus",
        "taxon": 9925,
    },
    "rabbit": {
        "name": "OryCun2.0",
        "accession": "GCF_000003625.3",
        "organism": "Oryctolagus cuniculus",
        "taxon": 9986,
    },
    "alpaca": {
        "name": "VicPac3.1",
        "accession": "GCF_000164845.3",
        "organism": "Vicugna pacos",
        "taxon": 30538,
    },
}


def download_genome_ncbi_datasets(species: str, output_dir: str) -> str:
    """
    Download reference genome using NCBI datasets CLI.

    Args:
        species: Species name
        output_dir: Output directory

    Returns:
        Path to downloaded genome FASTA file
    """
    assembly = GENOME_ASSEMBLIES[species]
    accession = assembly["accession"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    genome_file = output_path / f"{species}_genome.fasta"

    if genome_file.exists():
        logger.info(f"Genome already downloaded: {genome_file}")
        return str(genome_file)

    # Try using NCBI datasets CLI
    logger.info(f"Downloading {species} genome ({assembly['name']}, {accession})...")

    try:
        zip_file = output_path / f"{species}_genome.zip"

        cmd = [
            "datasets",
            "download",
            "genome",
            "accession", accession,
            "--include", "genome",
            "--no-progressbar",
            "--filename", str(zip_file),
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if result.returncode != 0:
            logger.warning(f"datasets CLI failed: {result.stderr[:200]}")
            return None

        # Extract the genome
        import zipfile
        if zip_file.exists():
            with zipfile.ZipFile(zip_file, 'r') as z:
                for name in z.namelist():
                    if name.endswith('_genomic.fna') or name.endswith('.fna'):
                        logger.info(f"Extracting {name}...")
                        z.extract(name, output_path)
                        extracted = output_path / name
                        extracted.rename(genome_file)
                        break

            # Clean up zip
            zip_file.unlink(missing_ok=True)

        if genome_file.exists():
            size_mb = genome_file.stat().st_size / (1024 * 1024)
            logger.info(f"Downloaded genome: {genome_file} ({size_mb:.1f} MB)")
            return str(genome_file)

    except (subprocess.TimeoutExpired, FileNotFoundError) as e:
        logger.warning(f"datasets CLI error: {e}")

    return None


def download_genome_from_ncbi_api(species: str, output_dir: str) -> str:
    """
    Download reference genome using NCBI API (fallback).

    Args:
        species: Species name
        output_dir: Output directory

    Returns:
        Path to downloaded genome FASTA file
    """
    assembly = GENOME_ASSEMBLIES[species]
    accession = assembly["accession"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    genome_file = output_path / f"{species}_genome.fasta"

    if genome_file.exists():
        logger.info(f"Genome already downloaded: {genome_file}")
        return str(genome_file)

    # Use NCBI efetch to download genome
    logger.info(f"Downloading genome via NCBI API ({accession})...")

    import requests

    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
    }

    try:
        response = requests.get(url, params=params, timeout=300)
        response.raise_for_status()

        with open(genome_file, 'w') as f:
            f.write(response.text)

        size_mb = genome_file.stat().st_size / (1024 * 1024)
        logger.info(f"Downloaded genome: {genome_file} ({size_mb:.1f} MB)")
        return str(genome_file)

    except requests.RequestException as e:
        logger.warning(f"NCBI API download failed: {e}")

    return None


def make_blastdb(genome_file: str, output_dir: str, species: str) -> str:
    """
    Create BLAST database from genome file.

    Args:
        genome_file: Path to genome FASTA file
        output_dir: Output directory
        species: Species name

    Returns:
        Path to BLAST database
    """
    output_path = Path(output_dir)
    db_path = output_path / f"{species}_genome_db"

    # Check if already exists
    if (db_path.with_suffix('.nsq')).exists():
        logger.info(f"BLAST database already exists: {db_path}")
        return str(db_path)

    logger.info(f"Creating BLAST database...")

    cmd = [
        "makeblastdb",
        "-in", genome_file,
        "-dbtype", "nucl",
        "-out", str(db_path),
        "-title", f"{species} genome",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            logger.info(f"Created BLAST database: {db_path}")
            return str(db_path)
        else:
            logger.error(f"makeblastdb failed: {result.stderr[:200]}")
    except FileNotFoundError:
        logger.error("makeblastdb not found. Install BLAST+: conda install -c bioconda blast")

    return None


def search_germlines_blast(
    query_file: str,
    db_path: str,
    output_file: str,
    gene_type: str,
    evalue: float = 1e-10,
    min_identity: float = 70,
    min_length: int = 50,
) -> int:
    """
    Search for germline sequences using BLAST.

    Args:
        query_file: Path to query FASTA file (mouse germlines)
        db_path: Path to BLAST database
        output_file: Path to output FASTA file
        gene_type: Gene type (IGHV, IGHD, IGHJ, IGKV, IGKJ, IGLV, IGLJ)
        evalue: E-value threshold
        min_identity: Minimum percent identity
        min_length: Minimum sequence length

    Returns:
        Number of sequences found
    """
    blast_output = output_file + ".blast"

    logger.info(f"Searching for {gene_type} germlines...")

    # Run BLAST
    cmd = [
        "blastn",
        "-query", query_file,
        "-db", db_path,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq",
        "-evalue", str(evalue),
        "-max_target_seqs", "500",
        "-word_size", 7,
        "-reward", 1,
        "-penalty", -1,
        "-gapopen", 2,
        "-gapextend", 1,
        "-out", blast_output,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        if result.returncode != 0:
            logger.warning(f"BLAST failed for {gene_type}: {result.stderr[:200]}")
            return 0
    except subprocess.TimeoutExpired:
        logger.warning(f"BLAST timed out for {gene_type}")
        return 0

    # Parse BLAST output and extract sequences
    sequences = {}

    with open(blast_output, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) >= 13:
                query_id = fields[0]
                subject_id = fields[1]
                pident = float(fields[2])
                length = int(fields[3])
                evalue = float(fields[10])
                sequence = fields[12]

                # Filter by identity and length
                if pident >= min_identity and length >= min_length and evalue <= evalue:
                    # Clean up sequence
                    sequence = sequence.replace('-', '').upper()

                    # Use subject ID as key to avoid duplicates
                    if subject_id not in sequences:
                        sequences[subject_id] = sequence

    # Write to FASTA
    with open(output_file, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

    logger.info(f"Found {len(sequences)} {gene_type} sequences")
    return len(sequences)


def extract_germlines_for_species(species: str, output_dir: str):
    """
    Extract all immunoglobulin germline genes for a species.

    Args:
        species: Species name
        output_dir: Output directory
    """
    if species not in GENOME_ASSEMBLIES:
        logger.error(f"Unknown species: {species}")
        return

    assembly = GENOME_ASSEMBLIES[species]
    logger.info(f"\n{'='*60}")
    logger.info(f"Species: {species}")
    logger.info(f"Organism: {assembly['organism']}")
    logger.info(f"Genome: {assembly['name']} ({assembly['accession']})")
    logger.info(f"{'='*60}")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: Download genome
    genome_file = download_genome_ncbi_datasets(species, output_dir)
    if not genome_file:
        genome_file = download_genome_from_ncbi_api(species, output_dir)
    if not genome_file:
        logger.error(f"Failed to download genome for {species}")
        return

    # Step 2: Create BLAST database
    db_path = make_blastdb(genome_file, output_dir, species)
    if not db_path:
        logger.error(f"Failed to create BLAST database for {species}")
        return

    # Step 3: Search for germline genes using mouse as query
    mouse_germline_dir = Path("data/germlines/imgt/Mus_musculus")
    if not mouse_germline_dir.exists():
        logger.error(f"Mouse germline directory not found: {mouse_germline_dir}")
        logger.info("Run IMGT download first: python scripts/download_imgt_germlines.py")
        return

    gene_types = {
        "IGHV": "Heavy chain V genes",
        "IGHD": "Heavy chain D genes",
        "IGHJ": "Heavy chain J genes",
        "IGKV": "Kappa light chain V genes",
        "IGKJ": "Kappa light chain J genes",
        "IGLV": "Lambda light chain V genes",
        "IGLJ": "Lambda light chain J genes",
    }

    total_found = 0

    for gene_type, description in gene_types.items():
        query_file = mouse_germline_dir / f"{gene_type}.fasta"

        if not query_file.exists():
            logger.warning(f"Mouse {gene_type} file not found: {query_file}")
            continue

        output_file = output_path / f"{species}_{gene_type}_genome.fasta"

        # Skip if already extracted
        if output_file.exists():
            count = sum(1 for _ in open(output_file) if _.startswith('>'))
            logger.info(f"{gene_type}: {count} sequences already extracted")
            total_found += count
            continue

        count = search_germlines_blast(
            query_file=str(query_file),
            db_path=db_path,
            output_file=str(output_file),
            gene_type=gene_type,
        )

        total_found += count

    logger.info(f"\nTotal germline sequences extracted: {total_found}")
    logger.info(f"Output directory: {output_dir}")

    # Step 4: Also search for light chain genes from other species if mouse doesn't have them
    # For species like cattle, we might need to use human or other ruminant germlines as queries
    if species in ["cattle", "sheep", "goat"]:
        logger.info(f"\nSearching for ruminant-specific light chain genes...")
        search_ruminant_light_chains(species, db_path, output_path)


def search_ruminant_light_chains(species: str, db_path: str, output_path: Path):
    """
    Search for ruminant-specific light chain genes.

    Ruminants (cattle, sheep, goat) have unique light chain genetics.
    We need to search using ruminant germlines as queries.

    Args:
        species: Species name
        db_path: Path to BLAST database
        output_path: Output directory
    """
    # Use human germlines as backup queries
    human_germline_dir = Path("data/germlines/imgt/Homo_sapiens")

    if not human_germline_dir.exists():
        logger.warning("Human germline directory not found")
        return

    for gene_type in ["IGLV", "IGLJ"]:
        query_file = human_germline_dir / f"{gene_type}.fasta"

        if not query_file.exists():
            continue

        output_file = output_path / f"{species}_{gene_type}_genome_human_query.fasta"

        if output_file.exists():
            continue

        logger.info(f"Searching for {gene_type} using human germlines...")

        count = search_germlines_blast(
            query_file=str(query_file),
            db_path=db_path,
            output_file=str(output_file),
            gene_type=gene_type,
        )

        logger.info(f"Found {count} {gene_type} sequences using human query")


def main():
    parser = argparse.ArgumentParser(
        description="Extract immunoglobulin germline genes from genome assemblies"
    )
    parser.add_argument(
        "--species", "-s",
        help="Species name (e.g., hamster, dog, cat)",
    )
    parser.add_argument(
        "--output", "-o",
        help="Output directory for germline sequences",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Extract germlines for all species",
    )

    args = parser.parse_args()

    if args.all:
        for species in GENOME_ASSEMBLIES:
            output_dir = f"data/germlines/genome/{species}"
            try:
                extract_germlines_for_species(species, output_dir)
            except Exception as e:
                logger.warning(f"Failed to process {species}: {e}")
    elif args.species and args.output:
        extract_germlines_for_species(args.species, args.output)
    else:
        parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
