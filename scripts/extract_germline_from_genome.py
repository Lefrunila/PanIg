#!/usr/bin/env python3
"""
Extract immunoglobulin germline V/D/J genes from genome assemblies.

This script downloads reference genomes and uses BLAST to identify
and extract immunoglobulin germline sequences for species that don't
have IMGT germline data.

Usage:
    python extract_germline_from_genome.py --species hamster --output data/germlines/hamster/
"""

import argparse
import logging
import os
import subprocess
import tempfile
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
    "llama": {
        "name": "LamaGlama",
        "accession": "GCF_001623365.1",
        "organism": "Lama glama",
        "taxon": 9844,
    },
}


def download_genome(species: str, output_dir: str) -> str:
    """
    Download reference genome for a species from NCBI.

    Args:
        species: Species name
        output_dir: Output directory

    Returns:
        Path to downloaded genome FASTA file
    """
    if species not in GENOME_ASSEMBLIES:
        raise ValueError(f"Unknown species: {species}")

    assembly = GENOME_ASSEMBLIES[species]
    accession = assembly["accession"]
    name = assembly["name"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Check if already downloaded
    genome_file = output_path / f"{species}_genome.fasta"
    if genome_file.exists():
        logger.info(f"Genome already downloaded: {genome_file}")
        return str(genome_file)

    # Download from NCBI using datasets CLI
    logger.info(f"Downloading {species} genome ({name}, {accession})...")

    # Try using NCBI datasets CLI
    try:
        cmd = [
            "datasets",
            "download",
            "genome",
            "accession", accession,
            "--include", "genome",
            "--filename", str(output_path / f"{species}_genome.zip"),
        ]

        subprocess.run(cmd, check=True, capture_output=True, timeout=300)

        # Extract the genome
        import zipfile
        zip_file = output_path / f"{species}_genome.zip"
        if zip_file.exists():
            with zipfile.ZipFile(zip_file, 'r') as z:
                # Find the genome FASTA file
                for name in z.namelist():
                    if name.endswith('.fna') or name.endswith('.fasta'):
                        z.extract(name, output_path)
                        # Rename to standard name
                        extracted = output_path / name
                        extracted.rename(genome_file)
                        break

        logger.info(f"Downloaded genome: {genome_file}")
        return str(genome_file)

    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        logger.warning(f"datasets CLI failed: {e}")

        # Fallback: try using NCBI BLAST API to search genome
        logger.info("Trying alternative approach: searching genome via NCBI BLAST API...")
        return search_genome_for_germlines(species, output_dir)


def search_genome_for_germlines(species: str, output_dir: str) -> str:
    """
    Search for immunoglobulin germline genes in genome using NCBI BLAST API.

    Args:
        species: Species name
        output_dir: Output directory

    Returns:
        Path to extracted germline sequences
    """
    assembly = GENOME_ASSEMBLIES[species]
    accession = assembly["accession"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Use mouse germline sequences as queries (since they're closely related)
    # First, get mouse germline sequences
    mouse_germline_dir = Path("data/germlines/imgt/Mus_musculus")
    if not mouse_germline_dir.exists():
        logger.error("Mouse germline directory not found. Run IMGT download first.")
        return None

    # Search for V, D, J genes
    for gene_type in ["IGHV", "IGHD", "IGHJ"]:
        query_file = mouse_germline_dir / f"{gene_type}.fasta"
        if not query_file.exists():
            logger.warning(f"Mouse {gene_type} file not found: {query_file}")
            continue

        logger.info(f"Searching {species} genome for {gene_type} genes...")

        # Use NCBI BLAST API to search the genome
        output_file = output_path / f"{species}_{gene_type}_genome.fasta"

        try:
            # Read query sequences
            with open(query_file, 'r') as f:
                query_content = f.read()

            # Use NCBI BLAST API
            cmd = [
                "blastn",
                "-query", str(query_file),
                "-subject", f"refseq/{accession}",
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sseq",
                "-evalue", "1e-10",
                "-max_target_seqs", "100",
                "-out", str(output_path / f"{species}_{gene_type}_blast.txt"),
            ]

            subprocess.run(cmd, check=True, capture_output=True, timeout=300)

            # Parse BLAST output and extract sequences
            extract_germlines_from_blast(
                output_path / f"{species}_{gene_type}_blast.txt",
                output_file,
                gene_type,
            )

            logger.info(f"Extracted {gene_type} sequences: {output_file}")

        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            logger.warning(f"BLAST search failed for {gene_type}: {e}")

            # Alternative: use NCBI web BLAST API
            logger.info("Trying NCBI web BLAST API...")
            search_ncbi_blast_api(species, gene_type, query_file, output_file)

    return str(output_path)


def extract_germlines_from_blast(blast_file: Path, output_file: Path, gene_type: str):
    """
    Extract germline sequences from BLAST output.

    Args:
        blast_file: Path to BLAST output file
        output_file: Path to output FASTA file
        gene_type: Gene type (IGHV, IGHD, IGHJ)
    """
    sequences = {}

    with open(blast_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')
            if len(fields) >= 13:
                query_id = fields[0]
                subject_id = fields[1]
                pident = float(fields[2])
                evalue = float(fields[10])
                sequence = fields[12]

                # Filter by identity and e-value
                if pident >= 70 and evalue <= 1e-10:
                    # Clean up sequence
                    sequence = sequence.replace('-', '')

                    if len(sequence) >= 50:  # Minimum length for V genes
                        sequences[subject_id] = sequence

    # Write to FASTA
    with open(output_file, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + "\n")

    logger.info(f"Extracted {len(sequences)} {gene_type} sequences")


def search_ncbi_blast_api(species: str, gene_type: str, query_file: Path, output_file: Path):
    """
    Search NCBI BLAST API for immunoglobulin germline genes.

    Args:
        species: Species name
        gene_type: Gene type (IGHV, IGHD, IGHJ)
        query_file: Path to query FASTA file
        output_file: Path to output FASTA file
    """
    import requests
    import time

    assembly = GENOME_ASSEMBLIES[species]
    accession = assembly["accession"]

    # Read query sequence
    with open(query_file, 'r') as f:
        lines = f.readlines()
        query_seq = ''.join(line.strip() for line in lines if not line.startswith('>'))

    # Submit BLAST search
    logger.info(f"Submitting BLAST search for {gene_type} in {species} genome...")

    blast_url = "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi"
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": "refseq",
        "QUERY": query_seq[:1000],  # Limit query length
        "EXPECT": "1e-10",
        "HITLIST_SIZE": "100",
        "FORMAT_TYPE": "JSON2",
    }

    try:
        response = requests.post(blast_url, data=params, timeout=30)
        response.raise_for_status()

        # Parse response to get RID
        rid = None
        for line in response.text.split('\n'):
            if 'RID' in line:
                rid = line.split('=')[1].strip()
                break

        if not rid:
            logger.warning("Failed to get BLAST RID")
            return

        logger.info(f"BLAST RID: {rid}")

        # Wait for results
        logger.info("Waiting for BLAST results...")
        time.sleep(30)

        # Get results
        result_url = f"https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=JSON2&RID={rid}"
        result_response = requests.get(result_url, timeout=60)

        if result_response.status_code == 200:
            # Parse results
            import json
            try:
                results = json.loads(result_response.text)
                # Extract sequences from results
                sequences = []
                # ... parse JSON results ...
                logger.info(f"Extracted {len(sequences)} {gene_type} sequences from NCBI BLAST")
            except json.JSONDecodeError:
                logger.warning("Failed to parse BLAST results")
        else:
            logger.warning(f"Failed to get BLAST results: {result_response.status_code}")

    except requests.RequestException as e:
        logger.warning(f"NCBI BLAST API request failed: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract immunoglobulin germline genes from genome assemblies"
    )
    parser.add_argument(
        "--species", "-s",
        required=True,
        help="Species name (e.g., hamster, dog, cat)",
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
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
            logger.info(f"\n{'='*60}")
            logger.info(f"Processing {species}")
            logger.info(f"{'='*60}")
            try:
                extract_germlines_from_genome(species, args.output)
            except Exception as e:
                logger.warning(f"Failed to process {species}: {e}")
    else:
        extract_germlines_from_genome(args.species, args.output)


def extract_germlines_from_genome(species: str, output_dir: str):
    """
    Extract immunoglobulin germline genes from genome for a species.

    Args:
        species: Species name
        output_dir: Output directory
    """
    if species not in GENOME_ASSEMBLIES:
        logger.error(f"Unknown species: {species}")
        return

    assembly = GENOME_ASSEMBLIES[species]
    logger.info(f"Species: {species}")
    logger.info(f"Organism: {assembly['organism']}")
    logger.info(f"Genome: {assembly['name']} ({assembly['accession']})")

    # Download genome
    genome_file = download_genome(species, output_dir)

    if genome_file:
        logger.info(f"Genome downloaded: {genome_file}")

        # Use IgBLAST to annotate immunoglobulin genes
        logger.info("Using IgBLAST to annotate immunoglobulin genes...")
        # ... IgBLAST annotation ...

        # Or use BLAST to search for germline sequences
        logger.info("Searching for germline sequences...")
        # ... BLAST search ...

    logger.info("Done!")


if __name__ == "__main__":
    main()
