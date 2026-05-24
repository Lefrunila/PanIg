#!/usr/bin/env python3
"""
Batch xenotypization example for PanIg.

Demonstrates how to xenotypize multiple sequences to different species.
"""

import logging
from pathlib import Path

from panig import Xenotypizer, Numberer, SpeciesProfile

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


def main():
    print("=" * 60)
    print("PanIg: Batch Xenotypization Example")
    print("=" * 60)

    # Example sequences
    sequences = {
        "Nb21": (
            "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS"
            "AISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
            "DFYDGSGYWGQGTQVTVSS"
        ),
        "Nb20": (
            "QVQLVESGGGLVQAGGSLRLSCAASGRTFSSYAMGWFRQAPGKEREFVA"
            "AITWSGGNTYYADSVKGRFTISRDNAKNTVYLQMNSLKPEDTAVYYCAA"
            "DRGYYGSGYWGQGTQVTVSS"
        ),
    }

    # Target species
    target_species = ["dog", "cat", "horse"]

    # Initialize xenotypizer
    xenotypizer = Xenotypizer(scheme="imgt", threshold=0.1)

    # Create output directory
    output_dir = Path("batch_results")
    output_dir.mkdir(exist_ok=True)

    # Process each species
    for species in target_species:
        print(f"\n{'=' * 40}")
        print(f"Target species: {species}")
        print(f"{'=' * 40}")

        try:
            # Load species profile
            # Note: This requires species profiles to be available
            # Run 'panig download --species <species>' for each species

            # Xenotypize batch
            results = xenotypizer.xenotypize_batch(
                sequences=sequences,
                target_species=species,
                chain_type="nanobody",
            )

            # Save results
            species_dir = output_dir / species
            species_dir.mkdir(exist_ok=True)

            for result in results:
                # Save animalized sequence
                fasta_path = species_dir / f"{result.original.name}_animalized.fasta"
                with open(fasta_path, "w") as f:
                    f.write(f">{result.animalized_name}\n{result.animalized_sequence}\n")

                # Save report
                report_path = species_dir / f"{result.original.name}_report.csv"
                result.to_csv(str(report_path))

                print(f"\n  {result.original.name}:")
                print(f"    Substitutions: {result.total_substitutions}")
                print(f"    Excluded: {result.total_excluded}")

        except FileNotFoundError as e:
            logger.warning(f"Could not process {species}: {e}")
            logger.info(f"Run 'panig download --species {species}' to get the profile")
            continue

    print(f"\n{'=' * 60}")
    print(f"Results saved to: {output_dir}/")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
