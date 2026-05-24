#!/usr/bin/env python3
"""
Basic usage example for PanIg.

Demonstrates how to xenotypize a nanobody sequence to a target species.
"""

from panig import Xenotypizer, SpeciesProfile


def main():
    # Example nanobody sequence (Nb21 from Llamanade test data)
    nanobody_seq = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS"
        "AISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
        "DFYDGSGYWGQGTQVTVSS"
    )

    print("=" * 60)
    print("PanIg: Basic Usage Example")
    print("=" * 60)

    # Initialize xenotypizer
    xenotypizer = Xenotypizer(scheme="imgt", threshold=0.1)

    # Example: Xenotypize to dog species
    # Note: This requires a dog species profile to be available
    # Run 'panig download --species dog' first
    try:
        result = xenotypizer.xenotypize(
            sequence=nanobody_seq,
            target_species="dog",
            name="Nb21",
            chain_type="nanobody",
        )

        print(f"\nOriginal sequence: {nanobody_seq[:50]}...")
        print(f"Xenotypized sequence: {result.xenotypized_sequence[:50]}...")
        print(f"\nTarget species: {result.target_species}")
        print(f"Substitutions: {result.total_substitutions}")
        print(f"Excluded positions: {result.total_excluded}")

        # Save report
        result.to_csv("nb21_dog_report.csv")
        print(f"\nReport saved to: nb21_dog_report.csv")

        # Save animalized sequence
        with open("nb21_dog.fasta", "w") as f:
            f.write(f">{result.animalized_name}\n{result.animalized_sequence}\n")
        print(f"Animalized sequence saved to: nb21_dog.fasta")

    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nPlease download the dog species profile first:")
        print("  panig download --species dog")


if __name__ == "__main__":
    main()
