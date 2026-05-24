#!/usr/bin/env python3
"""
Humanization example for PanIg.

Demonstrates how to humanize a nanobody sequence.
"""

from panig import Xenotypizer, SpeciesProfile


def main():
    # Example nanobody sequence (from camelid)
    nanobody_seq = (
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVS"
        "AISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK"
        "DFYDGSGYWGQGTQVTVSS"
    )

    print("=" * 60)
    print("PanIg: Humanization Example")
    print("=" * 60)

    # Initialize xenotypizer
    xenotypizer = Xenotypizer(scheme="imgt", threshold=0.1)

    # Humanize the nanobody
    # Note: This requires the human species profile to be available
    try:
        result = xenotypizer.humanize(
            sequence=nanobody_seq,
            name="camelid_nanobody",
            chain_type="nanobody",
        )

        print(f"\nOriginal sequence: {nanobody_seq[:50]}...")
        print(f"Humanized sequence: {result.modified_sequence[:50]}...")
        print(f"\nSubstitutions: {result.total_substitutions}")
        print(f"Excluded positions: {result.total_excluded}")

        # Save report
        result.to_csv("humanized_report.csv")
        print(f"\nReport saved to: humanized_report.csv")

        # Save humanized sequence
        with open("humanized.fasta", "w") as f:
            f.write(f">{result.modified_name}\n{result.modified_sequence}\n")
        print(f"Humanized sequence saved to: humanized.fasta")

    except FileNotFoundError as e:
        print(f"\nError: {e}")
        print("\nThe human species profile is included in the project.")
        print("If missing, check profiles/human_VH.json")


if __name__ == "__main__":
    main()
