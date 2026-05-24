"""
Command-line interface for PanIg.

Usage:
    panig xenotypize --input sequence.fasta --species dog --output results/
    panig animalize --input sequence.fasta --species dog --output results/  (alias)
    panig score --original orig.fasta --xenotypized xeno.fasta --species dog
    panig download --species dog
"""

import argparse
import logging
import sys
from pathlib import Path

from panig.xenotypizer import Xenotypizer
from panig.numbering import Numberer
from panig.scorer import Scorer
from panig.species_profiles import SpeciesProfile
from panig.structure import StructurePredictor
from panig.vhh_xenotypizer import VHHXenotypizer


def setup_logging(verbose: bool = False):
    """Configure logging."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )


def cmd_xenotypize(args):
    """Handle the 'xenotypize' command."""
    setup_logging(args.verbose)

    # Load input sequence
    sequences = Numberer._parse_fasta(args.input)
    if not sequences:
        print(f"Error: No sequences found in {args.input}")
        sys.exit(1)

    # Determine if we should use VHH xenotypizer
    use_vhh = getattr(args, "vhh", False) or args.chain_type == "nanobody"

    use_synthetic = getattr(args, "use_synthetic", False)

    if use_vhh:
        xenotypizer = VHHXenotypizer(
            threshold=args.threshold,
            scheme=args.scheme,
        )
    else:
        xenotypizer = Xenotypizer(
            threshold=args.threshold,
            scheme=args.scheme,
            use_synthetic=use_synthetic,
        )

    # Load species profile if provided
    species_profile = None
    if args.profile:
        species_profile = SpeciesProfile.load(args.profile)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process each sequence
    for name, sequence in sequences.items():
        print(f"\nProcessing: {name}")

        if use_vhh:
            # VHHXenotypizer expects a NumberedSequence
            numberer = Numberer(scheme=args.scheme)
            numbered = numberer.number_sequence(
                sequence, name, chain_type="nanobody"
            )
            result = xenotypizer.xenotypize(
                numbered_seq=numbered,
                target_species=args.species,
                species_profile=species_profile,
            )
        else:
            result = xenotypizer.xenotypize(
                sequence=sequence,
                target_species=args.species,
                name=name,
                chain_type=args.chain_type,
                species_profile=species_profile,
            )

        # Save xenotypized sequence
        fasta_path = output_dir / f"{name}_xenotypized.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">{result.modified_name}\n{result.modified_sequence}\n")
        print(f"  Xenotypized sequence: {fasta_path}")

        # Save report
        report_path = output_dir / f"{name}_report.csv"
        result.to_csv(str(report_path))
        print(f"  Report: {report_path}")

        # Print summary
        summary = result.summary()
        print(f"  Target species: {summary['target_species']}")
        print(f"  Substitutions: {summary['substitutions']}")
        print(f"  Excluded positions: {summary['excluded']}")
        if use_vhh:
            print(f"  VHH-locked positions: {summary.get('locked_vhh', 0)}")

    print(f"\nDone! Results saved to {output_dir}")


def cmd_humanize(args):
    """Handle the 'humanize' command."""
    setup_logging(args.verbose)

    # Load input sequence
    sequences = Numberer._parse_fasta(args.input)
    if not sequences:
        print(f"Error: No sequences found in {args.input}")
        sys.exit(1)

    # Initialize xenotypizer
    xenotypizer = Xenotypizer(
        threshold=args.threshold,
        scheme=args.scheme,
    )

    # Load human profile if provided
    species_profile = None
    if args.profile:
        species_profile = SpeciesProfile.load(args.profile)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process each sequence
    for name, sequence in sequences.items():
        print(f"\nProcessing: {name}")

        result = xenotypizer.humanize(
            sequence=sequence,
            name=name,
            chain_type=args.chain_type,
            species_profile=species_profile,
        )

        # Save humanized sequence
        fasta_path = output_dir / f"{name}_humanized.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">{result.modified_name}\n{result.modified_sequence}\n")
        print(f"  Humanized sequence: {fasta_path}")

        # Save report
        report_path = output_dir / f"{name}_report.csv"
        result.to_csv(str(report_path))
        print(f"  Report: {report_path}")

        # Print summary
        summary = result.summary()
        print(f"  Substitutions: {summary['substitutions']}")
        print(f"  Excluded positions: {summary['excluded']}")

    print(f"\nDone! Results saved to {output_dir}")



def cmd_predict(args):
    """Handle the 'predict' command."""
    setup_logging(args.verbose)

    # Load input sequence
    sequences = Numberer._parse_fasta(args.input)
    if not sequences:
        print(f"Error: No sequences found in {args.input}")
        sys.exit(1)

    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize structure predictor
    predictor = StructurePredictor(device=args.device)

    # Predict structure for each sequence
    for name, sequence in sequences.items():
        print(f"\nPredicting structure for: {name}")

        try:
            pdb_path = predictor.predict(
                sequence=sequence,
                chain_type=args.chain_type,
                output_path=str(output_dir),
                name=name,
            )
            print(f"  Structure saved: {pdb_path}")

        except ImportError:
            print("  Error: ImmuneBuilder is not installed")
            print("  Install with: pip install immunebuilder")
            sys.exit(1)
        except Exception as e:
            print(f"  Error predicting structure: {e}")
            sys.exit(1)

    print(f"\nDone! Structures saved to {output_dir}")


def cmd_score(args):
    """Handle the 'score' command."""
    setup_logging(args.verbose)

    # Load sequences
    original_seqs = Numberer._parse_fasta(args.original)
    xenotypized_seqs = Numberer._parse_fasta(args.xenotypized)

    if not original_seqs:
        print(f"Error: No sequences found in {args.original}")
        sys.exit(1)
    if not xenotypized_seqs:
        print(f"Error: No sequences found in {args.xenotypized}")
        sys.exit(1)

    # Initialize scorer
    scorer = Scorer(blastdb_path=args.blastdb)

    # Load species profile if species is provided
    species_profile = None
    if args.species:
        from panig.species_profiles import SpeciesProfile
        chain_type = args.chain_type or "heavy"
        chain_map = {"heavy": "VH", "light": "VL", "nanobody": "VHH"}
        suffix = chain_map.get(chain_type, "VH")
        profile_dir = Path(__file__).parent.parent / "profiles"
        profile_path = profile_dir / f"{args.species}_{suffix}.json"
        if profile_path.exists():
            species_profile = SpeciesProfile.load(str(profile_path))

    # Numberer for position coverage
    numberer = Numberer(scheme="imgt") if species_profile else None

    # Score each pair
    for name in original_seqs:
        if name not in xenotypized_seqs:
            print(f"Warning: {name} not found in xenotypized sequences")
            continue

        orig_score, xeno_score = scorer.score_both(
            original=original_seqs[name],
            xenotypized=xenotypized_seqs[name],
        )

        print(f"\n{name}:")
        print(f"  Original T20: {orig_score:.2f}")
        print(f"  Xenotypized T20: {xeno_score:.2f}")
        print(f"  Improvement: {xeno_score - orig_score:.2f}")

        # Compute position coverage if species profile is available
        if species_profile and numberer:
            chain_type = args.chain_type or "heavy"
            xeno_numbered = numberer.number_sequence(
                xenotypized_seqs[name], name, chain_type
            )
            if xeno_numbered:
                coverage, covered, total_fr = scorer.compute_position_coverage(
                    xeno_numbered, species_profile, chain_type
                )
                print(f"  Position Coverage: {coverage:.3f} ({covered}/{total_fr})")


def cmd_download(args):
    """Handle the 'download' command."""
    setup_logging(args.verbose)

    print(f"Downloading databases for {args.species}...")

    # Download BLAST database
    scorer = Scorer()
    try:
        scorer.download_database(
            species=args.species,
            chain_type=args.chain_type,
            force=args.force,
        )
        print(f"BLAST database downloaded successfully")
    except Exception as e:
        print(f"Warning: Could not download BLAST database: {e}")

    # Download species profile
    cache_dir = Path.home() / ".panig" / "cache" / "profiles"
    cache_dir.mkdir(parents=True, exist_ok=True)

    profile_name = f"{args.species}_{args.chain_type}.json"
    remote_path = f"gdrive:PanIg_databases/profiles/{profile_name}"

    try:
        import subprocess
        subprocess.run(
            ["rclone", "copy", remote_path, str(cache_dir)],
            check=True,
            capture_output=True,
        )
        print(f"Species profile downloaded successfully")
    except Exception as e:
        print(f"Warning: Could not download species profile: {e}")

    print("\nDownload complete!")


def cmd_list_species(args):
    """Handle the 'list-species' command."""
    setup_logging(args.verbose)

    # List local profiles
    profile_dir = Path(__file__).parent.parent / "profiles"
    cache_dir = Path.home() / ".panig" / "cache" / "profiles"

    print("Available species profiles:")
    print("-" * 40)

    for source, directory in [
        ("Local", profile_dir),
        ("Cached", cache_dir),
    ]:
        if directory.exists():
            for f in sorted(directory.glob("*.json")):
                species = f.stem.rsplit("_", 1)[0]
                chain_type = f.stem.rsplit("_", 1)[-1]
                print(f"  [{source}] {species} ({chain_type})")


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        prog="panig",
        description="PanIg: Pan-species Immunoglobulin Xenotypization Tool",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging",
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # xenotypize command (primary)
    xenotypize_parser = subparsers.add_parser(
        "xenotypize",
        help="Xenotypize antibody/nanobody sequences to a target species",
    )
    xenotypize_parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with antibody/nanobody sequence(s)",
    )
    xenotypize_parser.add_argument(
        "--species", "-s",
        required=True,
        help="Target species (e.g., dog, cat, horse, cattle, pig)",
    )
    xenotypize_parser.add_argument(
        "--output", "-o",
        default="results",
        help="Output directory (default: results/)",
    )
    xenotypize_parser.add_argument(
        "--chain-type",
        choices=["heavy", "nanobody", "light"],
        default=None,
        help="Chain type (auto-detected if not specified)",
    )
    xenotypize_parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "kabat", "chothia", "martin", "aho"],
        help="Numbering scheme (default: imgt)",
    )
    xenotypize_parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Frequency threshold for native residues (default: 0.1)",
    )
    xenotypize_parser.add_argument(
        "--profile",
        help="Path to species profile JSON file (optional)",
    )
    xenotypize_parser.add_argument(
        "--vhh",
        action="store_true",
        help="Use VHH-specific xenotypizer (preserves nanobody hallmarks)",
    )
    xenotypize_parser.add_argument(
        "--use-synthetic",
        action="store_true",
        help="Prefer synthetic profiles for species with sparse data (e.g. cat/goat VL). Experimental.",
    )
    xenotypize_parser.set_defaults(func=cmd_xenotypize)

    # animalize alias command (backward compatibility)
    animalize_parser = subparsers.add_parser(
        "animalize",
        help="[alias] Xenotypize antibody/nanobody sequences to a target species",
    )
    animalize_parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with antibody/nanobody sequence(s)",
    )
    animalize_parser.add_argument(
        "--species", "-s",
        required=True,
        help="Target species (e.g., dog, cat, horse, cattle, pig)",
    )
    animalize_parser.add_argument(
        "--output", "-o",
        default="results",
        help="Output directory (default: results/)",
    )
    animalize_parser.add_argument(
        "--chain-type",
        choices=["heavy", "nanobody", "light"],
        default=None,
        help="Chain type (auto-detected if not specified)",
    )
    animalize_parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "kabat", "chothia", "martin", "aho"],
        help="Numbering scheme (default: imgt)",
    )
    animalize_parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Frequency threshold for native residues (default: 0.1)",
    )
    animalize_parser.add_argument(
        "--profile",
        help="Path to species profile JSON file (optional)",
    )
    animalize_parser.add_argument(
        "--vhh",
        action="store_true",
        help="Use VHH-specific xenotypizer (preserves nanobody hallmarks)",
    )
    animalize_parser.add_argument(
        "--use-synthetic",
        action="store_true",
        help="Prefer synthetic profiles for species with sparse data (e.g. cat/goat VL). Experimental.",
    )
    animalize_parser.set_defaults(func=cmd_xenotypize)

    # humanize command
    humanize_parser = subparsers.add_parser(
        "humanize",
        help="Humanize antibody/nanobody sequences",
    )
    humanize_parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with antibody/nanobody sequence(s)",
    )
    humanize_parser.add_argument(
        "--output", "-o",
        default="results",
        help="Output directory (default: results/)",
    )
    humanize_parser.add_argument(
        "--chain-type",
        choices=["heavy", "nanobody", "light"],
        default=None,
        help="Chain type (auto-detected if not specified)",
    )
    humanize_parser.add_argument(
        "--scheme",
        default="imgt",
        choices=["imgt", "kabat", "chothia", "martin", "aho"],
        help="Numbering scheme (default: imgt)",
    )
    humanize_parser.add_argument(
        "--threshold",
        type=float,
        default=0.1,
        help="Frequency threshold for native residues (default: 0.1)",
    )
    humanize_parser.add_argument(
        "--profile",
        help="Path to human species profile JSON file (optional)",
    )
    humanize_parser.set_defaults(func=cmd_humanize)

    # score command
    score_parser = subparsers.add_parser(
        "score",
        help="Score original and xenotypized sequences",
    )
    score_parser.add_argument(
        "--original",
        required=True,
        help="Original sequence FASTA file",
    )
    score_parser.add_argument(
        "--xenotypized",
        required=True,
        help="Xenotypized sequence FASTA file",
    )
    score_parser.add_argument(
        "--blastdb",
        help="Path to BLAST database",
    )
    score_parser.add_argument(
        "--species",
        help="Target species for position coverage (e.g., dog, cat, horse)",
    )
    score_parser.add_argument(
        "--chain-type",
        choices=["heavy", "nanobody", "light"],
        default=None,
        help="Chain type for position coverage (default: heavy)",
    )
    score_parser.set_defaults(func=cmd_score)

    # download command
    download_parser = subparsers.add_parser(
        "download",
        help="Download species databases from Google Drive",
    )
    download_parser.add_argument(
        "--species", "-s",
        required=True,
        help="Species to download",
    )
    download_parser.add_argument(
        "--chain-type",
        default="VH",
        choices=["VH", "VHH"],
        help="Chain type (default: VH)",
    )
    download_parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download",
    )
    download_parser.set_defaults(func=cmd_download)

    # list-species command
    list_parser = subparsers.add_parser(
        "list-species",
        help="List available species profiles",
    )
    list_parser.set_defaults(func=cmd_list_species)

    # predict command
    predict_parser = subparsers.add_parser(
        "predict",
        help="Predict 3D structure for antibody/nanobody sequences",
    )
    predict_parser.add_argument(
        "--input", "-i",
        required=True,
        help="Input FASTA file with antibody/nanobody sequence(s)",
    )
    predict_parser.add_argument(
        "--chain-type",
        choices=["nanobody", "antibody"],
        default="nanobody",
        help="Chain type (default: nanobody)",
    )
    predict_parser.add_argument(
        "--output", "-o",
        default="structures",
        help="Output directory for PDB files (default: structures/)",
    )
    predict_parser.add_argument(
        "--device",
        help="PyTorch device (cpu, cuda, cuda:0, etc.)",
    )
    predict_parser.set_defaults(func=cmd_predict)

    # Parse and execute
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
