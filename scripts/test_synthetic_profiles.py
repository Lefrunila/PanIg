#!/usr/bin/env python3
"""
Test synthetic profiles by xenotypizing Trastuzumab VL to cat and goat.

Compares T20 scores between germline-only and synthetic profiles.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from panig.xenotypizer import Xenotypizer
from panig.scorer import Scorer
from panig.species_profiles import SpeciesProfile

PANIG_DIR = Path(__file__).parent.parent


def main():
    # Trastuzumab VL sequence
    trastuzumab_VL = 'DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK'

    print("=" * 80)
    print("Testing Synthetic Profiles for Cat and Goat VL")
    print("=" * 80)
    print(f"\nTrastuzumab VL: {trastuzumab_VL[:50]}...")
    print(f"Length: {len(trastuzumab_VL)}")

    # Test configurations
    tests = [
        {
            'species': 'cat',
            'chain_type': 'light',
            'blastdb': str(PANIG_DIR / 'blastdb' / 'cat_VL_blastdb' / 'cat_VL'),
            'germline_profile': str(PANIG_DIR / 'profiles' / 'cat_VL.json'),
            'synthetic_profile': str(PANIG_DIR / 'profiles' / 'cat_VL_synthetic.json'),
        },
        {
            'species': 'goat',
            'chain_type': 'light',
            'blastdb': str(PANIG_DIR / 'blastdb' / 'goat_VL_blastdb' / 'goat_VL'),
            'germline_profile': str(PANIG_DIR / 'profiles' / 'goat_VL.json'),
            'synthetic_profile': str(PANIG_DIR / 'profiles' / 'goat_VL_synthetic.json'),
        },
    ]
    
    for test in tests:
        species = test['species']
        print(f"\n{'='*80}")
        print(f"Testing {species.upper()} VL")
        print(f"{'='*80}")
        
        # Initialize components
        xeno = Xenotypizer(scheme='imgt', threshold=0.1)
        scorer = Scorer(blastdb_path=test['blastdb'])
        
        # 1. Baseline: score original Trastuzumab VL
        orig_score = scorer.score_sequence(trastuzumab_VL)
        print(f"\n1. Baseline (Trastuzumab VL):")
        print(f"   Original T20: {orig_score:.2f}")
        
        # 2. Xenotypize with germline-only profile
        germline_profile = SpeciesProfile.load(test['germline_profile'])
        result_germline = xeno.xenotypize(
            trastuzumab_VL,
            species,
            'Trastuzumab_VL',
            test['chain_type'],
            species_profile=germline_profile,
        )
        germline_score = scorer.score_sequence(result_germline.xenotypized_sequence)
        print(f"\n2. Germline-only profile:")
        print(f"   Xenotypized T20: {germline_score:.2f}")
        print(f"   Improvement: {germline_score - orig_score:.2f}")
        print(f"   Substitutions: {result_germline.total_substitutions}")
        
        # 3. Xenotypize with synthetic profile
        synthetic_profile = SpeciesProfile.load(test['synthetic_profile'])
        result_synthetic = xeno.xenotypize(
            trastuzumab_VL,
            species,
            'Trastuzumab_VL',
            test['chain_type'],
            species_profile=synthetic_profile,
        )
        synthetic_score = scorer.score_sequence(result_synthetic.xenotypized_sequence)
        print(f"\n3. Synthetic profile (500 seqs):")
        print(f"   Xenotypized T20: {synthetic_score:.2f}")
        print(f"   Improvement: {synthetic_score - orig_score:.2f}")
        print(f"   Substitutions: {result_synthetic.total_substitutions}")
        
        # 4. Compare
        print(f"\n4. Comparison:")
        print(f"   Germline improvement: {germline_score - orig_score:.2f}")
        print(f"   Synthetic improvement: {synthetic_score - orig_score:.2f}")
        print(f"   Difference: {(synthetic_score - orig_score) - (germline_score - orig_score):.2f}")


if __name__ == "__main__":
    main()
