#!/usr/bin/env python3
"""
Test IGoR/OLGA-based synthetic profiles vs old synthetic profiles.

Compares T20 scores between old and new synthetic profiles.
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
    print("Testing IGoR/OLGA-based Synthetic Profiles")
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
            'old_synthetic_profile': str(PANIG_DIR / 'profiles' / 'cat_VL_synthetic.json'),
            'new_synthetic_profile': str(PANIG_DIR / 'profiles' / 'cat_VL_synthetic_igor.json'),
        },
        {
            'species': 'goat',
            'chain_type': 'light',
            'blastdb': str(PANIG_DIR / 'blastdb' / 'goat_VL_blastdb' / 'goat_VL'),
            'germline_profile': str(PANIG_DIR / 'profiles' / 'goat_VL.json'),
            'old_synthetic_profile': str(PANIG_DIR / 'profiles' / 'goat_VL_synthetic.json'),
            'new_synthetic_profile': str(PANIG_DIR / 'profiles' / 'goat_VL_synthetic_igor.json'),
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
        
        # 3. Xenotypize with old synthetic profile
        old_synthetic_profile = SpeciesProfile.load(test['old_synthetic_profile'])
        result_old_synthetic = xeno.xenotypize(
            trastuzumab_VL,
            species,
            'Trastuzumab_VL',
            test['chain_type'],
            species_profile=old_synthetic_profile,
        )
        old_synthetic_score = scorer.score_sequence(result_old_synthetic.xenotypized_sequence)
        print(f"\n3. Old synthetic profile (random CDR3):")
        print(f"   Xenotypized T20: {old_synthetic_score:.2f}")
        print(f"   Improvement: {old_synthetic_score - orig_score:.2f}")
        print(f"   Substitutions: {result_old_synthetic.total_substitutions}")
        
        # 4. Xenotypize with new synthetic profile (IGoR/OLGA)
        new_synthetic_profile = SpeciesProfile.load(test['new_synthetic_profile'])
        result_new_synthetic = xeno.xenotypize(
            trastuzumab_VL,
            species,
            'Trastuzumab_VL',
            test['chain_type'],
            species_profile=new_synthetic_profile,
        )
        new_synthetic_score = scorer.score_sequence(result_new_synthetic.xenotypized_sequence)
        print(f"\n4. New synthetic profile (IGoR/OLGA CDR3):")
        print(f"   Xenotypized T20: {new_synthetic_score:.2f}")
        print(f"   Improvement: {new_synthetic_score - orig_score:.2f}")
        print(f"   Substitutions: {result_new_synthetic.total_substitutions}")
        
        # 5. Compare
        print(f"\n5. Comparison:")
        print(f"   Germline improvement: {germline_score - orig_score:.2f}")
        print(f"   Old synthetic improvement: {old_synthetic_score - orig_score:.2f}")
        print(f"   New synthetic improvement: {new_synthetic_score - orig_score:.2f}")
        print(f"   New vs Old difference: {(new_synthetic_score - orig_score) - (old_synthetic_score - orig_score):.2f}")


if __name__ == "__main__":
    main()
