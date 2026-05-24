#!/usr/bin/env python3
"""
CDR3 generation using IGoR/OLGA models with fallback to improved random generation.

This module provides realistic CDR3 sequence generation for antibody repertoires:
1. For human/mouse: Uses OLGA with pre-trained IGoR models
2. For veterinary species: Uses improved random generation with:
   - Realistic CDR3 length distributions (from IMGT J-gene analysis)
   - Conserved anchor residues (Cys start, J-gene-encoded C-terminal)
   - Position-specific amino acid biases from known CDR3 sequences

Usage:
    from generate_cdr3_igor import CDR3Generator

    gen = CDR3Generator(species='cat', germline_dir='data/germlines/imgt')
    cdr3 = gen.generate()
"""

import logging
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# ============================================================
# Species mapping for OLGA models
# ============================================================

# Map species to OLGA default models (if available)
OLGA_MODELS = {
    'human': {
        'kappa': 'human_B_kappa',
        'lambda': 'human_B_lambda',
        'heavy': 'human_B_heavy',
    },
    'mouse': {
        'kappa': 'mouse_B_kappa',
        'lambda': 'mouse_B_lambda',
        'heavy': 'mouse_B_heavy',
    },
}

# Map species to IMGT directory names
SPECIES_DIR_MAP = {
    'cat': 'Felis_catus',
    'goat': 'Capra_hircus',
    'dog': 'Canis_lupus_familiaris',
    'horse': 'Equus_caballus',
    'cattle': 'Bos_taurus',
    'pig': 'Sus_scrofa',
    'sheep': 'Ovis_aries',
    'rabbit': 'Oryctolagus_cuniculus',
    'alpaca': 'Vicugna_pacos',
    'llama': 'Vicugna_pacos',
    'human': 'Homo_sapiens',
    'mouse': 'Mus_musculus',
}


# ============================================================
# Realistic CDR3 length distributions
# ============================================================

# These are based on analysis of known antibody sequences
# Source: IMGT, OAS (Observed Antibody Space), and published literature
CDR3_LENGTH_DISTRIBUTIONS = {
    # Light chain (kappa + lambda) CDR3 lengths
    'light': {
        # Most common: 9-11 residues for light chains
        9: 0.15,
        10: 0.25,
        11: 0.20,
        8: 0.10,
        12: 0.10,
        7: 0.05,
        13: 0.05,
        6: 0.03,
        14: 0.03,
        5: 0.02,
        15: 0.02,
    },
    # Heavy chain CDR3 lengths
    'heavy': {
        # Most common: 12-16 residues for heavy chains
        12: 0.10,
        13: 0.12,
        14: 0.15,
        15: 0.15,
        16: 0.12,
        17: 0.08,
        18: 0.06,
        11: 0.06,
        19: 0.04,
        20: 0.03,
        10: 0.03,
        21: 0.02,
        9: 0.02,
        22: 0.01,
        8: 0.01,
    },
}


# ============================================================
# Position-specific amino acid biases for CDR3
# ============================================================

# Based on analysis of known CDR3 sequences
# Position 0 is always Cys (from V-gene), last position is from J-gene

# Amino acid frequencies for internal CDR3 positions (excluding anchors)
# These reflect the observed diversity in real CDR3 loops
CDR3_INTERNAL_AA_FREQ = {
    'Y': 0.12, 'G': 0.11, 'S': 0.10, 'T': 0.09, 'D': 0.08,
    'N': 0.07, 'A': 0.06, 'L': 0.05, 'F': 0.04, 'R': 0.04,
    'H': 0.03, 'P': 0.03, 'I': 0.03, 'V': 0.03, 'K': 0.03,
    'E': 0.02, 'W': 0.02, 'Q': 0.02, 'M': 0.01, 'C': 0.01,
}

# Position-specific biases (relative to CDR3 start)
# Position 1 (after Cys): often Tyr, Phe, or other aromatic
CDR3_POS1_BIAS = {
    'Y': 0.25, 'F': 0.15, 'S': 0.10, 'G': 0.10, 'D': 0.08,
    'N': 0.06, 'T': 0.06, 'A': 0.05, 'H': 0.04, 'R': 0.04,
    'L': 0.03, 'I': 0.02, 'V': 0.02,
}

# Middle positions: more diverse
CDR3_MIDDLE_AA_FREQ = CDR3_INTERNAL_AA_FREQ.copy()

# Position -2 (before J-gene anchor): often Gly, Ser, or Ala
CDR3_PRE_ANCHOR_BIAS = {
    'G': 0.20, 'S': 0.15, 'A': 0.10, 'T': 0.08, 'D': 0.08,
    'N': 0.06, 'Y': 0.06, 'R': 0.05, 'L': 0.04, 'F': 0.04,
    'H': 0.03, 'P': 0.03, 'I': 0.03, 'V': 0.03, 'K': 0.02,
    'E': 0.02, 'W': 0.01, 'Q': 0.01,
}


# ============================================================
# CDR3Generator class
# ============================================================

class CDR3Generator:
    """
    Generate realistic CDR3 sequences using IGoR/OLGA or improved random generation.

    For human/mouse: Uses OLGA with pre-trained IGoR models
    For veterinary species: Uses improved random generation with realistic biases
    """

    def __init__(
        self,
        species: str,
        germline_dir: str = 'data/germlines/imgt',
        chain_type: str = 'light',
        olga_model_dir: Optional[str] = None,
    ):
        """
        Initialize the CDR3 generator.

        Args:
            species: Species name (e.g., 'cat', 'goat', 'human')
            germline_dir: Path to IMGT germline directory
            chain_type: 'light' or 'heavy'
            olga_model_dir: Path to custom OLGA model directory (optional)
        """
        self.species = species.lower()
        self.germline_dir = Path(germline_dir)
        self.chain_type = chain_type
        self.olga_model_dir = olga_model_dir

        # Check if OLGA model is available
        self.use_olga = False
        self.olga_model = None

        if self.species in OLGA_MODELS:
            # Try to load OLGA model
            try:
                self._load_olga_model()
                self.use_olga = True
                logger.info(f"Using OLGA model for {self.species} {chain_type}")
            except Exception as e:
                logger.warning(f"Could not load OLGA model for {self.species}: {e}")
                logger.info("Falling back to improved random generation")

        if not self.use_olga:
            logger.info(f"Using improved random CDR3 generation for {self.species}")

        # Load J-gene sequences for anchor residues
        self._load_j_genes()

    def _load_olga_model(self):
        """Load OLGA model for the species."""
        import olga.load_model as load_model
        import olga.sequence_generation as seq_gen

        # Determine model folder
        if self.olga_model_dir:
            model_folder = self.olga_model_dir
        else:
            import olga
            olga_base = Path(olga.__file__).parent / 'default_models'

            if self.chain_type == 'heavy':
                model_name = OLGA_MODELS[self.species].get('heavy', 'human_B_heavy')
            elif self.chain_type == 'light':
                # Try kappa first, then lambda
                model_name = OLGA_MODELS[self.species].get('kappa', 'human_B_kappa')
            else:
                model_name = OLGA_MODELS[self.species].get('kappa', 'human_B_kappa')

            model_folder = str(olga_base / model_name)

        # Load model files
        params_file = Path(model_folder) / 'model_params.txt'
        marginals_file = Path(model_folder) / 'model_marginals.txt'
        V_anchor_file = Path(model_folder) / 'V_gene_CDR3_anchors.csv'
        J_anchor_file = Path(model_folder) / 'J_gene_CDR3_anchors.csv'

        # Check files exist
        for f in [params_file, marginals_file, V_anchor_file, J_anchor_file]:
            if not f.exists():
                raise FileNotFoundError(f"OLGA model file not found: {f}")

        # Load genomic data and generative model
        if self.chain_type == 'heavy':
            self.olga_genomic_data = load_model.GenomicDataVDJ()
            self.olga_genomic_data.load_igor_genomic_data(
                str(params_file), str(V_anchor_file), str(J_anchor_file)
            )
            self.olga_generative_model = load_model.GenerativeModelVDJ()
            self.olga_generative_model.load_and_process_igor_model(str(marginals_file))
            self.olga_seq_gen = seq_gen.SequenceGenerationVDJ(
                self.olga_generative_model, self.olga_genomic_data
            )
        else:
            self.olga_genomic_data = load_model.GenomicDataVJ()
            self.olga_genomic_data.load_igor_genomic_data(
                str(params_file), str(V_anchor_file), str(J_anchor_file)
            )
            self.olga_generative_model = load_model.GenerativeModelVJ()
            self.olga_generative_model.load_and_process_igor_model(str(marginals_file))
            self.olga_seq_gen = seq_gen.SequenceGenerationVJ(
                self.olga_generative_model, self.olga_genomic_data
            )

    def _load_j_genes(self):
        """Load J-gene sequences for anchor residue information."""
        species_dir = SPECIES_DIR_MAP.get(self.species)
        if species_dir is None:
            self.j_gene_anchors = {}
            return

        imgt_dir = self.germline_dir / species_dir
        if not imgt_dir.exists():
            self.j_gene_anchors = {}
            return

        self.j_gene_anchors = {}

        # Load J-genes from FASTA files
        for j_file in ['IGKJ.fasta', 'IGLJ.fasta', 'IGHJ.fasta']:
            j_path = imgt_dir / j_file
            if j_path.exists() and j_path.stat().st_size > 0:
                genes = self._parse_j_gene_fasta(str(j_path))
                self.j_gene_anchors.update(genes)

    def _parse_j_gene_fasta(self, fasta_path: str) -> Dict[str, str]:
        """
        Parse J-gene FASTA file to extract N-terminal anchor residues.

        The anchor residue is the FIRST amino acid of the J-gene protein sequence
        (typically F, W, or V), which is the conserved residue at the end of CDR3.

        Returns:
            Dict of {gene_name: anchor_amino_acid}
        """
        genes = {}
        current_name = None
        current_nuc = []

        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name is not None and current_nuc:
                        nuc_seq = ''.join(current_nuc).replace('.', '').replace('-', '')
                        if len(nuc_seq) >= 4:  # Need at least 4 nt (1 offset + 3 for codon)
                            # Skip first nucleotide (reading frame offset)
                            nuc_seq = nuc_seq[1:]
                            # Translate first codon to get anchor residue
                            first_codon = nuc_seq[:3].upper()
                            aa = self._translate_codon(first_codon)
                            if aa and aa != '*':
                                genes[current_name] = aa
                    current_name = line[1:].split()[0]
                    current_nuc = []
                else:
                    current_nuc.append(line)

        # Process last gene
        if current_name is not None and current_nuc:
            nuc_seq = ''.join(current_nuc).replace('.', '').replace('-', '')
            if len(nuc_seq) >= 4:  # Need at least 4 nt (1 offset + 3 for codon)
                # Skip first nucleotide (reading frame offset)
                nuc_seq = nuc_seq[1:]
                # Translate first codon to get anchor residue
                first_codon = nuc_seq[:3].upper()
                aa = self._translate_codon(first_codon)
                if aa and aa != '*':
                    genes[current_name] = aa

        return genes

    def _translate_codon(self, codon: str) -> Optional[str]:
        """Translate a single codon to amino acid."""
        codon_table = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
        }
        return codon_table.get(codon.upper(), 'X')

    def generate(
        self,
        j_gene_name: Optional[str] = None,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
    ) -> str:
        """
        Generate a single CDR3 sequence.

        Args:
            j_gene_name: Name of J-gene (for anchor residue)
            min_length: Minimum CDR3 length (overrides distribution)
            max_length: Maximum CDR3 length (overrides distribution)

        Returns:
            CDR3 amino acid sequence (without leading Cys)
        """
        if self.use_olga:
            return self._generate_olga()
        else:
            return self._generate_improved_random(
                j_gene_name=j_gene_name,
                min_length=min_length,
                max_length=max_length,
            )

    def generate_batch(
        self,
        n: int,
        j_gene_names: Optional[List[str]] = None,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
    ) -> List[str]:
        """
        Generate a batch of CDR3 sequences.

        Args:
            n: Number of sequences to generate
            j_gene_names: List of J-gene names (for anchor residues)
            min_length: Minimum CDR3 length
            max_length: Maximum CDR3 length

        Returns:
            List of CDR3 amino acid sequences
        """
        cdr3s = []
        for i in range(n):
            j_gene = j_gene_names[i] if j_gene_names and i < len(j_gene_names) else None
            cdr3 = self.generate(
                j_gene_name=j_gene,
                min_length=min_length,
                max_length=max_length,
            )
            cdr3s.append(cdr3)
        return cdr3s

    def _generate_olga(self) -> str:
        """
        Generate CDR3 using OLGA.

        Returns:
            CDR3 amino acid sequence (without leading Cys)
        """
        while True:
            try:
                ntseq, aaseq, v_idx, j_idx = self.olga_seq_gen.gen_rnd_prod_CDR3()
                # OLGA returns full CDR3 including leading Cys
                # Remove leading Cys for consistency
                if aaseq.startswith('C'):
                    return aaseq[1:]
                return aaseq
            except Exception as e:
                logger.warning(f"OLGA generation failed: {e}")
                # Fall back to random generation
                return self._generate_improved_random()

    def _generate_improved_random(
        self,
        j_gene_name: Optional[str] = None,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
    ) -> str:
        """
        Generate CDR3 using improved random generation with realistic biases.

        Args:
            j_gene_name: Name of J-gene (for anchor residue)
            min_length: Minimum CDR3 length
            max_length: Maximum CDR3 length

        Returns:
            CDR3 amino acid sequence (without leading Cys)
        """
        # Determine CDR3 length
        if min_length is not None and max_length is not None:
            length = random.randint(min_length, max_length)
        else:
            length = self._sample_cdr3_length()

        # Generate internal residues (excluding Cys and J-gene anchor)
        internal_length = length - 1  # Subtract 1 for the anchor residue from J-gene

        if internal_length <= 0:
            # Very short CDR3, just return anchor
            return self._get_j_anchor(j_gene_name)

        # Build CDR3 sequence
        cdr3 = []

        # Position 1 (after Cys): use position-specific bias
        if internal_length >= 1:
            cdr3.append(self._sample_aa(CDR3_POS1_BIAS))

        # Middle positions
        for i in range(1, internal_length - 1):
            if i == internal_length - 2:
                # Pre-anchor position
                cdr3.append(self._sample_aa(CDR3_PRE_ANCHOR_BIAS))
            else:
                cdr3.append(self._sample_aa(CDR3_INTERNAL_AA_FREQ))

        # Last internal position (before J-gene anchor)
        if internal_length >= 2:
            cdr3.append(self._sample_aa(CDR3_PRE_ANCHOR_BIAS))

        # Add J-gene anchor residue
        anchor = self._get_j_anchor(j_gene_name)
        cdr3.append(anchor)

        return ''.join(cdr3)

    def _sample_cdr3_length(self) -> int:
        """
        Sample CDR3 length from realistic distribution.

        Returns:
            CDR3 length
        """
        dist = CDR3_LENGTH_DISTRIBUTIONS.get(self.chain_type, CDR3_LENGTH_DISTRIBUTIONS['light'])
        lengths = list(dist.keys())
        weights = list(dist.values())
        return random.choices(lengths, weights=weights, k=1)[0]

    def _sample_aa(self, freq_dict: Dict[str, float]) -> str:
        """
        Sample an amino acid from a frequency distribution.

        Args:
            freq_dict: Dictionary of {amino_acid: frequency}

        Returns:
            Sampled amino acid
        """
        aas = list(freq_dict.keys())
        weights = list(freq_dict.values())
        return random.choices(aas, weights=weights, k=1)[0]

    def _get_j_anchor(self, j_gene_name: Optional[str] = None) -> str:
        """
        Get the J-gene anchor residue.

        Args:
            j_gene_name: Name of J-gene

        Returns:
            Anchor amino acid (typically F, W, or V)
        """
        if j_gene_name and j_gene_name in self.j_gene_anchors:
            return self.j_gene_anchors[j_gene_name]

        # Default anchor residues based on chain type
        if self.chain_type == 'heavy':
            # Heavy chain: typically F or W
            return random.choices(['F', 'W', 'V'], weights=[0.7, 0.25, 0.05], k=1)[0]
        else:
            # Light chain: typically F
            return random.choices(['F', 'W', 'V'], weights=[0.85, 0.10, 0.05], k=1)[0]

    def get_cdr3_length_distribution(self) -> Dict[int, float]:
        """
        Get the CDR3 length distribution for this chain type.

        Returns:
            Dictionary of {length: probability}
        """
        return CDR3_LENGTH_DISTRIBUTIONS.get(
            self.chain_type,
            CDR3_LENGTH_DISTRIBUTIONS['light']
        ).copy()

    def is_using_olga(self) -> bool:
        """Check if OLGA is being used for generation."""
        return self.use_olga


# ============================================================
# Convenience functions
# ============================================================

def generate_cdr3(
    species: str,
    germline_dir: str = 'data/germlines/imgt',
    chain_type: str = 'light',
    j_gene_name: Optional[str] = None,
) -> str:
    """
    Generate a single CDR3 sequence for the given species.

    Args:
        species: Species name
        germline_dir: Path to IMGT germline directory
        chain_type: 'light' or 'heavy'
        j_gene_name: J-gene name (for anchor residue)

    Returns:
        CDR3 amino acid sequence (without leading Cys)
    """
    gen = CDR3Generator(
        species=species,
        germline_dir=germline_dir,
        chain_type=chain_type,
    )
    return gen.generate(j_gene_name=j_gene_name)


def generate_cdr3_batch(
    species: str,
    n: int,
    germline_dir: str = 'data/germlines/imgt',
    chain_type: str = 'light',
    j_gene_names: Optional[List[str]] = None,
) -> List[str]:
    """
    Generate a batch of CDR3 sequences for the given species.

    Args:
        species: Species name
        n: Number of sequences
        germline_dir: Path to IMGT germline directory
        chain_type: 'light' or 'heavy'
        j_gene_names: List of J-gene names

    Returns:
        List of CDR3 amino acid sequences
    """
    gen = CDR3Generator(
        species=species,
        germline_dir=germline_dir,
        chain_type=chain_type,
    )
    return gen.generate_batch(n=n, j_gene_names=j_gene_names)


# ============================================================
# Main (for testing)
# ============================================================

if __name__ == '__main__':
    import argparse

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    parser = argparse.ArgumentParser(description="Generate CDR3 sequences")
    parser.add_argument('--species', '-s', required=True, help='Species name')
    parser.add_argument('--num-seq', '-n', type=int, default=10, help='Number of sequences')
    parser.add_argument('--chain-type', '-c', default='light', choices=['light', 'heavy'])
    parser.add_argument('--germline-dir', default='data/germlines/imgt')
    parser.add_argument('--seed', type=int, default=42)

    args = parser.parse_args()
    random.seed(args.seed)
    np.random.seed(args.seed)

    gen = CDR3Generator(
        species=args.species,
        germline_dir=args.germline_dir,
        chain_type=args.chain_type,
    )

    print(f"Species: {args.species}")
    print(f"Chain type: {args.chain_type}")
    print(f"Using OLGA: {gen.is_using_olga()}")
    print(f"\nGenerated CDR3 sequences:")

    for i in range(args.num_seq):
        cdr3 = gen.generate()
        print(f"  {i+1:3d}. {cdr3} (len={len(cdr3)})")
