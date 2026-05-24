#!/usr/bin/env python3
"""
Generate synthetic antibody VL repertoires from IMGT germline V(D)J recombination + SHM.

This addresses the problem that germline-only profiles (40-55 genes) don't capture
the diversity of mature B-cell repertoires. By programmatically simulating V(D)J
recombination and somatic hypermutation, we generate 500 representative sequences
per species that better reflect circulating antibodies.

Key insight: IMGT germline FASTA files contain aligned nucleotide sequences with gaps.
We must preserve the alignment structure to correctly map to IMGT positions.

Usage:
    python generate_synthetic_repertoire.py --species cat --num-seq 500 --output data/synthetic/cat_VL_synthetic.fasta
"""

import argparse
import json
import logging
import random
import sys
from collections import defaultdict
from pathlib import Path

# Add scripts directory to path for importing CDR3Generator
sys.path.insert(0, str(Path(__file__).parent))
from generate_cdr3_igor import CDR3Generator

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)


# ============================================================
# Constants
# ============================================================

STANDARD_CODON_TABLE = {
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

# CDR3 amino acid composition (enriched in real CDR3 loops)
CDR3_AA_WEIGHTS = {
    'Y': 15, 'G': 14, 'S': 13, 'T': 11, 'D': 10, 'N': 9,
    'A': 7, 'L': 6, 'F': 5, 'R': 5, 'H': 4, 'P': 4,
    'I': 3, 'V': 3, 'K': 3, 'E': 3, 'W': 2, 'C': 2,
    'Q': 2, 'M': 1,
}

# AID hotspot motifs: WRCY and RGYW
# W = A/T, R = A/G, Y = C/T
WRCY_MOTIFS = set()
RGYW_MOTIFS = set()
for w in 'AT':
    for r in 'AG':
        for y in 'CT':
            WRCY_MOTIFS.add(f'{w}{r}C{y}')
            RGYW_MOTIFS.add(f'{r}G{y}{w}')

# AID transition mutations: C->T, G->A (in context of WRCY/RGYW)
AID_TRANSITIONS = {'C': 'T', 'G': 'A', 'A': 'G', 'T': 'C'}

# Species name to IMGT directory mapping
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
# IMGT alignment-aware parsing
# ============================================================

def translate_codon(codon: str) -> str:
    """Translate a single codon to amino acid."""
    return STANDARD_CODON_TABLE.get(codon.upper(), 'X')


def parse_germline_v_genes_aligned(fasta_path: str) -> dict:
    """
    Parse IMGT germline V-gene FASTA file preserving IMGT alignment.
    
    IMGT germline files contain nucleotide sequences aligned to IMGT positions.
    Gaps are represented as dots. We need to:
    1. Walk through the alignment 3 nucleotides at a time (codon by codon)
    2. For each codon, track which IMGT position it corresponds to
    3. Skip codons that contain gaps
    4. Translate valid codons to amino acids
    
    Returns:
        Dict of {gene_name: {
            'protein': protein_sequence,
            'imgt_positions': list of (imgt_pos, amino_acid) tuples,
            'last_imgt_pos': last IMGT position covered
        }}
    """
    genes = {}
    current_name = None
    current_nuc_aligned = []
    
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    # Process the previous gene
                    nuc_aligned = "".join(current_nuc_aligned)
                    gene_data = _process_v_gene_alignment(nuc_aligned)
                    if gene_data:
                        genes[current_name] = gene_data
                current_name = line[1:].split()[0]
                current_nuc_aligned = []
            else:
                current_nuc_aligned.append(line)
    
    if current_name is not None:
        nuc_aligned = "".join(current_nuc_aligned)
        gene_data = _process_v_gene_alignment(nuc_aligned)
        if gene_data:
            genes[current_name] = gene_data
    
    return genes


def _process_v_gene_alignment(nuc_aligned: str) -> dict:
    """
    Process an aligned V-gene nucleotide sequence to extract protein and IMGT positions.
    
    The alignment has dots (.) for gaps. We need to:
    1. Walk through the alignment 3 nucleotides at a time (codon by codon)
    2. For each codon, track which IMGT position it corresponds to
    3. Skip codons that contain gaps
    4. Translate valid codons to amino acids
    """
    # Remove newlines and spaces
    nuc_aligned = nuc_aligned.replace('\n', '').replace(' ', '')
    
    # Walk through the alignment
    protein = []
    imgt_positions = []
    imgt_pos = 1  # IMGT positions start at 1
    
    i = 0
    while i < len(nuc_aligned) - 2:
        codon = nuc_aligned[i:i+3].upper()
        
        # Check if codon contains a gap
        if '.' in codon or '-' in codon:
            # This codon has a gap - skip it but don't advance IMGT position
            i += 3
            continue
        
        # Check if it's a valid DNA codon
        if all(c in 'ACGT' for c in codon):
            aa = translate_codon(codon)
            if aa == '*':
                break  # Stop codon
            protein.append(aa)
            imgt_positions.append((imgt_pos, aa))
            imgt_pos += 1
        
        i += 3
    
    if len(protein) < 50:
        return None
    
    return {
        'protein': ''.join(protein),
        'imgt_positions': imgt_positions,
        'last_imgt_pos': imgt_pos - 1,
    }


def parse_germline_j_genes_aligned(fasta_path: str) -> dict:
    """
    Parse IMGT germline J-gene FASTA file preserving IMGT alignment.
    
    J-genes cover IMGT positions 118-128 (FR4).
    The J-gene nucleotide sequence has a leading nucleotide before the reading frame.
    
    Returns:
        Dict of {gene_name: {
            'protein': protein_sequence,
            'imgt_positions': list of (imgt_pos, amino_acid) tuples
        }}
    """
    genes = {}
    current_name = None
    current_nuc_aligned = []
    
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    # Process the previous gene
                    nuc_aligned = "".join(current_nuc_aligned)
                    gene_data = _process_j_gene_alignment(nuc_aligned)
                    if gene_data:
                        genes[current_name] = gene_data
                current_name = line[1:].split()[0]
                current_nuc_aligned = []
            else:
                current_nuc_aligned.append(line)
    
    if current_name is not None:
        nuc_aligned = "".join(current_nuc_aligned)
        gene_data = _process_j_gene_alignment(nuc_aligned)
        if gene_data:
            genes[current_name] = gene_data
    
    return genes


def _process_j_gene_alignment(nuc_aligned: str) -> dict:
    """
    Process an aligned J-gene nucleotide sequence.
    
    J-genes have a leading nucleotide before the reading frame starts.
    We need to skip the first nucleotide to get the correct reading frame.
    """
    # Remove newlines and spaces
    nuc_aligned = nuc_aligned.replace('\n', '').replace(' ', '')
    
    # Skip the first nucleotide (reading frame offset)
    nuc_aligned = nuc_aligned[1:]
    
    # Walk through the alignment
    protein = []
    imgt_positions = []
    imgt_pos = 118  # J-gene starts at IMGT position 118
    
    i = 0
    while i < len(nuc_aligned) - 2:
        codon = nuc_aligned[i:i+3].upper()
        
        # Check if codon contains a gap
        if '.' in codon or '-' in codon:
            i += 3
            continue
        
        # Check if it's a valid DNA codon
        if all(c in 'ACGT' for c in codon):
            aa = translate_codon(codon)
            if aa == '*':
                break  # Stop codon
            protein.append(aa)
            imgt_positions.append((imgt_pos, aa))
            imgt_pos += 1
        
        i += 3
    
    if len(protein) < 5:
        return None
    
    return {
        'protein': ''.join(protein),
        'imgt_positions': imgt_positions,
    }


# ============================================================
# CDR3 generation
# ============================================================

# Cache for CDR3Generator instances to avoid recreating them
_cdr3_generator_cache = {}
_cdr3_rng = random.Random()  # Separate RNG for CDR3 generation

def generate_cdr3(min_len: int = 5, max_len: int = 15, species: str = 'human', j_gene_name: str = None) -> str:
    """
    Generate a CDR3 loop sequence using IGoR/OLGA or improved random generation.
    
    Uses a separate random state to avoid perturbing the main sequence generation.
    
    For human/mouse: Uses OLGA with pre-trained IGoR models
    For other species: Uses improved random generation with realistic biases
    
    Args:
        min_len: Minimum CDR3 length (default 5)
        max_len: Maximum CDR3 length (default 15)
        species: Species name for OLGA model selection
        j_gene_name: J-gene name for anchor residue
    
    Returns:
        CDR3 amino acid sequence (without leading Cys)
    """
    # Use cached generator if available
    cache_key = f"{species}_light"
    if cache_key not in _cdr3_generator_cache:
        gen = CDR3Generator(
            species=species,
            germline_dir='data/germlines/imgt',
            chain_type='light',
        )
        _cdr3_generator_cache[cache_key] = gen
    gen = _cdr3_generator_cache[cache_key]
    
    # Use separate RNG for CDR3 generation to avoid perturbing main random state
    old_state = random.getstate()
    random.setstate(_cdr3_rng.getstate())
    try:
        result = gen.generate(j_gene_name=j_gene_name, min_length=min_len, max_length=max_len)
    finally:
        _cdr3_rng.setstate(random.getstate())
        random.setstate(old_state)
    
    return result


# ============================================================
# Somatic hypermutation (SHM)
# ============================================================

def has_aid_motif(seq: str, pos: int, motif_len: int = 4) -> bool:
    """
    Check if position is within an AID hotspot motif (WRCY or RGYW).
    
    AID (Activation-Induced Cytidine Deaminase) preferentially targets
    C residues in WRCY/RGYW motifs, causing C->T transitions.
    """
    for start in range(max(0, pos - motif_len + 1), min(len(seq) - motif_len + 1, pos + 1)):
        window = seq[start:start + motif_len].upper()
        if window in WRCY_MOTIFS or window in RGYW_MOTIFS:
            return True
    return False


def apply_shm(
    sequence: str,
    v_gene_len: int,
    cdr3_len: int,
    base_rate: float = 0.03,
    hotspot_boost: float = 3.0,
) -> str:
    """
    Apply somatic hypermutation to an amino acid sequence.
    
    SHM introduces point mutations biased toward AID hotspot motifs.
    Framework positions get lower mutation rates than CDR positions.
    
    Args:
        sequence: Full amino acid sequence (V + CDR3 + J)
        v_gene_len: Length of V-gene portion (FR1-CDR1-FR2-CDR2-FR3)
        cdr3_len: Length of CDR3 portion
        base_rate: Base mutation probability per position (2-5%)
        hotspot_boost: Multiplier for mutation rate in AID hotspots
    
    Returns:
        Mutated amino acid sequence
    """
    mutated = list(sequence)
    
    # Find conserved positions in the V-gene to protect
    conserved_indices = set()
    
    # First Cys is always at position 23 in the V-gene (0-indexed 22)
    if len(sequence) > 22 and sequence[22] == 'C':
        conserved_indices.add(22)
    
    # Second Cys: find the last Cys in the V-gene portion
    v_portion = sequence[:v_gene_len]
    last_cys_in_v = v_portion.rfind('C')
    if last_cys_in_v > 0:
        conserved_indices.add(last_cys_in_v)
    
    # Trp at position ~41 in V-gene (0-indexed 40)
    if len(sequence) > 40 and sequence[40] in ('W', 'Y', 'G', 'F'):
        conserved_indices.add(40)
    
    for i, aa in enumerate(sequence):
        # Skip conserved positions
        if i in conserved_indices:
            continue
        
        # Determine mutation rate based on region
        if i < v_gene_len:
            # In V-gene: framework vs CDR
            # CDR1 is roughly positions 27-38 (0-indexed 26-37)
            # CDR2 is roughly positions 56-65 (0-indexed 55-64)
            if (26 <= i <= 37) or (55 <= i <= 64):
                rate = base_rate * 1.5  # CDR: elevated
            else:
                rate = base_rate  # Framework: base rate
        elif i < v_gene_len + cdr3_len:
            # CDR3: highest mutation rate
            rate = base_rate * 2.0
        else:
            # FR4: base rate
            rate = base_rate
        
        # Boost rate if in AID hotspot motif
        if has_aid_motif(sequence, i):
            rate *= hotspot_boost
        
        # Apply mutation
        if random.random() < rate:
            # Try AID-biased transition first
            if aa in AID_TRANSITIONS and random.random() < 0.6:
                new_aa = AID_TRANSITIONS[aa]
            else:
                # Random substitution (excluding current AA)
                other_aas = [a for a in 'ACDEFGHIKLMNPQRSTVWY' if a != aa]
                new_aa = random.choice(other_aas)
            
            mutated[i] = new_aa
    
    return ''.join(mutated)


# ============================================================
# Sequence validation
# ============================================================

def validate_sequence(seq: str, v_gene_len: int) -> bool:
    """
    Validate a synthetic antibody sequence.
    
    Checks:
    - No stop codons
    - Correct length (80-130 aa for light chain)
    - Conserved Cys at start of V-gene (position 23)
    - Conserved Cys at end of FR3 (in V-gene)
    - Second Cys present somewhere after first Cys
    """
    # No stop codons
    if '*' in seq:
        return False
    
    # Correct length
    if not (80 <= len(seq) <= 130):
        return False
    
    # First Cys at position 23 (0-indexed 22)
    if len(seq) > 22 and seq[22] != 'C':
        return False
    
    # Second Cys: find the last Cys in the V-gene portion
    v_portion = seq[:v_gene_len]
    last_cys_in_v = v_portion.rfind('C')
    if last_cys_in_v <= 22:
        return False  # No second Cys found
    
    return True


# ============================================================
# Main repertoire generation
# ============================================================

def generate_synthetic_repertoire(
    species: str,
    germline_dir: str,
    num_sequences: int = 500,
    seed: int = 42,
) -> tuple:
    """
    Generate a synthetic antibody VL repertoire for a species.
    
    Args:
        species: Species name (e.g., 'cat', 'goat')
        germline_dir: Path to IMGT germline directory
        num_sequences: Number of synthetic sequences to generate
        seed: Random seed for reproducibility
    
    Returns:
        Tuple of (sequences, imgt_position_mapping) where:
        - sequences: List of (name, sequence) tuples
        - imgt_position_mapping: Dict mapping sequence index to list of (imgt_pos, aa) tuples
    """
    random.seed(seed)
    
    species_dir = SPECIES_DIR_MAP.get(species.lower())
    if species_dir is None:
        raise ValueError(f"Unknown species: {species}. Available: {list(SPECIES_DIR_MAP.keys())}")
    
    imgt_dir = Path(germline_dir) / species_dir
    if not imgt_dir.exists():
        raise FileNotFoundError(f"Germline directory not found: {imgt_dir}")
    
    # Load V-genes (both kappa and lambda) with alignment
    v_genes = {}
    for v_file in ['IGKV.fasta', 'IGLV.fasta']:
        v_path = imgt_dir / v_file
        if v_path.exists() and v_path.stat().st_size > 0:
            loaded = parse_germline_v_genes_aligned(str(v_path))
            logger.info(f"Loaded {len(loaded)} V-genes from {v_file}")
            v_genes.update(loaded)
    
    if not v_genes:
        raise ValueError(f"No V-genes found for {species}")
    
    # Load J-genes (both kappa and lambda) with alignment
    j_genes = {}
    for j_file in ['IGKJ.fasta', 'IGLJ.fasta']:
        j_path = imgt_dir / j_file
        if j_path.exists() and j_path.stat().st_size > 0:
            loaded = parse_germline_j_genes_aligned(str(j_path))
            logger.info(f"Loaded {len(loaded)} J-genes from {j_file}")
            j_genes.update(loaded)
    
    if not j_genes:
        raise ValueError(f"No J-genes found for {species}")
    
    logger.info(f"Total: {len(v_genes)} V-genes, {len(j_genes)} J-genes")
    
    # Generate synthetic sequences
    v_gene_names = list(v_genes.keys())
    j_gene_names = list(j_genes.keys())
    
    synthetic_sequences = []
    imgt_mappings = []
    attempts = 0
    max_attempts = num_sequences * 10  # Allow 10x retries for validation failures
    
    while len(synthetic_sequences) < num_sequences and attempts < max_attempts:
        attempts += 1
        
        # 1. Pick random V and J genes
        v_name = random.choice(v_gene_names)
        j_name = random.choice(j_gene_names)
        v_data = v_genes[v_name]
        j_data = j_genes[j_name]
        
        # 2. Generate random CDR3
        # CDR3 starts after the V-gene ends and before the J-gene starts
        # V-gene covers positions 1 to v_gene.last_imgt_pos
        # J-gene covers positions 118 onwards
        # CDR3 fills the gap between v_gene.last_imgt_pos+1 and 117
        cdr3_start_pos = v_data['last_imgt_pos'] + 1
        cdr3_end_pos = 117  # CDR3 ends at position 117
        cdr3_len = random.randint(5, 15)
        
        # Generate CDR3 sequence using IGoR/OLGA or improved random generation
        cdr3 = generate_cdr3(min_len=cdr3_len, max_len=cdr3_len, species=species, j_gene_name=j_name)
        
        # 3. Concatenate: V(FR1-CDR1-FR2-CDR2-FR3) + CDR3 + J(FR4)
        full_seq = v_data['protein'] + cdr3 + j_data['protein']
        
        # 4. Apply SHM
        mutated_seq = apply_shm(
            full_seq,
            v_gene_len=len(v_data['protein']),
            cdr3_len=len(cdr3),
            base_rate=0.03,
            hotspot_boost=3.0,
        )
        
        # 5. Validate
        if validate_sequence(mutated_seq, v_gene_len=len(v_data['protein'])):
            v_short = v_name.split('|')[1] if '|' in v_name else v_name
            j_short = j_name.split('|')[1] if '|' in j_name else j_name
            seq_name = f"syn_{species}_{len(synthetic_sequences)+1:04d}|V={v_short}|J={j_short}|CDR3len={len(cdr3)}"
            synthetic_sequences.append((seq_name, mutated_seq))
            
            # Build IMGT position mapping for this sequence
            imgt_mapping = []
            
            # V-gene positions (from the alignment)
            for imgt_pos, aa in v_data['imgt_positions']:
                # Apply SHM to the amino acid at this position
                seq_idx = imgt_pos - 1  # 0-indexed
                if seq_idx < len(mutated_seq):
                    imgt_mapping.append((imgt_pos, mutated_seq[seq_idx]))
            
            # CDR3 positions (from cdr3_start_pos to 117)
            cdr3_start_idx = len(v_data['protein'])
            for i, aa in enumerate(cdr3):
                imgt_pos = cdr3_start_pos + i
                if imgt_pos <= 117:
                    seq_idx = cdr3_start_idx + i
                    if seq_idx < len(mutated_seq):
                        imgt_mapping.append((imgt_pos, mutated_seq[seq_idx]))
            
            # J-gene positions (from the alignment)
            j_start_idx = len(v_data['protein']) + len(cdr3)
            for imgt_pos, aa in j_data['imgt_positions']:
                seq_idx = j_start_idx + (imgt_pos - 118)
                if seq_idx < len(mutated_seq):
                    imgt_mapping.append((imgt_pos, mutated_seq[seq_idx]))
            
            imgt_mappings.append(imgt_mapping)
    
    if len(synthetic_sequences) < num_sequences:
        logger.warning(
            f"Only generated {len(synthetic_sequences)}/{num_sequences} "
            f"sequences after {max_attempts} attempts"
        )
    
    logger.info(f"Generated {len(synthetic_sequences)} synthetic sequences")
    return synthetic_sequences, imgt_mappings


def build_frequency_profile(imgt_mappings: list) -> dict:
    """
    Build an amino acid frequency profile from IMGT-positioned sequences.
    
    Args:
        imgt_mappings: List of lists of (imgt_position, amino_acid) tuples
    
    Returns:
        Dict of {position: {amino_acid: frequency}}
    """
    position_counts = defaultdict(lambda: defaultdict(int))
    
    for mapping in imgt_mappings:
        for imgt_pos, aa in mapping:
            if aa not in ('-', '.', 'X', '*'):
                position_counts[imgt_pos][aa] += 1
    
    # Convert counts to frequencies
    profile = {}
    for position, counts in position_counts.items():
        total = sum(counts.values())
        if total > 0:
            profile[position] = {
                aa: count / total
                for aa, count in counts.items()
            }
    
    return profile


def save_fasta(sequences: list, output_path: str):
    """Save sequences to a FASTA file."""
    with open(output_path, 'w') as f:
        for name, seq in sequences:
            f.write(f">{name}\n")
            # Write sequence in 60-char lines
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')


def save_profile_json(
    species: str,
    chain_type: str,
    profile: dict,
    num_sequences: int,
    output_path: str,
):
    """Save frequency profile as JSON (compatible with SpeciesProfile.load)."""
    data = {
        "species": species,
        "chain_type": chain_type,
        "profile": {str(k): v for k, v in profile.items()},
        "source": "synthetic_repertoire",
        "scheme": "IMGT",
        "num_synthetic_sequences": num_sequences,
    }
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic antibody VL repertoires from IMGT germline V(D)J + SHM"
    )
    parser.add_argument(
        '--species', '-s',
        required=True,
        help='Species name (e.g., cat, goat)',
    )
    parser.add_argument(
        '--num-seq', '-n',
        type=int,
        default=500,
        help='Number of synthetic sequences to generate (default: 500)',
    )
    parser.add_argument(
        '--output-fasta', '-o',
        help='Output FASTA file path (default: data/synthetic/{species}_VL_synthetic.fasta)',
    )
    parser.add_argument(
        '--output-profile', '-p',
        help='Output profile JSON path (default: profiles/{species}_VL_synthetic.json)',
    )
    parser.add_argument(
        '--germline-dir',
        default='data/germlines/imgt',
        help='Directory with IMGT germline FASTA files',
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)',
    )
    parser.add_argument(
        '--chain-type',
        default='VL',
        choices=['VL', 'VH'],
        help='Chain type (default: VL)',
    )
    
    args = parser.parse_args()
    
    # Set default paths
    if args.output_fasta is None:
        args.output_fasta = f'data/synthetic/{args.species}_VL_synthetic.fasta'
    if args.output_profile is None:
        args.output_profile = f'profiles/{args.species}_VL_synthetic.json'
    
    # Create output directories
    Path(args.output_fasta).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output_profile).parent.mkdir(parents=True, exist_ok=True)
    
    # Generate repertoire
    logger.info(f"Generating {args.num_seq} synthetic {args.species} {args.chain_type} sequences...")
    sequences, imgt_mappings = generate_synthetic_repertoire(
        species=args.species,
        germline_dir=args.germline_dir,
        num_sequences=args.num_seq,
        seed=args.seed,
    )
    
    if not sequences:
        logger.error("No sequences generated!")
        sys.exit(1)
    
    # Save FASTA
    save_fasta(sequences, args.output_fasta)
    logger.info(f"Saved {len(sequences)} sequences to {args.output_fasta}")
    
    # Build and save frequency profile
    profile = build_frequency_profile(imgt_mappings)
    save_profile_json(
        species=args.species,
        chain_type=args.chain_type,
        profile=profile,
        num_sequences=len(sequences),
        output_path=args.output_profile,
    )
    logger.info(f"Saved profile with {len(profile)} positions to {args.output_profile}")
    
    # Print summary
    lengths = [len(seq) for _, seq in sequences]
    print(f"\n{'='*60}")
    print(f"Synthetic Repertoire Generation Summary")
    print(f"{'='*60}")
    print(f"Species: {args.species}")
    print(f"Chain type: {args.chain_type}")
    print(f"Sequences generated: {len(sequences)}")
    print(f"Profile positions: {len(profile)}")
    print(f"Sequence lengths: {min(lengths)}-{max(lengths)} (mean: {sum(lengths)/len(lengths):.1f})")
    print(f"FASTA output: {args.output_fasta}")
    print(f"Profile output: {args.output_profile}")
    
    # Show some example sequences
    print(f"\nExample sequences:")
    for name, seq in sequences[:3]:
        print(f"  {name}: {seq[:50]}...")


if __name__ == "__main__":
    main()
