# BLAST Databases

This directory contains species-specific BLAST databases for T20-like scoring.

## Storage

Large BLAST databases are stored on Google Drive to conserve local disk space:
- **Remote location**: `gdrive:PanIg_databases/blastdb/`
- **Local cache**: `~/.panig/cache/blastdb/`

## Downloading

Use the PanIg CLI to download databases:

```bash
# Download dog BLAST database
panig download --species dog

# Download with specific chain type
panig download --species dog --chain-type VHH
```

## Creating Custom Databases

To create a custom BLAST database from antibody sequences:

```bash
# 1. Prepare a FASTA file with antibody sequences
# 2. Create the BLAST database
makeblastdb -in sequences.fasta -dbtype prot -out my_species_VH_blastdb

# 3. Place in this directory or cache directory
```

## Available Databases

| Species | Chain Type | Status       | Notes |
|---------|------------|--------------|-------|
| Human   | VH         | Available    | |
| Dog     | VH         | Available    | |
| Cat     | VH         | Available    | |
| Cat     | VL         | Missing      | Germline-based profile only; no expressed sequences available |
| Cattle  | VH         | Available    | |
| Cattle  | VL         | Available    | |
| Goat    | VH         | Available    | 5 sequences (very sparse) |
| Goat    | VL         | Missing      | Germline-based profile only; no expressed sequences available |
| Hamster | VH         | Available    | 31 sequences |
| Horse   | VH         | Available    | |
| Pig     | VH         | Available    | |
| Pig     | VL         | Available    | |
| Rabbit  | VH         | Available    | |
| Rabbit  | VL         | Available    | |
| Sheep   | VH         | Available    | |
| Sheep   | VL         | Available    | |
| VHH     | VHH        | Available    | |
