"""
Configuration management for PanIg.

Handles paths, Google Drive integration, and default settings.
"""

import os
from pathlib import Path


# Project root directory
PROJECT_ROOT = Path(__file__).parent.parent

# Default directories
PROFILES_DIR = PROJECT_ROOT / "profiles"
BLASTDB_DIR = PROJECT_ROOT / "blastdb"
TESTS_DIR = PROJECT_ROOT / "tests"
EXAMPLES_DIR = PROJECT_ROOT / "examples"

# Google Drive configuration
GDRIVE_REMOTE = "gdrive"
GDRIVE_DATABASE_FOLDER = "PanIg_databases"

# Local cache directory
CACHE_DIR = Path.home() / ".panig" / "cache"
CACHE_PROFILES_DIR = CACHE_DIR / "profiles"
CACHE_BLASTDB_DIR = CACHE_DIR / "blastdb"

# Default settings
DEFAULT_SCHEME = "imgt"
DEFAULT_THRESHOLD = 0.1
DEFAULT_CHAIN_TYPE = "VH"  # Heavy chain by default

# Supported species
SUPPORTED_SPECIES = [
    "human",
    "dog",
    "cat",
    "horse",
    "cattle",
    "pig",
    "sheep",
    "goat",
    "rabbit",
    "hamster",
    "mouse",
    "rat",
]

# Supported chain types
SUPPORTED_CHAIN_TYPES = ["VH", "VHH", "VL"]


def get_profile_path(species: str, chain_type: str = "VH") -> Path:
    """
    Get the path to a species profile file.

    Searches in order:
    1. Local profiles directory
    2. Cache directory

    Args:
        species: Species name
        chain_type: Chain type ('VH' or 'VHH')

    Returns:
        Path to profile file
    """
    filename = f"{species}_{chain_type}.json"

    # Check local first
    local_path = PROFILES_DIR / filename
    if local_path.exists():
        return local_path

    # Check cache
    cache_path = CACHE_PROFILES_DIR / filename
    if cache_path.exists():
        return cache_path

    # Return default location
    return local_path


def get_blastdb_path(species: str, chain_type: str = "VH") -> Path:
    """
    Get the path to a BLAST database.

    Args:
        species: Species name
        chain_type: Chain type

    Returns:
        Path to BLAST database directory
    """
    db_name = f"{species}_{chain_type}_blastdb"

    # Check local first
    local_path = BLASTDB_DIR / db_name
    if local_path.exists():
        return local_path

    # Check cache
    cache_path = CACHE_BLASTDB_DIR / db_name
    if cache_path.exists():
        return cache_path

    return local_path


def ensure_cache_dirs():
    """Create cache directories if they don't exist."""
    CACHE_PROFILES_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_BLASTDB_DIR.mkdir(parents=True, exist_ok=True)


def get_gdrive_remote_path(subpath: str) -> str:
    """
    Get the full Google Drive remote path.

    Args:
        subpath: Sub-path within the database folder

    Returns:
        Full rclone remote path
    """
    return f"{GDRIVE_REMOTE}:{GDRIVE_DATABASE_FOLDER}/{subpath}"
