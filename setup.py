"""Setup script for PanIg."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="panig",
    version="0.1.0",
    author="PanIg Contributors",
    description="Pan-species Immunoglobulin Xenotypization Tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Lefrunila/PanIg",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "anarcii>=2.0.0",
        "immunebuilder>=1.1.0",
        "biopython>=1.80",
        "numpy>=1.24",
        "pandas>=2.0",
        "torch>=2.0",
        "gdown>=4.7",
    ],
    extras_require={
        "interactions": ["protinter>=1.0.0"],
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black",
            "flake8",
        ],
    },
    entry_points={
        "console_scripts": [
            "panig=panig.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
