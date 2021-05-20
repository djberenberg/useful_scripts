# Useful scripts

- "General" purpose utility scripts 

# Requirements
```
pip install numpy scipy matplotlib torch Biopython tqdm
```

# What's inside
- `mkdmap.py` - make a distance map from pdb file
- `filter-fasta.py` - split and/or filter sequences by length from a fasta file
- `plot_map.py` - plots a contact map

```
usage: filter-fasta.py [-h] -i INPUT [-o OUTPUT] [-d DOMAIN_FILE]
                       [--split SPLIT] [-s] [-v]
                       [--assert ASSERTION [ASSERTION ...]]

Split or filter a

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT              Input filename
  -o OUTPUT             Output fasta file
  -d DOMAIN_FILE        Filter by provided domains
  --split SPLIT         sequence header delimiter to split on
  -s, --include-stops   Include sequences with stop codons
  -v, --verbose         Verbose output
  --assert ASSERTION [ASSERTION ...]
                        Condition for sequences of the form '[>|<|>=|<=]\d+'
```
