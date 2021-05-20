#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Split or filter a 
"""

import re
import io
import sys
import operator
import argparse
import itertools
from pathlib import Path

from tqdm import tqdm
from nard.utils import FASTAReader

clear = f"\r{100 * ' '}\r"
FASTA_STOP_CODON = '*'

def _valid_condition(cond):
    condition_structure = "(<|>|<=|>=)(\d+)" 
    match = re.match(condition_structure, cond)
    if match:
        op, operand = match.groups()
        return op, int(operand)
    else:
        raise ValueError(f"Invalid condition: {cond}")

def _construct_conditional(conditions, domain_file):
    if domain_file is not None:
        with open(domain_file, 'r') as domf:
            domains = set(map(lambda line: line.strip(), domf))
        def seq_id_is_good(seq_id):
            return seq_id in domains
    else:
        def seq_id_is_good(seq_id):
            return True

    if not conditions:
        def conditional(seq_id, sequence):
            return seq_id_is_good(seq_id)
    else:
        op2op = {'>': operator.gt, '<': operator.lt, '>=': operator.ge, '<=': operator.le} 
        operations, operands = zip(*conditions)
        operations = list(map(lambda op: op2op[op], operations))
        operands   = list(map(int, operands))
    
        def conditional(seq_id, sequence):
            evaluated = [op(len(sequence), operand) for op, operand in zip(operations, operands)]
            evaluated.append(seq_id_is_good(seq_id))
            return all(evaluated)

    return conditional

def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", help="Input filename", type=Path, metavar="INPUT", dest='input', required=True)
    parser.add_argument("-o", help="Output fasta file", type=Path, metavar="OUTPUT", dest='output')
    parser.add_argument("-d", help="Filter by provided domains", dest='domain_file') 
    parser.add_argument("--split", default=" ")
    parser.add_argument("-s", "--include-stops", help="Include sequences with stop codons",
                        default=False, action='store_true', dest='allow_stop_codons')
    parser.add_argument("-v", "--verbose", action='store_true', default=False, help="Verbose output")
    parser.add_argument("--assert", dest='assertion',
                        type=_valid_condition,
                        nargs='+',
                        help="Condition for sequences of the form '[>|<|>=|<=]\d+'",
                        default=None)

    args = parser.parse_args()
    args.condition = _construct_conditional(args.assertion, args.domain_file)
    
    args.input  = open(args.input, 'r')
    args.output = sys.stdout if args.output is None else open(args.output, 'w')

    return args

if __name__ == "__main__":
    args = arguments()
    spliterator = FASTAReader(args.input, preprocess_header=lambda h: h.lstrip(">").split(args.split)[0])
    if not args.allow_stop_codons:
        spliterator = itertools.filterfalse(lambda tup: FASTA_STOP_CODON in tup[1], spliterator)

    if args.verbose:
        spliterator = tqdm(spliterator, desc="filter-fasta", ascii=True)

    total = 0
    satisfied = 0
    for header, sequence in spliterator:
        total += 1 
        if args.condition(header, sequence):
            record = f">{header}\n{sequence}\n"
            args.output.write(record)
            satisfied += 1

    args.output.close()
    
    print(f"Extracted {satisfied}/{total} sequences", file=sys.stderr)
