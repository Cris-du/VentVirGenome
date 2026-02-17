#!/usr/bin/env python3

import argparse
import os

def parse_crispr_identifier(record_title):
    """
    从记录标题中解析CRISPR标识符。
    """
    crispr_end_index = record_title.find('CRISPR') + record_title[record_title.find('CRISPR'):].find('_')
    crispr_identifier = record_title[:crispr_end_index]
    return crispr_identifier

def process_file(input_file, output_file):
    crispr_counts = {}
    sequences = {}

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                crispr_identifier = parse_crispr_identifier(line[1:].strip())
                crispr_counts[crispr_identifier] = crispr_counts.get(crispr_identifier, 0) + 1
                if crispr_identifier not in sequences:
                    sequences[crispr_identifier] = []
                sequences[crispr_identifier].append(line.strip() + '\n' + next(f).strip())

    sequences_to_write = [seqs for seqs in sequences.values() if len(seqs) >= 3]

    if sequences_to_write:
        with open(output_file, 'w') as output_f:
            for seq_group in sequences_to_write:
                for seq in seq_group:
                    output_f.write(f'{seq}\n')

def main():
    parser = argparse.ArgumentParser(description='Process CRISPR spacer sequences.')
    parser.add_argument('-i', '--input', required=True, help='Input fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output fasta file')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Input file {args.input} does not exist.")
        return

    process_file(args.input, args.output)

if __name__ == "__main__":
    main()
