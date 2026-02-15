#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description='Split sequences based on provirus keyword and TSV annotation.')
    parser.add_argument('-ifg', '--input_raw_fna', required=True, help='Input raw sequences (FASTA).')
    parser.add_argument('-ifc', '--input_provirus_fna', required=True, help='Input provirus candidate sequences (FASTA).')
    parser.add_argument('-it', '--input_tsv', required=True, help='TSV file: col1=seqID, col3=Yes/No.')
    parser.add_argument('-o', '--output_non_provirus', required=True, help='Output: raw seqs WITHOUT "provirus" in name.')
    parser.add_argument('-og', '--output_provirus_no', required=True, help='Output: raw seqs WITH "provirus" in name AND TSV says "No".')
    parser.add_argument('-ogc', '--output_provirus_confirmed', required=True, help='Output: provirus file seqs WITH "provirus" in name.')

    args = parser.parse_args()

    # Step 1: Load TSV annotation into dict: {seq_id: status}
    tsv_status = {}
    with open(args.input_tsv, 'r') as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                seq_id = parts[0]
                status = parts[2]
                if status in ('Yes', 'No'):
                    tsv_status[seq_id] = status

    # Step 2: Process -ifg (raw sequences)
    with open(args.output_non_provirus, 'w') as out_o, \
         open(args.output_provirus_no, 'w') as out_og:

        for record in SeqIO.parse(args.input_raw_fna, 'fasta'):
            seq_id = record.id
            has_provirus = 'provirus' in seq_id.lower()

            if not has_provirus:
                # Rule 1: output to -o
                SeqIO.write(record, out_o, 'fasta')
            else:
                # Rule 2: only output to -og if TSV says "No"
                if tsv_status.get(seq_id) == 'No':
                    SeqIO.write(record, out_og, 'fasta')
                # If not in TSV or status is "Yes", discard

    # Step 3: Process -ifc (provirus sequences)
    with open(args.output_provirus_confirmed, 'w') as out_ogc:
        for record in SeqIO.parse(args.input_provirus_fna, 'fasta'):
            seq_id = record.id
            if 'provirus' in seq_id.lower():
                # Rule 3: output to -ogc if name contains "provirus"
                SeqIO.write(record, out_ogc, 'fasta')
            # else: skip

if __name__ == '__main__':
    main()
