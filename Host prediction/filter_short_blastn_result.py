#!/usr/bin/env python3

import argparse
import csv

def main():
    parser = argparse.ArgumentParser(description="Filter BLAST TSV results with strict near-perfect match criteria.")
    parser.add_argument("-i", "--input", required=True, help="Input BLAST TSV file")
    parser.add_argument("-o", "--output", required=True, help="Output filtered TSV file (only written if valid rows exist)")
    args = parser.parse_args()

    valid_rows = []

    with open(args.input, 'r') as infile:
        reader = csv.reader(infile, delimiter='\t')
        for row in reader:
            # 跳过行数不足的异常行
            if len(row) < 13:
                continue

            try:
                col_4 = int(row[3])   # alignment length
                col_5 = int(row[4])   # mismatches
                col_6 = int(row[5])   # gaps
                col_13 = int(row[12]) # query sequence length
            except ValueError:
                continue  # 跳过无法转为整数的行

            # 严格过滤条件
            condition1 = (col_13 - col_4 == 1) and (col_5 + col_6 == 0)
            condition2 = (col_13 - col_4 == 0) and (col_5 + col_6 <= 1)

            if condition1 or condition2:
                valid_rows.append(row)

    # 仅当有符合条件的行时才写入输出文件
    if valid_rows:
        with open(args.output, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter='\t')
            writer.writerows(valid_rows)
        print(f"Filtered results saved to {args.output}")
    else:
        print("No rows met the filtering criteria. Output file not created.")

if __name__ == "__main__":
    main()
