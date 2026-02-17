#!/usr/bin/env python3

import argparse

def process_file(input_file, output_file):
    with open(input_file, 'r') as src_file:
        lines = src_file.readlines()

    with open(output_file, 'w') as dest_file:
        for line in lines:
            # 检查行是否包含"ORGANISM"，如果是，则在"ORGANISM"之前插入一个换行符
            if 'ORGANISM' in line:
                parts = line.split('ORGANISM', 1)  # 仅分割第一次出现的"ORGANISM"
                line = parts[0] + '\nORGANISM' + parts[1]
            # 直接写入包含"CRISPR"或"ORGANISM"的行，或者不包含排除关键词的行
            if 'CRISPR' in line or 'ORGANISM' in line or (line.strip() != '' and 'Time' not in line and '-' not in line and 'Repeats' not in line and 'POSITION' not in line):
                dest_file.write(line)

def main():
    parser = argparse.ArgumentParser(description='Process a single CRISPR report file.')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')

    args = parser.parse_args()

    process_file(args.input, args.output)
    print("File has been processed successfully.")

if __name__ == "__main__":
    main()
