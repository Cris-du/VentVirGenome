#!/usr/bin/env python3

import argparse

def filter_blast_results(input_file, output_file):
    filtered_lines = []  # 用于暂存符合条件的行
    with open(input_file, 'r') as infile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) < 14:
                continue  # 确保行有足够的列数
            
            # 解析必要的字段
            query_name, subject_name, identity, alignment_length, evalue, bit_score = (
                parts[0], parts[1], float(parts[2]), int(parts[3]), float(parts[10]), float(parts[11])
            )
            
            # 应用筛选条件
            if (query_name != subject_name and 
                alignment_length >= 2500 and 
                identity >= 70.0 and 
                evalue <= 0.001 and 
                bit_score >= 50):
                filtered_lines.append(line)  # 只有符合条件的行才被添加到列表中
    
    # 只有当有符合条件的行时，才创建输出文件并写入数据
    if filtered_lines:
        with open(output_file, 'w') as outfile:
            for line in filtered_lines:
                outfile.write(line)

def main():
    parser = argparse.ArgumentParser(description='Filter BLAST .tsv results based on specific criteria.')
    parser.add_argument('-i', '--input', required=True, help='Input file path.')
    parser.add_argument('-o', '--output', required=True, help='Output file path.')

    args = parser.parse_args()

    filter_blast_results(args.input, args.output)

if __name__ == '__main__':
    main()
