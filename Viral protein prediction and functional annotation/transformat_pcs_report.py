#!/usr/bin/env python3

import argparse

def process_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        cluster = ""
        for line in infile:
            if line.startswith('>Cluster'):
                cluster = line.strip().split()[1]  # 获取 Cluster 号码
            else:
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    number = parts[0]  # 获取跟随数字
                    details = parts[1]  # 第二列的信息

                    # 给每个数字前加上 Cluster 并保留数字编号
                    new_cluster_info = f"Cluster_{cluster}"
                    new_number = number

                    # 替换 "aa,>" 和 "..." 为 tab 并分割
                    details = details.replace("aa, >", "\t").replace("...", "\t")

                    # 生成处理后的行并写入
                    new_line = f"{new_cluster_info}\t{new_number}\t{details}\n"
                    outfile.write(new_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process tsv file.")
    parser.add_argument('-i', '--input', required=True, help="Input tsv file")
    parser.add_argument('-o', '--output', required=True, help="Output tsv file")
    
    args = parser.parse_args()
    process_file(args.input, args.output)
