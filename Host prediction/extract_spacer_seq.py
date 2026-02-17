#!/usr/bin/env python3

import argparse

def extract_spacers_to_fasta(file_path, output_file_path):
    spacers = []
    current_organism = None
    crispr_count = 0  # 将在每次新 organism 时重置

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('ORGANISM:'):
                # 遇到新 organism：提取名称并重置 CRISPR 计数器
                current_organism = line.split('ORGANISM:')[1].strip().replace(' ', '_').replace(',', '')
                crispr_count = 0  # ✅ 修复点：按 organism 重置
            elif current_organism and 'CRISPR' in line and 'Range:' in line:
                crispr_count += 1
            elif current_organism and line and not line.startswith('Bases:') and '[' in line:
                parts = line.split('\t')
                if len(parts) >= 4:
                    spacer_seq = parts[3].strip()
                    if spacer_seq:
                        base_pos = int(parts[0])
                        pos_info = parts[-1].strip('[]').split(',')
                        if len(pos_info) >= 2:
                            try:
                                offset = int(pos_info[0])
                                length = int(pos_info[1])
                                start = base_pos + offset
                                end = start + length - 1
                                header = f">{current_organism}_CRISPR{crispr_count}_{start}_{end}"
                                spacers.append(f"{header}\n{spacer_seq}\n")
                            except ValueError:
                                continue  # 跳过无法解析的位置信息

    # 写入输出
    if spacers:
        with open(output_file_path, 'w') as out:
            out.writelines(spacers)
        print(f"SPACER序列已以FASTA格式保存到 {output_file_path}")
    else:
        print("未找到SPACER序列，不创建文件。")

def main():
    parser = argparse.ArgumentParser(description='Extract CRISPR spacers to FASTA format.')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')

    args = parser.parse_args()

    extract_spacers_to_fasta(args.input, args.output)

if __name__ == "__main__":
    main()
