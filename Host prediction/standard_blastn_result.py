#!/usr/bin/env python3

import argparse

def main():
    parser = argparse.ArgumentParser(description="Extract unique host-virus pairs")
    parser.add_argument("-i", "--input", required=True, help="Input txt file")
    parser.add_argument("-o", "--output", required=True, help="Output txt file")
    args = parser.parse_args()

    unique_pairs = set()

    # 读取输入文件
    with open(args.input, "r") as infile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            host_seq = parts[0]
            virus = parts[1]

            # 提取宿主名称：第三个下划线之前的所有内容
            host_parts = host_seq.split("_")
            if len(host_parts) >= 3:
                host_name = "_".join(host_parts[:3])
            else:
                host_name = host_seq  # 如果不够三个下划线，则原样保留

            unique_pairs.add((host_name, virus))

    # 输出去重后的结果
    with open(args.output, "w") as outfile:
        outfile.write("Host\tVirus\n")
        for host, virus in sorted(unique_pairs):
            outfile.write(f"{host}\t{virus}\n")


if __name__ == "__main__":
    main()
