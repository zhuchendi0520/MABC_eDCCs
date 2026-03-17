#!/usr/bin/env python3

import sys
from collections import defaultdict
import os

def load_cc_positions(cc_file):
    cc_dict = defaultdict(set)
    try:
        with open(cc_file, 'r') as f:
            next(f)  
            for line in f:
                line = line.strip().split("\t")
                if len(line) != 3:
                    continue
                try:
                    subtype, position, ref_base = line[0], int(line[1]), line[2]
                    cc_dict[subtype].add((position, ref_base))
                except ValueError:
                    print(f"Warning: Skipping malformed line in {cc_file}: {line}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: File '{cc_file}' not found.", file=sys.stderr)
        sys.exit(1)
    return cc_dict

def load_snp_file(snp_file):
    snp_set = set()
    try:
        with open(snp_file, 'r') as f:
            for line in f:
                cols = line.strip().split()
                if len(cols) < 3:
                    continue
                try:
                    position, var_base = int(cols[0]), cols[2]
                    snp_set.add((position, var_base))
                except ValueError:
                    print(f"Warning: Skipping malformed line in {snp_file}: {cols}", file=sys.stderr)
    except FileNotFoundError:
        print(f"Error: File '{snp_file}' not found.", file=sys.stderr)
        sys.exit(1)
    return snp_set

def calculate_similarity(snp_set, cc_dict):
    results = []
    for cc_name, cc_positions in cc_dict.items():
        if not cc_positions:
            continue
        common_snps = snp_set & cc_positions  
        match_count = len(common_snps)
        total_count = len(cc_positions)
        percentage = (match_count / total_count * 100) if total_count else 0 
        results.append((cc_name, percentage, match_count, total_count))

    if not results:
        return None 

    # 按匹配度排序（降序）
    results.sort(key=lambda x: x[1], reverse=True)
    
    # 获取最高匹配度
    max_match = results[0][1]  # 最高匹配百分比
    best_matches = [r for r in results if r[1] == max_match]  # 找到所有匹配度最高的 CC
    return best_matches

def main():
    if len(sys.argv) < 3:
        print("Usage: script.py <snp_file> <cc_position_file>", file=sys.stderr)
        sys.exit(1)

    snp_file = sys.argv[2]  
    cc_file = sys.argv[1]   
    snp_filename = os.path.splitext(os.path.basename(snp_file))[0]  

    snp_set = load_snp_file(snp_file)
    cc_dict = load_cc_positions(cc_file)

    best_matches = calculate_similarity(snp_set, cc_dict)


    if best_matches:
        for cc_name, percentage, match_count, total_count in best_matches:
            if percentage >= 95:
                category = cc_name
            else:
                category = "Non_DCC"
                
            print(f"{snp_filename}\t{cc_name}\t{percentage:.2f}%\t({match_count}/{total_count})\t{category}")
    else:
        print(f"{snp_filename}\tNo matching CC found\t0.00%\t(0/0)\teDCC16")

if __name__ == "__main__":
    main()

