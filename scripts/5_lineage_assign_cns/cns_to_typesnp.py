#!/usr/bin/env python3
import sys
import os

if len(sys.argv) != 2:
    print("usage: python cns_to_typesnp.py <input.cns>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = input_file + ".typesnp"

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    outfile.write("Position\tRef\tVar\n")

    for line in infile:
        if not line.strip():
            continue

        parts = line.strip().split()

        if parts[0].lower() in ["chrom", "chromosome", "position"]:
            continue

        if len(parts) < 4:
            print(f"skip: {line.strip()}")
            continue

        position = parts[1]
        ref = parts[2]
        var = parts[3]

        if var == ".":
            var = ref

        outfile.write(f"{position}\t{ref}\t{var}\n")