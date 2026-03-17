import sys

if len(sys.argv) != 3:
    print("usage: python replace_lines.py <input.txt> <output.txt>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for line in infile:
        if not line.strip():
            continue

        parts = line.strip().split()

        if parts[0].lower() == "position":
            outfile.write("\t".join(parts) + "\n")
            continue

        if len(parts) != 3:
            print(f"skip: {line.strip()}")
            continue

        if parts[2] == ".":
            parts[2] = parts[1]

        outfile.write("\t".join(parts) + "\n")

