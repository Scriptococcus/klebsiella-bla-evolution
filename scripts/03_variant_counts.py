#!/usr/bin/env python3
from collections import Counter
import argparse
import csv

def count_variants(fasta_file):
    counts = Counter()

    with open(fasta_file, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]  # remove >
                header = header.split()[0]  # keep only first part of header
                variant = header.split("|")[0]  # take SHV-1 from SHV-1|CP166160
                counts[variant] += 1

    return counts

def main():
    parser = argparse.ArgumentParser(
        description="Count number of sequences for each variant in a FASTA file"
    )
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", default="variant_counts.csv",
                        help="Output CSV file (default: variant_counts.csv)")
    args = parser.parse_args()

    counts = count_variants(args.input)

    print("Variant counts:")
    for variant, count in counts.most_common():
        print(f"{variant}\t{count}")

    with open(args.output, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Variant", "Count"])
        for variant, count in counts.most_common():
            writer.writerow([variant, count])

    print(f"\nSaved to {args.output}")

if __name__ == "__main__":
    main()