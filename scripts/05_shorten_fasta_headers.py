#!/usr/bin/env python3
from pathlib import Path
import sys

def shorten_fasta_headers(input_fasta: str, output_fasta: str) -> None:
    input_path = Path(input_fasta)
    output_path = Path(output_fasta)

    with input_path.open("r", encoding="utf-8") as fin, \
         output_path.open("w", encoding="utf-8") as fout:

        for line in fin:
            line = line.rstrip("\n")

            if line.startswith(">"):
                header = line[1:]  # remove >
                parts = header.split("|")

                # keep only the first two parts if available
                if len(parts) >= 2:
                    new_header = "|".join(parts[:2])
                else:
                    new_header = header

                fout.write(f">{new_header}\n")
            else:
                fout.write(line + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python shorten_fasta.py input.fasta output.fasta")
        sys.exit(1)

    shorten_fasta_headers(sys.argv[1], sys.argv[2])