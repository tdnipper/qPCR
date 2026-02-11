#!/usr/bin/env python3

import sys
from qpcr_importer import FileImporter


def usage():
    print("Usage: python filter_results.py <input.csv> <output.csv>")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        usage()
        sys.exit(2)

    infile, outfile = sys.argv[1], sys.argv[2]
    try:
        importer = FileImporter(infile)
        data = importer.import_file()
    except Exception as e:
        print(f"Error importing file: {e}")
        sys.exit(1)

    try:
        data.to_csv(outfile, index=False)
    except Exception as e:
        print(f"Error writing output file '{outfile}': {e}")
        sys.exit(1)
    