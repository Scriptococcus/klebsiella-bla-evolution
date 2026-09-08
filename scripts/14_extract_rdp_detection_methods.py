import os
import csv
import re
import pandas as pd
import argparse


# ============================================================
# FOLDER CONTAINING RDP4 CSV FILES
# ============================================================

CSV_FOLDER = None


# ============================================================
# DETECTION METHODS TO EXTRACT
# ============================================================

METHODS = [
    "RDP",
    "GENECONV",
    "Bootscan",
    "Maxchi",
    "Chimaera",
    "SiSscan",
    "PhylPro",
    "LARD",
    "3Seq"
]


# ============================================================
# CLEAN TEXT
# ============================================================

def clean(value):
    if value is None:
        return ""

    return str(value).strip().replace("\xa0", " ")


# ============================================================
# COMMAND-LINE ARGUMENTS / INPUT FILE
# ============================================================

parser = argparse.ArgumentParser(
    description="Extract recombination-detection method results from an RDP4 CSV export."
)
parser.add_argument("--input", required=True, help="RDP4 CSV file to process.")
parser.add_argument(
    "--output", default=None,
    help="Output Excel path. Defaults to <input_stem>_detection_methods.xlsx."
)
args = parser.parse_args()

input_file = os.path.abspath(args.input)
selected_file = os.path.basename(input_file)
CSV_FOLDER = os.path.dirname(input_file)

print("\nSelected:")
print(selected_file)


# ============================================================
# READ THE RDP CSV
#
# IMPORTANT:
# RDP4 exports have several introductory rows before
# the actual column header.
#
# We therefore use csv.reader instead of pandas.read_csv().
# ============================================================

print("\nReading RDP4 CSV...")


encodings = [
    "utf-8-sig",
    "utf-8",
    "cp1252",
    "latin1"
]

rows = None

for encoding in encodings:

    try:

        with open(
            input_file,
            "r",
            encoding=encoding,
            newline=""
        ) as f:

            reader = csv.reader(f)

            rows = list(reader)

        print(f"Encoding: {encoding}")
        break

    except UnicodeDecodeError:
        continue


if rows is None:

    print("\nERROR: Could not read the CSV.")

    input("\nPress Enter to exit...")
    raise SystemExit


# ============================================================
# FIND THE REAL HEADER ROW
# ============================================================

header_index = None

for i, row in enumerate(rows):

    cleaned_row = [
        clean(x).lower()
        for x in row
    ]

    if (
        any(
            "recombination event number" in x
            for x in cleaned_row
        )
        and
        any(
            x == "rdp"
            for x in cleaned_row
        )
        and
        any(
            x == "genecconv" or x == "geneconv"
            for x in cleaned_row
        )
    ):

        header_index = i
        break


# ============================================================
# FALLBACK HEADER DETECTION
# ============================================================

if header_index is None:

    for i, row in enumerate(rows):

        cleaned_row = [
            clean(x).lower()
            for x in row
        ]

        if (
            "recombination event number" in cleaned_row
            and "rdp" in cleaned_row
        ):

            header_index = i
            break


if header_index is None:

    print("\nERROR:")
    print("Could not locate the RDP4 table header.")

    input("\nPress Enter to exit...")
    raise SystemExit


print(
    f"\nRDP4 table header found on CSV row "
    f"{header_index + 1}"
)


# ============================================================
# GET HEADER
# ============================================================

header = [
    clean(x)
    for x in rows[header_index]
]


print("\nDetected RDP4 columns:")

for i, column in enumerate(header):

    print(
        f"  {i}: {column}"
    )


# ============================================================
# FIND IMPORTANT COLUMN POSITIONS
# ============================================================

def find_column(name):

    target = name.lower().replace(" ", "")

    for i, column in enumerate(header):

        current = (
            clean(column)
            .lower()
            .replace(" ", "")
        )

        if current == target:
            return i

    return None


event_col = find_column(
    "Recombination Event Number"
)


if event_col is None:

    print(
        "\nERROR: Recombination Event Number column "
        "was not found."
    )

    input("\nPress Enter to exit...")
    raise SystemExit


method_columns = {}

for method in METHODS:

    col = find_column(method)

    if col is not None:
        method_columns[method] = col


# ============================================================
# SHOW METHOD COLUMNS
# ============================================================

print("\nDetection methods:")

for method in METHODS:

    if method in method_columns:

        print(
            f"  {method}: column "
            f"{method_columns[method]}"
        )

    else:

        print(
            f"  {method}: NOT FOUND"
        )


# ============================================================
# EXTRACT EVENTS
#
# VERY IMPORTANT:
#
# In the RDP4 CSV, each recombination event can occupy
# several rows.
#
# Example:
#
# Row 1:
# event = 1
# detection methods = values
#
# Row 2:
# event = 1
# detection methods = empty
#
# Row 3:
# event = 1
# detection methods = empty
#
# ...
#
# Therefore:
#
# We ONLY extract a row when the detection-method section
# contains actual values.
# ============================================================

records = []

current_event = None


for row in rows[header_index + 1:]:

    if not row:
        continue


    # --------------------------------------------------------
    # Make sure the row has enough fields
    # --------------------------------------------------------

    # RDP rows can have different lengths.
    # We simply pad short rows with empty values.

    if len(row) < len(header):

        row = row + [""] * (
            len(header) - len(row)
        )


    # --------------------------------------------------------
    # Event number
    # --------------------------------------------------------

    event_value = clean(
        row[event_col]
        if event_col < len(row)
        else ""
    )


    # --------------------------------------------------------
    # Ignore completely empty rows
    # --------------------------------------------------------

    if event_value == "":
        continue


    # --------------------------------------------------------
    # Extract numeric event number
    # --------------------------------------------------------

    match = re.search(
        r"\d+",
        event_value
    )

    if not match:
        continue


    event_number = int(
        match.group()
    )


    # --------------------------------------------------------
    # Check whether this is the MAIN event row
    #
    # Main event rows contain detection-method results.
    # Additional sequence rows do not.
    # --------------------------------------------------------

    method_values = []

    for method in METHODS:

        if method in method_columns:

            col = method_columns[method]

            if col < len(row):

                value = clean(
                    row[col]
                )

            else:

                value = ""

            method_values.append(value)


    # If all detection method cells are empty,
    # this is NOT the main event row.
    if not any(
        value != ""
        for value in method_values
    ):

        continue


    # --------------------------------------------------------
    # Create one row
    # --------------------------------------------------------

    record = {
        "Recombination Event Number": event_number
    }


    for method in METHODS:

        if method in method_columns:

            col = method_columns[method]

            if col < len(row):

                record[method] = clean(
                    row[col]
                )

            else:

                record[method] = ""

        else:

            record[method] = ""


    records.append(record)


# ============================================================
# CREATE DATAFRAME
# ============================================================

result = pd.DataFrame(
    records,
    columns=[
        "Recombination Event Number"
    ] + METHODS
)


# ============================================================
# REMOVE DUPLICATE EVENT NUMBERS
#
# Normally there should already be one main row per event.
# This is an additional safety step.
# ============================================================

result = result.drop_duplicates(
    subset=["Recombination Event Number"],
    keep="first"
)


# ============================================================
# SORT EVENTS
# ============================================================

result = result.sort_values(
    "Recombination Event Number"
).reset_index(
    drop=True
)


# ============================================================
# OUTPUT FILE
# ============================================================

base_name = os.path.splitext(
    selected_file
)[0]


output_file = os.path.abspath(
    args.output if args.output else os.path.join(
        CSV_FOLDER,
        base_name + "_detection_methods.xlsx"
    )
)


# ============================================================
# WRITE EXCEL
# ============================================================

result.to_excel(
    output_file,
    index=False,
    sheet_name="Detection Methods"
)


# ============================================================
# FINISHED
# ============================================================

print("\n==============================================")
print("DONE")
print("==============================================")

print(
    f"\nRecombination events extracted: "
    f"{len(result)}"
)

print("\nExcel created:")

print(output_file)

print("\nPreview:\n")

print(
    result.head(10).to_string(
        index=False
    )
)

input("\nPress Enter to exit...")