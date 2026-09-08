import csv
import re
import os
import glob
import math
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator


# ============================================================
# USER SETTINGS
# ============================================================

# Folder containing all RDP CSV files
CSV_FOLDER = None

# ------------------------------------------------------------
# Fixed physical size of EACH individual event-map figure
# ------------------------------------------------------------

INDIVIDUAL_WIDTH = 11
INDIVIDUAL_HEIGHT = 6

# ------------------------------------------------------------
# Fixed physical size of EACH panel in combined figure
# ------------------------------------------------------------

PANEL_WIDTH = 11
PANEL_HEIGHT = 3.0

# ------------------------------------------------------------
# Output resolution
# ------------------------------------------------------------

DPI = 1200

# ------------------------------------------------------------
# Approximate number of histogram bins
# ------------------------------------------------------------

TARGET_BINS = 25


# ============================================================
# COMMAND-LINE ARGUMENTS
# ============================================================

parser = argparse.ArgumentParser(
    description="Generate recombination breakpoint/event maps from RDP4 CSV files."
)
parser.add_argument("--input-dir", required=True, help="Directory containing RDP4 CSV files.")
parser.add_argument(
    "--region", required=True, choices=["upstream", "downstream"],
    help="Region to analyse."
)
args = parser.parse_args()

CSV_FOLDER = os.path.abspath(args.input_dir)
region = args.region.lower()


# ============================================================
# FIND ALL CSV FILES
# ============================================================

all_csv_files = glob.glob(
    os.path.join(
        CSV_FOLDER,
        "*.csv"
    )
)


selected_files = []


for filepath in all_csv_files:

    filename = os.path.basename(
        filepath
    ).lower()

    if region == "upstream":

        if (
            "upstream" in filename
            or "_up" in filename
        ):

            selected_files.append(
                filepath
            )

    elif region == "downstream":

        if (
            "downstream" in filename
            or "_down" in filename
        ):

            selected_files.append(
                filepath
            )


# ============================================================
# CHECK FILES
# ============================================================

if len(selected_files) == 0:

    raise FileNotFoundError(
        f"\nNo {region} CSV files were found in:\n"
        f"{CSV_FOLDER}\n\n"
        "Make sure the filenames contain "
        "'upstream' or 'downstream'."
    )


selected_files = sorted(
    selected_files
)


print("\n==============================================")
print(f"{region.upper()} CSV FILES FOUND")
print("==============================================")


for filepath in selected_files:

    print(
        os.path.basename(filepath)
    )


print(
    f"\nTotal CSV files: {len(selected_files)}"
)


# ============================================================
# FUNCTION: GET GENE NAME
# ============================================================

def get_gene_name(filepath):

    filename = os.path.splitext(
        os.path.basename(filepath)
    )[0]

    name = filename.lower()


    # --------------------------------------------------------
    # Remove upstream/downstream terms
    # --------------------------------------------------------

    name = re.sub(
        r"[_\-]?(upstream|downstream)",
        "",
        name,
        flags=re.IGNORECASE
    )


    name = re.sub(
        r"[_\-]?(up|down)$",
        "",
        name,
        flags=re.IGNORECASE
    )


    # --------------------------------------------------------
    # Recognize common beta-lactamase families
    # --------------------------------------------------------

    patterns = [

        (r"oxa[-_]?48", "blaOXA-48"),

        (r"imp", "blaIMP"),

        (r"kpc", "blaKPC"),

        (r"ndm", "blaNDM"),

        (r"shv", "blaSHV"),

        (r"tem", "blaTEM"),

        (r"oxa", "blaOXA")
    ]


    for pattern, gene in patterns:

        if re.search(
            pattern,
            name,
            flags=re.IGNORECASE
        ):

            return gene


    # --------------------------------------------------------
    # Fallback
    # --------------------------------------------------------

    name = re.sub(
        r"[^A-Za-z0-9\-]+",
        "",
        name
    )

    return name


# ============================================================
# FUNCTION: READ RDP CSV
# ============================================================

def read_rdp_csv(filepath):

    with open(
        filepath,
        encoding="utf-8-sig",
        newline=""
    ) as f:

        rows = list(
            csv.reader(f)
        )


    events = []


    # --------------------------------------------------------
    # Extract first integer from a field
    # --------------------------------------------------------

    def extract_number(value):

        match = re.search(
            r"\d+",
            str(value)
        )

        if match:

            return int(
                match.group()
            )

        return np.nan


    # --------------------------------------------------------
    # RDP event table
    #
    # Column 0 = event number
    # Column 2 = beginning coordinate
    # Column 3 = ending coordinate
    # --------------------------------------------------------

    for row in rows[3:]:

        if len(row) <= 10:
            continue


        if not row[0].strip().isdigit():
            continue


        if not row[2].strip():
            continue


        event_number = int(
            row[0].strip()
        )


        begin = extract_number(
            row[2]
        )


        end = extract_number(
            row[3]
        )


        if (
            np.isnan(begin)
            or np.isnan(end)
        ):

            continue


        begin = int(begin)
        end = int(end)


        # Make sure beginning < ending
        if end < begin:

            begin, end = end, begin


        events.append(
            {
                "event": event_number,
                "begin": begin,
                "end": end
            }
        )


    # --------------------------------------------------------
    # Keep only first occurrence of each event
    # --------------------------------------------------------

    unique_events = []

    seen = set()


    for event in events:

        if event["event"] not in seen:

            unique_events.append(
                event
            )

            seen.add(
                event["event"]
            )


    return sorted(
        unique_events,
        key=lambda x: (
            x["begin"],
            x["end"]
        )
    )


# ============================================================
# READ ALL DATASETS
# ============================================================

datasets = []


for filepath in selected_files:

    try:

        events = read_rdp_csv(
            filepath
        )

    except Exception as error:

        print(
            f"\nERROR reading "
            f"{os.path.basename(filepath)}"
        )

        print(error)

        continue


    if len(events) == 0:

        print(
            f"\nWARNING: No RDP events found in "
            f"{os.path.basename(filepath)}"
        )

        continue


    gene = get_gene_name(
        filepath
    )


    # ========================================================
    # ALL BREAKPOINT COORDINATES
    # ========================================================

    breakpoints = []


    for event in events:

        breakpoints.append(
            event["begin"]
        )

        breakpoints.append(
            event["end"]
        )


    breakpoints = np.array(
        breakpoints,
        dtype=float
    )


    # ========================================================
    # X-AXIS RANGE FROM ACTUAL CSV DATA
    # ========================================================

    xmin_data = float(
        np.min(
            breakpoints
        )
    )


    xmax_data = float(
        np.max(
            breakpoints
        )
    )


    data_range = (
        xmax_data -
        xmin_data
    )


    if data_range <= 0:

        data_range = 1


    # --------------------------------------------------------
    # Small automatic visual margin
    #
    # This is NOT a fixed alignment length.
    # --------------------------------------------------------

    x_margin = (
        data_range * 0.025
    )


    xmin_plot = (
        xmin_data -
        x_margin
    )


    xmax_plot = (
        xmax_data +
        x_margin
    )


    # ========================================================
    # HISTOGRAM BINNING
    # ========================================================

    bin_count = min(
        TARGET_BINS,
        max(
            8,
            int(
                np.ceil(
                    data_range / 50
                )
            )
        )
    )


    if xmax_data == xmin_data:

        bins = np.array(
            [
                xmin_data - 1,
                xmax_data + 1
            ]
        )

    else:

        bins = np.linspace(
            xmin_data,
            xmax_data,
            bin_count + 1
        )


    counts, edges = np.histogram(
        breakpoints,
        bins=bins
    )


    centers = (
        edges[:-1] +
        edges[1:]
    ) / 2


    widths = np.diff(
        edges
    )


    # ========================================================
    # STORE DATASET
    # ========================================================

    datasets.append(
        {
            "file": filepath,

            "gene": gene,

            "events": events,

            "breakpoints": breakpoints,

            "xmin_data": xmin_data,

            "xmax_data": xmax_data,

            "xmin_plot": xmin_plot,

            "xmax_plot": xmax_plot,

            "counts": counts,

            "edges": edges,

            "centers": centers,

            "widths": widths
        }
    )


# ============================================================
# CHECK DATA
# ============================================================

if len(datasets) == 0:

    raise ValueError(
        "No valid RDP datasets were found."
    )


# ============================================================
# SORT GENES
# ============================================================

preferred_order = [
    "blaIMP",
    "blaKPC",
    "blaTEM",
    "blaNDM",
    "blaSHV",
    "blaOXA-48"
]


def sort_key(dataset):

    gene = dataset["gene"]


    if gene in preferred_order:

        return (
            0,
            preferred_order.index(gene)
        )


    return (
        1,
        gene
    )


datasets = sorted(
    datasets,
    key=sort_key
)


# ============================================================
# PANEL NUMBER
# ============================================================

def panel_number(index):

    return str(index + 1)


# ============================================================
# ============================================================
#
# PART 1
# INDIVIDUAL RECOMBINATION-EVENT FIGURES
#
# ============================================================
# ============================================================

print(
    "\nCreating individual recombination-event maps..."
)


individual_output_dir = os.path.join(
    CSV_FOLDER,
    f"{region}_individual_event_maps"
)


os.makedirs(
    individual_output_dir,
    exist_ok=True
)


# ============================================================
# EVENT COLOR PALETTE
# ============================================================

event_cmap = plt.get_cmap(
    "tab20"
)


# ============================================================
# CREATE ONE EVENT MAP FOR EACH GENE
# ============================================================

for index, dataset in enumerate(
    datasets
):

    gene = dataset["gene"]

    events = dataset["events"]

    n_events = len(
        events
    )


    # ========================================================
    # FIXED PHYSICAL SIZE
    # ========================================================

    fig, ax = plt.subplots(

        figsize=(
            INDIVIDUAL_WIDTH,
            INDIVIDUAL_HEIGHT
        ),

        dpi=DPI
    )


    # ========================================================
    # EVENT BAR HEIGHT
    # ========================================================

    bar_height = 0.62


    # ========================================================
    # DRAW RECOMBINATION EVENTS
    # ========================================================

    for i, event in enumerate(
        events
    ):

        y = (
            n_events -
            i
        )


        color = event_cmap(
            i % 20
        )


        # ----------------------------------------------------
        # Recombination interval
        # ----------------------------------------------------

        ax.add_patch(

            Rectangle(

                (
                    event["begin"],

                    y -
                    bar_height / 2
                ),

                event["end"] -
                event["begin"],

                bar_height,

                facecolor=color,

                edgecolor="black",

                linewidth=0.8,

                alpha=0.90
            )
        )


        # ----------------------------------------------------
        # Breakpoint markers
        # ----------------------------------------------------

        ax.scatter(

            [
                event["begin"],
                event["end"]
            ],

            [
                y,
                y
            ],

            s=28,

            facecolors="white",

            edgecolors="black",

            linewidths=0.9,

            zorder=5
        )


        # ----------------------------------------------------
        # Event label
        # ----------------------------------------------------

        label_offset = max(

            (
                dataset["xmax_data"] -
                dataset["xmin_data"]
            ) * 0.015,

            5
        )


        ax.text(

            event["begin"] -
            label_offset,

            y,

            f"E{event['event']}",

            ha="right",

            va="center",

            fontsize=10,

            fontweight="normal"
        )


    # ========================================================
    # DATA-DRIVEN X AXIS
    # ========================================================

    ax.set_xlim(

        dataset["xmin_plot"],

        dataset["xmax_plot"]
    )


    # ========================================================
    # Y AXIS
    # ========================================================

    ax.set_ylim(

        0.3,

        n_events + 0.7
    )


    ax.set_yticks([])


    # ========================================================
    # AXIS LABELS
    # ========================================================

    ax.set_xlabel(

        "Alignment coordinate",

        fontsize=12
    )


    ax.set_ylabel(

        "RDP4 recombination events",

        fontsize=12
    )


    # ========================================================
    # AUTOMATIC X TICKS
    # ========================================================

    ax.xaxis.set_major_locator(

        MaxNLocator(

            nbins=7,

            integer=True
        )
    )


    # ========================================================
    # GRID
    # ========================================================

    ax.grid(

        axis="x",

        alpha=0.18,

        linewidth=0.8
    )


    # ========================================================
    # SPINES
    # ========================================================

    ax.spines[
        "top"
    ].set_visible(False)


    ax.spines[
        "right"
    ].set_visible(False)


    ax.spines[
        "left"
    ].set_visible(False)


    # ========================================================
    # GENE LABEL
    #
    # MOVED LEFT + UP
    #
    # 1.blaIMP
    #
    # Entire label is normal upright text.
    # ========================================================

    number = panel_number(
        index
    )


    ax.text(

        -0.015,

        1.12,

        f"{number}.{gene}",

        transform=ax.transAxes,

        ha="left",

        va="bottom",

        fontsize=18,

        fontweight="normal",

        fontstyle="normal",

        clip_on=False
    )


    # ========================================================
    # EXTRA TOP SPACE
    # ========================================================

    fig.subplots_adjust(

        left=0.10,

        right=0.98,

        top=0.88,

        bottom=0.12
    )


    # ========================================================
    # SAFE GENE NAME
    # ========================================================

    safe_gene = re.sub(

        r"[^A-Za-z0-9\-]+",

        "_",

        gene
    )


    # ========================================================
    # OUTPUT FILES
    # ========================================================

    individual_png = os.path.join(

        individual_output_dir,

        f"{safe_gene}_recombination_events.png"
    )


    individual_jpg = os.path.join(

        individual_output_dir,

        f"{safe_gene}_recombination_events.jpg"
    )


    individual_svg = os.path.join(

        individual_output_dir,

        f"{safe_gene}_recombination_events.svg"
    )


    # ========================================================
    # SAVE
    # ========================================================

    fig.savefig(

        individual_png,

        dpi=DPI,

        bbox_inches="tight"
    )


    fig.savefig(

        individual_jpg,

        dpi=DPI,

        bbox_inches="tight"
    )


    fig.savefig(

        individual_svg,

        bbox_inches="tight"
    )


    plt.close(
        fig
    )


    print(
        f"  Created individual map: {gene}"
    )


# ============================================================
# ============================================================
#
# PART 2
# COMBINED BREAKPOINT-FREQUENCY PANEL
#
# ============================================================
# ============================================================

print(
    "\nCreating combined breakpoint-frequency panel..."
)


n_panels = len(
    datasets
)


# One gene per row
n_cols = 1

n_rows = n_panels


# ============================================================
# FIXED PHYSICAL SIZE
#
# Every panel has identical physical width and height.
# ============================================================

fig_width = PANEL_WIDTH

fig_height = (
    PANEL_HEIGHT *
    n_rows
)


fig, axes = plt.subplots(

    n_rows,

    n_cols,

    figsize=(
        fig_width,
        fig_height
    ),

    dpi=DPI,

    squeeze=False
)


axes = axes.flatten()


# ============================================================
# DRAW EACH BREAKPOINT PANEL
# ============================================================

for index, dataset in enumerate(
    datasets
):

    ax = axes[index]


    # ========================================================
    # BREAKPOINT HISTOGRAM
    # ========================================================

    ax.bar(

        dataset["centers"],

        dataset["counts"],

        width=dataset["widths"],

        align="center",

        edgecolor="black",

        linewidth=0.5,

        alpha=0.90
    )


    # ========================================================
    # X RANGE FROM CSV
    #
    # NOT FIXED TO 1500 / 2000 / ETC.
    # ========================================================

    ax.set_xlim(

        dataset["xmin_plot"],

        dataset["xmax_plot"]
    )


    # ========================================================
    # Y RANGE FROM CSV
    # ========================================================

    max_count = max(

        dataset["counts"]
    )


    if max_count < 1:

        max_count = 1


    ax.set_ylim(

        0,

        max_count * 1.12
    )


    # ========================================================
    # AUTOMATIC X TICKS
    # ========================================================

    ax.xaxis.set_major_locator(

        MaxNLocator(

            nbins=8,

            integer=True
        )
    )


    # ========================================================
    # AUTOMATIC Y TICKS
    # ========================================================

    ax.yaxis.set_major_locator(

        MaxNLocator(

            nbins=5,

            integer=True
        )
    )


    # ========================================================
    # REMOVE INDIVIDUAL AXIS LABELS
    # ========================================================

    ax.set_xlabel("")
    ax.set_ylabel("")


    # ========================================================
    # GENE LABEL
    #
    # IMPORTANT:
    #
    # x = -0.015
    #     slightly toward LEFT
    #
    # y = 1.12
    #     clearly ABOVE DATA
    #
    # ALL TEXT NORMAL / UPRIGHT
    # ========================================================

    number = panel_number(
        index
    )

    gene = dataset["gene"]


    ax.text(

        -0.015,

        1.12,

        f"{number}.{gene}",

        transform=ax.transAxes,

        ha="left",

        va="bottom",

        fontsize=18,

        fontweight="normal",

        fontstyle="normal",

        clip_on=False
    )


    # ========================================================
    # SPINES
    # ========================================================

    ax.spines[
        "top"
    ].set_visible(False)


    ax.spines[
        "right"
    ].set_visible(False)


    # ========================================================
    # TICK STYLE
    # ========================================================

    ax.tick_params(

        axis="both",

        labelsize=10,

        direction="out"
    )


# ============================================================
# COMMON X AXIS LABEL
# ============================================================

fig.supxlabel(

    "Alignment coordinate",

    fontsize=19,

    fontweight="normal",

    y=0.012
)


# ============================================================
# COMMON Y AXIS LABEL
# ============================================================

fig.supylabel(

    "Breakpoint count",

    fontsize=19,

    fontweight="normal",

    x=0.015
)


# ============================================================
# PANEL SPACING
#
# Extra vertical space keeps the gene labels away from
# the previous panel.
# ============================================================

fig.subplots_adjust(

    left=0.09,

    right=0.985,

    top=0.96,

    bottom=0.07,

    hspace=0.90
)


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

panel_output_dir = os.path.join(

    CSV_FOLDER,

    f"{region}_RDP_breakpoint_panel"
)


os.makedirs(

    panel_output_dir,

    exist_ok=True
)


# ============================================================
# OUTPUT FILES
# ============================================================

panel_png = os.path.join(

    panel_output_dir,

    f"{region}_breakpoint_frequency_panel.png"
)


panel_jpg = os.path.join(

    panel_output_dir,

    f"{region}_breakpoint_frequency_panel.jpg"
)


panel_svg = os.path.join(

    panel_output_dir,

    f"{region}_breakpoint_frequency_panel.svg"
)


# ============================================================
# SAVE COMBINED PANEL
# ============================================================

fig.savefig(

    panel_png,

    dpi=DPI,

    bbox_inches="tight"
)


fig.savefig(

    panel_jpg,

    dpi=DPI,

    bbox_inches="tight"
)


fig.savefig(

    panel_svg,

    bbox_inches="tight"
)


plt.show()


# ============================================================
# FINAL SUMMARY
# ============================================================

print("\n")
print("======================================================")
print("RDP RECOMBINATION FIGURES COMPLETE")
print("======================================================")


print(
    f"Region: {region.upper()}"
)


print(
    f"Genes plotted: {len(datasets)}"
)


print(
    "\nCoordinate ranges calculated directly from CSV:"
)


for index, dataset in enumerate(
    datasets
):

    print(

        f"{panel_number(index)}. "
        f"{dataset['gene']} : "
        f"{dataset['xmin_data']:.0f} - "
        f"{dataset['xmax_data']:.0f} "
        f"| Events = "
        f"{len(dataset['events'])}"
    )


print(
    "\nIndividual event maps:"
)

print(
    individual_output_dir
)


print(
    "\nCombined breakpoint panel:"
)

print(
    panel_output_dir
)


print(
    "\nCombined panel files:"
)

print(
    panel_png
)

print(
    panel_jpg
)

print(
    panel_svg
)


print(
    "\n======================================================"
)