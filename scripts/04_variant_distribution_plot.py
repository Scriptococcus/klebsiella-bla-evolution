#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def detect_columns(df):
    cols_lower = [c.lower().strip() for c in df.columns]

    variant_candidates = ["variant", "variants", "gene", "genename", "name"]
    count_candidates = ["count", "counts", "accessions",
                        "accession_count", "number", "n"]

    variant_col = None
    count_col = None

    for c in variant_candidates:
        if c in cols_lower:
            variant_col = df.columns[cols_lower.index(c)]
            break

    for c in count_candidates:
        if c in cols_lower:
            count_col = df.columns[cols_lower.index(c)]
            break

    if variant_col is None:
        variant_col = df.columns[0]

    if count_col is None:
        count_col = df.columns[1]

    return variant_col, count_col


def clean_x_label(text):
    """
    IMP-4    -> 4
    KPC-123  -> 123
    NDM-16b  -> 16b
    SHV-27   -> 27
    """
    text = str(text).strip()

    if "-" in text:
        return text.split("-", 1)[1]

    return text


def read_plot_data(csv_file):

    df = pd.read_csv(csv_file)

    variant_col, count_col = detect_columns(df)

    plot_df = df[[variant_col, count_col]].copy()
    plot_df.columns = ["Variant", "Count"]

    plot_df["Variant"] = plot_df["Variant"].astype(str).str.strip()
    plot_df["Count"] = pd.to_numeric(
        plot_df["Count"],
        errors="coerce"
    )

    plot_df = plot_df.dropna(subset=["Count"])
    plot_df = plot_df[plot_df["Count"] > 0]

    # preserve CSV order
    plot_df = plot_df.reset_index(drop=True)

    plot_df["XLabel"] = plot_df["Variant"].apply(clean_x_label)

    return plot_df


def make_panel_plot(inputs, families, output_png):

    data_list = [read_plot_data(f) for f in inputs]

    plt.rcParams.update({
        "font.family": "Arial",
        "font.size": 10,
        "axes.labelsize": 11,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "savefig.dpi": 600
    })

    # A4 manuscript width
    fig_width = 7.0
    fig_height = 10.5

    fig, axes = plt.subplots(
        5,
        1,
        figsize=(fig_width, fig_height),
        sharey=False
    )

    panel_letters = ["A", "B", "C", "D", "E"]

    bar_color = "#2B6EA6"

    for i, ax in enumerate(axes):

        plot_df = data_list[i]

        x = range(len(plot_df))

        ax.bar(
            x,
            plot_df["Count"],
            color=bar_color,
            edgecolor="black",
            linewidth=0.4
        )

        # remove left/right gaps
        ax.margins(x=0)
        ax.set_xlim(-0.5, len(plot_df) - 0.5)

        ax.set_xticks(list(x))

        # crowded panels rotate automatically
        if len(plot_df) > 25:
            rotation = 90
            tick_size = 6
        else:
            rotation = 0
            tick_size = 8

        ax.set_xticklabels(
            plot_df["XLabel"],
            rotation=rotation,
            ha="center"
        )

        ax.tick_params(
            axis="x",
            labelsize=tick_size,
            pad=1
        )

        ax.tick_params(
            axis="y",
            labelsize=9
        )

        ax.set_yscale("log")

        ax.grid(False)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

        ax.yaxis.set_major_locator(
            mticker.LogLocator(base=10)
        )

        ax.yaxis.set_minor_locator(
            mticker.LogLocator(
                base=10,
                subs=range(2, 10)
            )
        )

        ax.yaxis.set_minor_formatter(
            mticker.NullFormatter()
        )

        ymin = plot_df["Count"].min()
        ymax = plot_df["Count"].max()

        ax.set_ylim(
            max(0.8, ymin * 0.8),
            ymax * 1.3
        )

        # panel label
        ax.text(
            -0.06,
            1.03,
            f"{panel_letters[i]}. bla{families[i]}",
            transform=ax.transAxes,
            fontsize=12,
            fontweight="bold",
            va="bottom",
            ha="left"
        )

    # common labels
    fig.supylabel(
        "Accession Count (log scale)",
        fontsize=14,
        fontweight="bold",
        x=0.015
    )

    fig.supxlabel(
        "Variants",
        fontsize=14,
        fontweight="bold",
        y=0.02
    )

    plt.subplots_adjust(
        left=0.12,
        right=0.98,
        top=0.97,
        bottom=0.08,
        hspace=0.85
    )

    plt.savefig(
        output_png,
        dpi=600,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Panel plot saved to: {output_png}")


def main():

    parser = argparse.ArgumentParser(
        description="Create a 5-panel publication-quality variant plot."
    )

    parser.add_argument(
        "-i",
        "--inputs",
        nargs=5,
        required=True,
        help="Five CSV files"
    )

    parser.add_argument(
        "-f",
        "--families",
        nargs=5,
        required=True,
        help="Family names (e.g. IMP KPC TEM NDM SHV)"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="variant_panel.png",
        help="Output image file"
    )

    args = parser.parse_args()

    make_panel_plot(
        args.inputs,
        args.families,
        args.output
    )


if __name__ == "__main__":
    main()