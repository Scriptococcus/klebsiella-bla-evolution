#!/usr/bin/env python3
"""
compile_png_maps_by_gene.py

READ EXISTING PNG MAPS AND COMPILE THEM WITHOUT REDRAWING THE MAPS.

This script is intentionally separate from snapgene_feature_map_compiled.py.

It does NOT read .dna files.
It does NOT redraw arrows, genes, labels, or leader lines.

It reads the PNG files that you already generated, keeps each PNG as one
complete image, proportionally scales it, and stacks the images vertically
inside five fixed gene columns:

    blaIMP | blaNDM | blaKPC | blaSHV | blaTEM

Each original PNG is therefore visually unchanged except for proportional
scaling/cropping of only the OUTER empty margin.

Usage:

    python scripts/13_compile_genomic_feature_maps.py path/to/snapgene_images1

Optional output directory:

    python scripts/13_compile_genomic_feature_maps.py \
        path/to/snapgene_images1 \
        path/to/compiled_output

Outputs:

    compiled_gene_maps.png
    compiled_gene_maps.svg
    compiled_gene_maps.pdf
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager


# ============================================================
# SETTINGS
# ============================================================

GENES = ["blaIMP", "blaNDM", "blaKPC", "blaSHV", "blaTEM"]

SETTINGS = {
    # Final master image width.
    "canvas_width_px": 5000,

    # White space inside each column around the imported PNG.
    "column_padding_px": 45,

    # Vertical spacing between imported PNG maps in the same column.
    "map_gap_px": 55,

    # Top area reserved for the gene headings.
    "heading_height_px": 95,

    # Heading font size.
    "heading_font_size": 36,

    # Separator lines.
    "separator_width_px": 2,
    "separator_color": (130, 130, 130, 255),

    # White background.
    "background": (255, 255, 255, 255),

    # Crop only outer near-white/transparent margin.
    "trim_outer_margin": True,

    # Pixel threshold used to decide whether a pixel is background white.
    "white_threshold": 250,

    # Output filenames.
    "png_name": "compiled_gene_maps.png",
    "svg_name": "compiled_gene_maps.svg",
    "pdf_name": "compiled_gene_maps.pdf",
}


# ============================================================
# ROBUST GENE DETECTION FROM PNG FILENAMES
# ============================================================

def detect_gene(filename: str) -> str | None:
    """
    Detect gene family from names such as:

        IMP-1 AP018454_map.png
        IMP1CP141568_map.png
        KPC2CP015387_map.png
        ndm15mk372392_map.png
        SHV12CP170854_map.png
        tem1cp144513_map.png

    We deliberately use substring-style matching after removing punctuation
    because names like IMP1CP141568 do not contain word boundaries.
    """
    stem = Path(filename).stem.lower()

    # Remove the common trailing _map.
    stem = re.sub(r"_?map$", "", stem)

    # Keep letters and digits only.
    s = re.sub(r"[^a-z0-9]", "", stem)

    # Match the distinctive beta-lactamase family token.
    # Order is explicit and exact for the five requested columns.
    patterns = [
        ("blaIMP", r"(?:bla)?imp"),
        ("blaNDM", r"(?:bla)?ndm"),
        ("blaKPC", r"(?:bla)?kpc"),
        ("blaSHV", r"(?:bla)?shv"),
        ("blaTEM", r"(?:bla)?tem"),
    ]

    for gene, pattern in patterns:
        if re.search(pattern, s):
            return gene

    return None


# ============================================================
# IMAGE PROCESSING
# ============================================================

def trim_outer_margin(img: Image.Image) -> Image.Image:
    """
    Remove only unnecessary OUTER whitespace.

    The actual map and all labels remain in the PNG.
    """
    img = img.convert("RGBA")

    # --------------------------------------------------------
    # Transparent PNG
    # --------------------------------------------------------
    alpha = img.getchannel("A")

    # If there is any transparency, use the visible alpha bounding box.
    amin, amax = alpha.getextrema()
    if amin < 255:
        bbox = alpha.getbbox()
        if bbox is not None:
            return img.crop(bbox)

    # --------------------------------------------------------
    # White-background PNG
    # --------------------------------------------------------
    # Detect pixels that are meaningfully darker than white.
    rgb = img.convert("RGB")
    px = rgb.load()

    left = img.width
    top = img.height
    right = -1
    bottom = -1

    threshold = SETTINGS["white_threshold"]

    for y in range(img.height):
        for x in range(img.width):
            r, g, b = px[x, y]

            # Pixel is considered content when it is not near-white.
            if (
                r < threshold
                or g < threshold
                or b < threshold
            ):
                if x < left:
                    left = x
                if x > right:
                    right = x
                if y < top:
                    top = y
                if y > bottom:
                    bottom = y

    if right >= left and bottom >= top:
        return img.crop((left, top, right + 1, bottom + 1))

    return img


def load_map(path: Path) -> Image.Image:
    img = Image.open(path).convert("RGBA")

    if SETTINGS["trim_outer_margin"]:
        img = trim_outer_margin(img)

    return img


def scale_proportionally(
    img: Image.Image,
    target_width: int,
) -> Image.Image:
    """
    Scale the WHOLE PNG proportionally.
    """
    if img.width <= 0 or img.height <= 0:
        raise ValueError(f"Invalid PNG dimensions for {img}")

    ratio = target_width / img.width

    new_width = max(1, int(round(img.width * ratio)))
    new_height = max(1, int(round(img.height * ratio)))

    return img.resize(
        (new_width, new_height),
        Image.Resampling.LANCZOS,
    )


# ============================================================
# FIND INPUT PNG FILES
# ============================================================

def find_pngs(input_path: Path) -> list[Path]:
    if not input_path.exists():
        raise FileNotFoundError(
            f"Input folder/file not found: {input_path}"
        )

    if input_path.is_file():
        if input_path.suffix.lower() != ".png":
            raise ValueError(
                f"Input file must be PNG: {input_path}"
            )
        return [input_path]

    # Ignore previously generated compiled outputs.
    pngs = [
        p
        for p in input_path.iterdir()
        if p.is_file()
        and p.suffix.lower() == ".png"
        and not p.name.lower().startswith("compiled_gene_maps")
    ]

    pngs.sort(key=lambda p: p.name.lower())

    if not pngs:
        raise FileNotFoundError(
            f"No individual PNG map files found in: {input_path}"
        )

    return pngs


def group_by_gene(
    pngs: list[Path],
) -> tuple[dict[str, list[Path]], list[Path]]:
    groups = {gene: [] for gene in GENES}
    unassigned = []

    for png in pngs:
        gene = detect_gene(png.name)

        if gene is None:
            unassigned.append(png)
        else:
            groups[gene].append(png)

    for gene in GENES:
        groups[gene].sort(key=lambda p: p.name.lower())

    return groups, unassigned


# ============================================================
# FONT
# ============================================================

def load_heading_font(size: int):
    """Return a bold-italic serif font without relying on machine-specific paths."""
    families = [
        "Times New Roman",
        "Times",
        "Liberation Serif",
        "DejaVu Serif",
    ]
    for family in families:
        try:
            path = font_manager.findfont(
                font_manager.FontProperties(family=family, weight="bold", style="italic"),
                fallback_to_default=False,
            )
            if path and Path(path).exists():
                return ImageFont.truetype(path, size)
        except Exception:
            continue
    return ImageFont.load_default()


# ============================================================
# CREATE COMPILED IMAGE
# ============================================================

def create_compiled_image(
    groups: dict[str, list[Path]],
) -> Image.Image:
    canvas_width = SETTINGS["canvas_width_px"]
    ncols = len(GENES)
    column_width = canvas_width / ncols

    usable_column_width = int(
        column_width - 2 * SETTINGS["column_padding_px"]
    )

    # --------------------------------------------------------
    # Load and scale complete PNGs.
    # --------------------------------------------------------
    scaled = {gene: [] for gene in GENES}

    for gene in GENES:
        for png in groups[gene]:
            img = load_map(png)
            img = scale_proportionally(
                img,
                usable_column_width,
            )
            scaled[gene].append((png, img))

    # --------------------------------------------------------
    # Determine required canvas height from tallest column.
    # --------------------------------------------------------
    heading_height = SETTINGS["heading_height_px"]
    map_gap = SETTINGS["map_gap_px"]

    column_heights = []

    for gene in GENES:
        h = heading_height

        if scaled[gene]:
            h += sum(
                img.height
                for _, img in scaled[gene]
            )
            h += map_gap * (len(scaled[gene]) - 1)

        column_heights.append(h)

    canvas_height = max(column_heights) + 50

    canvas = Image.new(
        "RGBA",
        (
            canvas_width,
            canvas_height,
        ),
        SETTINGS["background"],
    )

    # --------------------------------------------------------
    # Paste each WHOLE PNG down its appropriate column.
    # --------------------------------------------------------
    for i, gene in enumerate(GENES):
        col_left = int(round(i * column_width))
        col_right = int(round((i + 1) * column_width))
        col_w = col_right - col_left

        y = heading_height

        for _, img in scaled[gene]:
            x = col_left + (col_w - img.width) // 2

            # Entire map PNG is inserted as one intact image.
            canvas.alpha_composite(img, (x, y))

            y += img.height + map_gap

    # --------------------------------------------------------
    # Draw separator lines between columns.
    # --------------------------------------------------------
    draw = ImageDraw.Draw(canvas)

    for i in range(1, ncols):
        x = int(round(i * column_width))

        draw.rectangle(
            [
                x - SETTINGS["separator_width_px"] // 2,
                10,
                x + SETTINGS["separator_width_px"] // 2,
                canvas_height - 10,
            ],
            fill=SETTINGS["separator_color"],
        )

    # --------------------------------------------------------
    # Column headings: italic gene names.
    # --------------------------------------------------------
    font = load_heading_font(
        SETTINGS["heading_font_size"]
    )

    for i, gene in enumerate(GENES):
        col_left = int(round(i * column_width))
        col_right = int(round((i + 1) * column_width))
        col_w = col_right - col_left

        bbox = draw.textbbox(
            (0, 0),
            gene,
            font=font,
        )
        text_w = bbox[2] - bbox[0]
        text_h = bbox[3] - bbox[1]

        x = col_left + (col_w - text_w) // 2
        y = max(5, (heading_height - text_h) // 2 - 2)

        draw.text(
            (x, y),
            gene,
            fill=(0, 0, 0, 255),
            font=font,
        )

    return canvas


# ============================================================
# SAVE PNG + SVG + PDF
# ============================================================

def save_outputs(
    canvas: Image.Image,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    png_path = output_dir / SETTINGS["png_name"]
    svg_path = output_dir / SETTINGS["svg_name"]
    pdf_path = output_dir / SETTINGS["pdf_name"]

    # --------------------------------------------------------
    # PNG
    # --------------------------------------------------------
    canvas.save(
        png_path,
        "PNG",
        optimize=True,
    )

    # --------------------------------------------------------
    # SVG + PDF
    #
    # Important: we place the already-composed master PNG as a
    # single image. Therefore matplotlib cannot reflow labels.
    # --------------------------------------------------------
    rgb = Image.new(
        "RGB",
        canvas.size,
        (255, 255, 255),
    )
    rgb.paste(
        canvas,
        mask=canvas.getchannel("A"),
    )

    dpi = 150
    fig_width = rgb.width / dpi
    fig_height = rgb.height / dpi

    fig = plt.figure(
        figsize=(fig_width, fig_height),
        dpi=dpi,
        facecolor="white",
    )

    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(
        rgb,
        interpolation="none",
    )
    ax.set_xlim(0, rgb.width)
    ax.set_ylim(rgb.height, 0)
    ax.axis("off")

    fig.savefig(
        svg_path,
        format="svg",
        dpi=dpi,
        facecolor="white",
        edgecolor="none",
        bbox_inches=None,
        pad_inches=0,
    )

    fig.savefig(
        pdf_path,
        format="pdf",
        dpi=dpi,
        facecolor="white",
        edgecolor="none",
        bbox_inches=None,
        pad_inches=0,
    )

    plt.close(fig)

    print("\nCreated:")
    print(f"  {png_path}")
    print(f"  {svg_path}")
    print(f"  {pdf_path}")


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    if len(sys.argv) < 2:
        raw = input(
            "Folder containing existing PNG maps: "
        ).strip().strip('"')
        input_path = Path(raw)
    else:
        input_path = Path(sys.argv[1])

    if len(sys.argv) >= 3:
        output_dir = Path(sys.argv[2])
    else:
        output_dir = input_path.parent / "compiled_output"

    pngs = find_pngs(input_path)

    print(f"\nFound {len(pngs)} individual PNG file(s).\n")

    groups, unassigned = group_by_gene(pngs)

    print("Files assigned to columns:")
    for gene in GENES:
        print(f"\n{gene}: {len(groups[gene])}")
        for p in groups[gene]:
            print(f"    {p.name}")

    if unassigned:
        print("\nWARNING - could not identify a gene for:")
        for p in unassigned:
            print(f"    {p.name}")
        print(
            "\nThose files are excluded from the compiled figure."
        )

    total_assigned = sum(
        len(groups[g]) for g in GENES
    )

    if total_assigned == 0:
        raise RuntimeError(
            "No PNG files could be assigned to the five gene columns."
        )

    canvas = create_compiled_image(groups)
    save_outputs(canvas, output_dir)

    print("\nDone.")
    print(
        f"Final canvas size: "
        f"{canvas.width} x {canvas.height} px"
    )


if __name__ == "__main__":
    main()
