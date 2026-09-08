#!/usr/bin/env python3
"""
SnapGene .dna -> publication-style genomic feature map.

Drawing rules
-------------
1. CDS/gene features and misc features occupy ONE common genomic band.
2. Every misc feature is a rectangle of the same height, whether or not it
   overlaps a CDS. Misc rectangles are drawn first, so an overlapping CDS is
   visibly on top of them.
3. Promoters are green CDS-style arrows in the same band.
4. blaIMP, blaNDM, blaKPC and blaSHV have fixed contrasting colours, but only
   when the feature itself is classified as a CDS/gene.
5. Parenthetical text inside CDS labels is normal/non-italic; the gene portion
   remains italic.
6. All labels placed below the map are on ONE level, preserve genomic left-to-right
   order, and their leader arrows cannot cross.
"""

from pathlib import Path
import re
import sys
import textwrap

from Bio import SeqIO
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle


# ============================================================
# EDITABLE PARAMETERS
# ============================================================

STYLE = {
    # FONT
    "font_family": "Times New Roman",

    # FIGURE
    "figure_width": 18,
    "figure_height": 7.5,
    "dpi": 600,
    # PNG/JPEG are raster outputs. Keep their renderer size moderate so
    # large plasmids do not cause a Windows MemoryError. SVG remains vector.
    "raster_dpi": 300,
    "transparent": True,

    # GENOMIC BACKBONE
    "backbone_color": "black",
    "backbone_linewidth": 2.5,
    "backbone_y": 0.0,

    # CDS / GENE ARROWS
    "gene_body_height": 0.16,
    "gene_head_height": 0.45,
    "gene_head_fraction": 0.16,
    "gene_color": "#E53935",
    "gene_edge_color": "black",
    "gene_edge_width": 1.0,

    # FIXED COLOURS — ONLY FOR CDS/GENE FEATURES
    "special_cds_colors": {
        # Carbapenemase / beta-lactamase CDS families: fixed green
        "blaIMP": "#008000",
        "blaNDM": "#008000",
        "blaKPC": "#008000",
        "blaSHV": "#008000",
    },

    # CDS LABELS
    "gene_label_color": "black",
    "gene_label_fontsize": 24,
    "gene_label_weight": "bold",
    "gene_label_style": "italic",
    "bracket_label_style": "normal",
    "gene_min_fraction_internal_label": 0.06,

    # MISC FEATURES — ONE COMMON BAND FOR ALL MISC FEATURES
    "misc_height": 0.28,          # slightly taller than CDS body
    "misc_min_fraction": 0.012,
    "misc_color": "#4A90E2",
    "misc_edge_color": "black",
    "misc_edge_width": 0.8,
    "misc_zorder": 2,
    "cds_zorder": 6,
    "promoter_zorder": 5,

    # LABELS / LEADER ARROWS
    "small_label_fontsize": 22,
    "small_label_color": "black",
    "small_label_weight": "bold",
    "small_label_fontstyle": "normal",
    "small_label_wrap_chars": 28,
    "leader_color": "black",
    "leader_width": 0.8,
    "leader_head_size": 5,
    "leader_gap": 0.015,

    # Bottom label layout — ONE LEVEL ONLY
    "label_lane_start": -0.42,
    "label_lane_step": -0.34,
    "label_lane_count": 1,
    "label_horizontal_gap": 0.035,
    "label_min_width_fraction": 0.012,

    # PROMOTER
    "promoter_color": "#90EE90",
    "promoter_head_fraction": 0.16,

    # OUTPUTS
    "write_png": True,
    "write_svg": True,
    "write_jpeg": True,
}


# ============================================================
# LABEL PROCESSING / CLASSIFICATION
# ============================================================

def clean_text(value):
    if isinstance(value, list):
        value = value[0] if value else ""
    value = str(value)
    value = value.replace("\\n", " ")
    return re.sub(r"\s+", " ", value).strip()


def get_feature_label(feature):
    q = feature.qualifiers
    preferred = ["label", "gene", "name", "locus_tag", "product", "note"]
    for key in preferred:
        if key in q:
            value = clean_text(q[key])
            if value:
                return value
    return clean_text(feature.type)


def remove_promoter_word(label):
    text = clean_text(label)
    text = re.sub(
        r"\s*\(\s*promoter(?:\s+region)?\s*\)",
        "",
        text,
        flags=re.IGNORECASE,
    )
    return re.sub(r"\s+", " ", text).strip()


def is_promoter(feature, original_label):
    ftype = str(feature.type).strip().lower()
    if ftype in {"promoter", "promoter_region", "promoter region"}:
        return True
    return bool(
        re.search(
            r"\(\s*promoter(?:\s+region)?\s*\)",
            clean_text(original_label),
            flags=re.IGNORECASE,
        )
    )


def is_gene_like(feature, original_label):
    if is_promoter(feature, original_label):
        return False

    label_lower = clean_text(original_label).lower()

    # These are feature annotations, not CDS/gene arrows.
    if "attc" in label_lower or "core" in label_lower or "att" in label_lower:
        return False

    ftype = str(feature.type).strip().lower()
    if ftype in {"cds", "gene", "orf"}:
        return True

    gene_words = [
        "bla", "aac", "aad", "imp", "ndm", "kpc", "oxa", "shv",
        "tem", "intl", "qnr", "mcr",
    ]
    return any(word in label_lower for word in gene_words)


def get_coordinates(feature):
    start = int(feature.location.start) + 1
    end = int(feature.location.end)
    if end < start:
        start, end = end, start
    return start, end


def get_strand(feature):
    try:
        return -1 if feature.location.strand == -1 else 1
    except Exception:
        return 1


def intervals_overlap(a_start, a_end, b_start, b_end):
    return a_start <= b_end and a_end >= b_start


# ============================================================
# COLOUR / LABEL HELPERS
# ============================================================

def get_special_cds_color(label):
    """Special colour is applied only after the feature is known to be CDS/gene."""
    text = clean_text(label).lower()
    for gene_name, color in STYLE["special_cds_colors"].items():
        if re.search(
            rf"(?<![a-z0-9]){re.escape(gene_name.lower())}(?![a-z0-9])",
            text,
        ):
            return color
    return STYLE["gene_color"]


def split_bracket_label(label):
    """Split a CDS label while keeping every parenthetical block together."""
    pieces = []
    for part in re.split(r"(\([^)]*\))", clean_text(label)):
        if part:
            pieces.append((part, part.startswith("(")))
    return pieces


def draw_mixed_gene_label(ax, x, y, label, fontsize, fontweight):
    """
    Draw one continuous CDS label.

    Text before/after parentheses is italic. Parenthetical text is normal.
    Every piece is placed edge-to-edge, so the parentheses cannot sit on top
    of the neighbouring text.
    """
    pieces = split_bracket_label(label)
    if not pieces:
        return

    if not any(is_bracket for _, is_bracket in pieces):
        ax.text(
            x, y, label,
            ha="center", va="center",
            fontsize=fontsize,
            color=STYLE["gene_label_color"],
            fontweight=fontweight,
            fontstyle=STYLE["gene_label_style"],
            fontfamily=STYLE["font_family"],
            zorder=STYLE["cds_zorder"] + 2,
            clip_on=True,
        )
        return

    # Measure each piece in pixels.
    temp = []
    for text, is_bracket in pieces:
        obj = ax.text(
            0, 0, text,
            ha="left", va="center",
            fontsize=fontsize,
            color=STYLE["gene_label_color"],
            fontweight=fontweight,
            fontstyle=(
                STYLE["bracket_label_style"]
                if is_bracket else STYLE["gene_label_style"]
            ),
            fontfamily=STYLE["font_family"],
            alpha=0,
        )
        temp.append(obj)

    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    widths = [obj.get_window_extent(renderer=renderer).width for obj in temp]
    for obj in temp:
        obj.remove()

    total = sum(widths)
    center_px = ax.transData.transform((x, y))[0]
    y_px = ax.transData.transform((x, y))[1]
    cursor = center_px - total / 2

    for (text, is_bracket), width_px in zip(pieces, widths):
        piece_center = cursor + width_px / 2
        piece_x = ax.transData.inverted().transform(
            (piece_center, y_px)
        )[0]

        ax.text(
            piece_x, y, text,
            ha="center", va="center",
            fontsize=fontsize,
            color=STYLE["gene_label_color"],
            fontweight=fontweight,
            fontstyle=(
                STYLE["bracket_label_style"]
                if is_bracket else STYLE["gene_label_style"]
            ),
            fontfamily=STYLE["font_family"],
            zorder=STYLE["cds_zorder"] + 2,
            clip_on=True,
        )
        cursor += width_px


# ============================================================
# FEATURE DRAWING
# ============================================================

def draw_backbone(ax, seq_len):
    y = STYLE["backbone_y"]
    ax.plot(
        [1, seq_len], [y, y],
        color=STYLE["backbone_color"],
        linewidth=STYLE["backbone_linewidth"],
        solid_capstyle="round",
        zorder=1,
    )


def make_arrow_points(start, end, strand, body_height, head_height,
                      head_fraction):
    """Return points for a CDS-style arrow with a stable, non-spiky head."""
    y = STYLE["backbone_y"]
    width = max(end - start + 1, 1)

    # The head must be visible but cannot consume the complete feature.
    head_width = width * head_fraction
    head_width = min(head_width, width * 0.32)
    head_width = max(head_width, width * 0.12)
    head_width = min(head_width, max(width - 1, width * 0.5))

    if strand >= 0:
        body_end = end - head_width
        points = [
            (start, y - body_height / 2),
            (body_end, y - body_height / 2),
            (body_end, y - head_height / 2),
            (end, y),
            (body_end, y + head_height / 2),
            (body_end, y + body_height / 2),
            (start, y + body_height / 2),
        ]
    else:
        body_start = start + head_width
        points = [
            (end, y - body_height / 2),
            (body_start, y - body_height / 2),
            (body_start, y - head_height / 2),
            (start, y),
            (body_start, y + head_height / 2),
            (body_start, y + body_height / 2),
            (end, y + body_height / 2),
        ]

    return points


def draw_arrow(ax, start, end, strand, facecolor, zorder):
    points = make_arrow_points(
        start, end, strand,
        STYLE["gene_body_height"],
        STYLE["gene_head_height"],
        STYLE["gene_head_fraction"],
    )
    patch = Polygon(
        points,
        closed=True,
        facecolor=facecolor,
        edgecolor=STYLE["gene_edge_color"],
        linewidth=STYLE["gene_edge_width"],
        joinstyle="miter",
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def draw_gene(ax, item, seq_len, labels):
    start, end, strand, label = (
        item["start"], item["end"], item["strand"], item["label"]
    )

    draw_arrow(
        ax, start, end, strand,
        get_special_cds_color(label),
        STYLE["cds_zorder"],
    )

    width = max(end - start + 1, 1)
    fraction = width / seq_len
    required_fraction = len(label) * 0.008

    if fraction >= STYLE["gene_min_fraction_internal_label"] and fraction >= required_fraction:
        draw_mixed_gene_label(
            ax,
            (start + end) / 2,
            STYLE["backbone_y"],
            label,
            STYLE["gene_label_fontsize"],
            STYLE["gene_label_weight"],
        )
    else:
        labels.append({
            "x": (start + end) / 2,
            "feature_y": STYLE["backbone_y"] - STYLE["gene_body_height"] / 2,
            "label": label,
            "fontsize": STYLE["gene_label_fontsize"],
            "fontweight": STYLE["gene_label_weight"],
            "fontstyle": STYLE["gene_label_style"],
            "kind": "gene",
        })


def draw_misc(ax, item, seq_len, labels):
    """
    ALL misc features are drawn at y=backbone_y.

    This is the important correction: there is NO alternate lower row for
    misc features. Overlapping and non-overlapping misc features have exactly
    the same vertical position and height.
    """
    base_y = STYLE["backbone_y"]
    width = max(item["end"] - item["start"] + 1, 1)
    display_width = max(width, seq_len * STYLE["misc_min_fraction"])
    center = (item["start"] + item["end"]) / 2
    x_left = center - display_width / 2

    rect = Rectangle(
        (x_left, base_y - STYLE["misc_height"] / 2),
        display_width,
        STYLE["misc_height"],
        facecolor=STYLE["misc_color"],
        edgecolor=STYLE["misc_edge_color"],
        linewidth=STYLE["misc_edge_width"],
        zorder=STYLE["misc_zorder"],
    )
    ax.add_patch(rect)

    labels.append({
        "x": center,
        "feature_y": base_y - STYLE["misc_height"] / 2,
        "label": item["label"],
        "fontsize": STYLE["small_label_fontsize"],
        "fontweight": STYLE["small_label_weight"],
        "fontstyle": STYLE["small_label_fontstyle"],
        "kind": "misc",
    })


def draw_promoter(ax, item, labels):
    """Promoter uses the same arrow geometry as CDS, only green."""
    draw_arrow(
        ax,
        item["start"], item["end"], item["strand"],
        STYLE["promoter_color"],
        STYLE["promoter_zorder"],
    )

    labels.append({
        "x": (item["start"] + item["end"]) / 2,
        "feature_y": STYLE["backbone_y"] - STYLE["gene_body_height"] / 2,
        "label": item["label"],
        "fontsize": STYLE["small_label_fontsize"],
        "fontweight": STYLE["small_label_weight"],
        "fontstyle": STYLE["small_label_fontstyle"],
        "kind": "promoter",
    })


# ============================================================
# BOTTOM LABEL LAYOUT
# ============================================================

def estimate_label_width(ax, label, fontsize, fontweight, fontstyle):
    """Measure the actual rendered label width in x/data coordinates."""
    label = clean_text(label)
    wrapped = textwrap.fill(
        label,
        width=STYLE["small_label_wrap_chars"],
        break_long_words=False,
        break_on_hyphens=False,
    ) or label

    dummy = ax.text(
        0, 0, wrapped,
        ha="left", va="top",
        fontsize=fontsize,
        fontweight=fontweight,
        fontstyle=fontstyle,
        fontfamily=STYLE["font_family"],
        alpha=0,
    )
    ax.figure.canvas.draw()
    renderer = ax.figure.canvas.get_renderer()
    bbox = dummy.get_window_extent(renderer=renderer)
    dummy.remove()

    # Correct conversion: pixels -> genomic x units using the axes transform.
    px0 = ax.transData.transform((0, 0))[0]
    px1 = px0 + bbox.width
    x0 = ax.transData.inverted().transform((px0, 0))[0]
    x1 = ax.transData.inverted().transform((px1, 0))[0]

    return max(abs(x1 - x0), STYLE["label_min_width_fraction"])

def intervals_overlap_with_gap(left1, right1, left2, right2, gap):
    return not (right1 + gap <= left2 or right2 + gap <= left1)


def _place_in_lane(desired, half_width, occupied, gap):
    """Return a collision-free x position in one lane if possible."""
    left = desired - half_width
    right = desired + half_width

    if not any(
        intervals_overlap_with_gap(
            left, right, old_left, old_right, gap
        )
        for old_left, old_right, _ in occupied
    ):
        return desired

    # Search only to the RIGHT of existing labels.
    # This preserves the left-to-right order of feature labels and prevents
    # leader arrows from crossing.
    candidates = []
    for old_left, old_right, _ in sorted(occupied):
        candidates.append(old_right + gap + half_width)

    for candidate in sorted(candidates, key=lambda v: abs(v - desired)):
        c_left = candidate - half_width
        c_right = candidate + half_width
        if not any(
            intervals_overlap_with_gap(
                c_left, c_right, old_left, old_right, gap
            )
            for old_left, old_right, _ in occupied
        ):
            return candidate

    return None


def layout_bottom_labels(ax, labels, seq_len):
    """
    Place all bottom annotations on ONE single level.

    Labels preserve the left-to-right genomic order of their features.
    They are first placed at their true feature positions and then shifted
    only to the right as needed to avoid overlap. Finally, the complete set
    is translated as a group if necessary so that the labels remain inside
    the plotting range. This prevents crossed leader arrows and also prevents
    labels from expanding the saved figure to an enormous bounding box.
    """
    if not labels:
        return

    pending = []
    for original in labels:
        item = dict(original)
        item["label_width"] = estimate_label_width(
            ax,
            item["label"],
            item["fontsize"],
            item["fontweight"],
            item["fontstyle"],
        )
        item["desired_x"] = item["x"]
        pending.append(item)

    # Genomic left-to-right order is the only ordering used.
    pending.sort(key=lambda z: (z["desired_x"], z["label"]))

    gap = STYLE["label_horizontal_gap"]
    xmin, xmax = ax.get_xlim()

    # Forward pass: never move a later label to the left of an earlier label.
    cursor = xmin
    for item in pending:
        half = item["label_width"] / 2
        desired = item["desired_x"]
        left_limit = xmin + half
        chosen_x = max(desired, left_limit, cursor + gap + half)
        item["label_x"] = chosen_x
        cursor = chosen_x + half

    # If the ordered labels run past the right edge, shift the whole layout
    # left by the smallest amount necessary. This preserves their order.
    rightmost = max(
        item["label_x"] + item["label_width"] / 2
        for item in pending
    )
    overflow = rightmost - xmax
    if overflow > 0:
        for item in pending:
            item["label_x"] -= overflow

    # Draw all labels at exactly the same y-coordinate.
    y_text = STYLE["label_lane_start"]
    for item in pending:
        draw_below_label(
            ax,
            x=item["label_x"],
            feature_x=item["x"],
            feature_y=item["feature_y"],
            label=item["label"],
            y_text=y_text,
            fontsize=item["fontsize"],
            fontweight=item["fontweight"],
            fontstyle=item["fontstyle"],
        )

def draw_below_label(ax, x, feature_x, feature_y, label, y_text, fontsize,
                     fontweight, fontstyle):
    label = clean_text(label)
    wrapped = textwrap.fill(
        label,
        width=STYLE["small_label_wrap_chars"],
        break_long_words=False,
        break_on_hyphens=False,
    )

    # The text can move sideways to avoid another label, but the leader line
    # always starts at the true feature coordinate.
    ax.annotate(
        wrapped,
        xy=(feature_x, feature_y),
        xytext=(x, y_text),
        ha="center",
        va="top",
        fontsize=fontsize,
        color=STYLE["small_label_color"],
        fontweight=fontweight,
        fontstyle=fontstyle,
        fontfamily=STYLE["font_family"],
        arrowprops=dict(
            arrowstyle="-|>",
            color=STYLE["leader_color"],
            lw=STYLE["leader_width"],
            mutation_scale=STYLE["leader_head_size"],
            shrinkA=0,
            shrinkB=0,
            connectionstyle="arc3,rad=0",
        ),
        zorder=30,
    )


# ============================================================
# FEATURE RELATIONSHIPS
# ============================================================

def assign_overlap(feature_list):
    """Mark misc features that overlap a CDS/gene, without changing their y."""
    genes = [item for item in feature_list if item["gene"]]

    for item in feature_list:
        item["overlap_gene"] = False
        if item["promoter"] or item["gene"]:
            continue

        item["overlap_gene"] = any(
            intervals_overlap(
                item["start"], item["end"],
                gene["start"], gene["end"],
            )
            for gene in genes
        )

    return feature_list


# ============================================================
# MAIN PIPELINE
# ============================================================

def create_map(dna_file, output_base):
    dna_file = Path(dna_file)
    if not dna_file.exists():
        raise FileNotFoundError(f"DNA file not found: {dna_file}")

    record = SeqIO.read(str(dna_file), "snapgene")
    seq_len = len(record.seq)
    feature_list = []

    for feature in record.features:
        try:
            start, end = get_coordinates(feature)
        except Exception:
            continue

        if start < 1 or end > seq_len:
            continue

        original_label = get_feature_label(feature)
        promoter = is_promoter(feature, original_label)
        display_label = (
            remove_promoter_word(original_label)
            if promoter else clean_text(original_label)
        )

        if not display_label:
            display_label = "promoter"

        gene = is_gene_like(feature, original_label)

        feature_list.append({
            "feature": feature,
            "start": start,
            "end": end,
            "strand": get_strand(feature),
            "label": display_label,
            "promoter": promoter,
            "gene": gene,
        })

    feature_list.sort(key=lambda x: (x["start"], x["end"]))
    feature_list = assign_overlap(feature_list)

    fig, ax = plt.subplots(
        figsize=(STYLE["figure_width"], STYLE["figure_height"])
    )
    fig.patch.set_alpha(0)
    ax.set_facecolor("none")

    # Room for the single bottom annotation level; the feature map itself
    # remains a single compact horizontal band.
    ax.set_xlim(-0.02 * seq_len, 1.02 * seq_len)
    ax.set_ylim(
        STYLE["label_lane_start"]
        + (STYLE["label_lane_count"] - 1) * STYLE["label_lane_step"]
        - 0.26,
        0.60,
    )

    draw_backbone(ax, seq_len)

    # ========================================================
    # CRITICAL DRAW ORDER
    # ========================================================
    # MISC FIRST: all misc rectangles are at y=0 and therefore form ONE band.
    # CDS SECOND: CDS arrows cover the middle of an overlapping misc rectangle.
    # PROMOTER: same common band, green CDS-style arrow.
    # ========================================================
    misc_items = [
        item for item in feature_list
        if not item["promoter"] and not item["gene"]
    ]
    gene_items = [item for item in feature_list if item["gene"]]
    promoter_items = [item for item in feature_list if item["promoter"]]

    bottom_labels = []

    # ALL misc features — overlapping AND non-overlapping — same y and height.
    for item in misc_items:
        draw_misc(ax, item, seq_len, bottom_labels)

    # Promoters are also in the common genomic band.
    for item in promoter_items:
        draw_promoter(ax, item, bottom_labels)

    # CDS arrows are on top of misc rectangles.
    for item in gene_items:
        draw_gene(ax, item, seq_len, bottom_labels)

    # Place annotations only after ALL features have been drawn, so labels can
    # be globally separated and no two labels are allowed to overlap.
    fig.canvas.draw()
    layout_bottom_labels(ax, bottom_labels, seq_len)

    ax.axis("off")
    plt.subplots_adjust(
        left=0.015,
        right=0.985,
        top=0.98,
        bottom=0.04,
    )

    output_base = Path(output_base)
    if output_base.suffix:
        output_base = output_base.with_suffix("")

    if STYLE["write_png"]:
        fig.savefig(
            output_base.with_suffix(".png"),
            dpi=STYLE["raster_dpi"],
            transparent=True,
            bbox_inches="tight",
            pad_inches=0.08,
        )

    if STYLE["write_svg"]:
        fig.savefig(
            output_base.with_suffix(".svg"),
            transparent=True,
            bbox_inches="tight",
            pad_inches=0.08,
        )

    if STYLE["write_jpeg"]:
        fig.savefig(
            output_base.with_suffix(".jpeg"),
            dpi=STYLE["raster_dpi"],
            facecolor="white",
            bbox_inches="tight",
            pad_inches=0.08,
        )

    plt.close(fig)

    print("\n============================================")
    print(" SnapGene genomic feature map created")
    print("============================================")


def main():
    if len(sys.argv) >= 2:
        dna_file = sys.argv[1]
    else:
        dna_file = input("SnapGene .dna file path: ").strip().strip('"')

    if not dna_file:
        sys.exit(1)

    if len(sys.argv) >= 3:
        output_base = sys.argv[2]
    else:
        default_name = Path(dna_file).stem + "_map"
        output_base = input(f"Output base name [{default_name}]: ").strip()
        if not output_base:
            output_base = default_name

    create_map(dna_file, output_base)


if __name__ == "__main__":
    main()
