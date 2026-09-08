import pandas as pd
import numpy as np
import random
import re
import argparse
from Bio import Phylo

# Try importing geopy for optional location resolving
try:
    from geopy.geocoders import Nominatim
    from geopy.exc import GeocoderTimedOut, GeocoderUnavailable
    GEOPY_AVAILABLE = True
except ImportError:
    GEOPY_AVAILABLE = False

# High-contrast, vibrant palette for variants.
# EXCLUDED: Bright Blue (#4363d8), Bright Yellow (#ffe119), and Black (#000000)
DISTINCT_COLORS = [
    "#e6194b", "#3cb44b", "#f58231", "#911eb4", "#46f0f0", 
    "#f032e6", "#bcf60c", "#008080", "#e6beff", "#9a6324", 
    "#800000", "#aaffc3", "#808000", "#ff7f50", "#ff1493", 
    "#8a2be2", "#ff8c00", "#00ff7f", "#dc143c", "#40e0d0", 
    "#dda0dd", "#fa8072", "#2f4f4f", "#a0522d", "#8b0000"
]

# Strict colors for Replicon Types
REPLICON_COLORS = {
    "plasmid": "#ffe119",    # Yellow
    "chromosome": "#4363d8", # Blue
    "unknown": "#000000"     # Black
}

def extract_variant_from_tip(tip_name):
    """Extracts the variant name directly from the tree leaf label."""
    tip_str = str(tip_name).strip("'\" ")
    match = re.search(r'\b([A-Za-z]{2,6}-\d+[a-zA-Z0-9]*)\b', tip_str)
    if match:
        return match.group(1)
    return "Unknown"

def get_color_map(ordered_values, singletons=None, is_replicon=False):
    if singletons is None: 
        singletons = []
        
    color_map = {}
    color_idx = 0
    
    sorted_vals = [str(v) for v in ordered_values if str(v).lower() != 'unknown']
    if 'unknown' in [str(v).lower() for v in ordered_values]:
        sorted_vals.append('unknown')

    for val in sorted_vals:
        val_lower = str(val).lower()
        if is_replicon:
            if val_lower == "plasmid":
                color_map[val] = REPLICON_COLORS["plasmid"]
            elif val_lower == "chromosome":
                color_map[val] = REPLICON_COLORS["chromosome"]
            else:
                color_map[val] = REPLICON_COLORS["unknown"] 
        else:
            if val_lower == 'unknown' or val in singletons:
                color_map[val] = "#000000" # Black for singletons & unknown
            else:
                # Assign high-contrast distinct color
                if color_idx < len(DISTINCT_COLORS):
                    color_map[val] = DISTINCT_COLORS[color_idx]
                    color_idx += 1
                else:
                    # Fallback random distinct hex color if variants exceed pre-defined list
                    color_map[val] = "#{:06x}".format(random.randint(0, 0xFFFFFF))
    return color_map

def get_core_accession(id_str):
    """Strips version numbers (.1) and NZ_ prefixes to find core accession IDs."""
    s = str(id_str)
    s = re.sub(r'\.\d+', '', s) 
    s = s.replace("NZ_", "")     
    return s.strip()

def main():
    parser = argparse.ArgumentParser(description="Generate SHV iTOL annotations with ultra-aggressive ID matching.")
    parser.add_argument("--tree", required=True, help="Path to the input .treefile")
    parser.add_argument("--metadata", required=True, help="Path to the metadata Excel file")
    parser.add_argument("--gene", required=True, help="The exact sheet name in Excel to read (e.g., SHV)")
    parser.add_argument("--id_col", default="accession", help="The column name containing tree leaf IDs")
    args = parser.parse_args()

    print("--- Step 1: Parsing Tree Leaf Names for Variants ---")
    try:
        tree = Phylo.read(args.tree, "newick")
        tree_taxa = [leaf.name for leaf in tree.get_terminals()]
        print(f"Success: Found {len(tree_taxa)} taxa in the tree file.")
    except Exception as e:
        print(f"CRITICAL ERROR reading tree file: {e}")
        return

    tree_variants = {}
    for taxon in tree_taxa:
        tree_variants[taxon] = extract_variant_from_tip(taxon)
    
    print(f"\n--- Step 2: Reading Excel Sheet: {args.gene} for Replicon Types ---")
    try:
        df = pd.read_excel(args.metadata, sheet_name=args.gene)
        print(f"Success: Loaded sheet '{args.gene}' with {len(df)} rows.")
    except Exception as e:
        print(f"CRITICAL ERROR reading Excel sheet '{args.gene}': {e}")
        return
    
    tree_name_mapping = {name.strip("'\" "): name for name in tree_taxa}
    
    def find_tree_match(excel_id):
        excel_id_clean = str(excel_id).strip("'\" ")
        if excel_id_clean == 'nan' or not excel_id_clean:
            return np.nan
        if excel_id_clean in tree_name_mapping:
            return tree_name_mapping[excel_id_clean]
        for t_id, original_t_name in tree_name_mapping.items():
            if excel_id_clean in t_id or t_id in excel_id_clean:
                return original_t_name
        ex_core = get_core_accession(excel_id_clean)
        for t_id, original_t_name in tree_name_mapping.items():
            t_core = get_core_accession(t_id)
            if len(ex_core) > 4 and len(t_core) > 4: 
                if ex_core in t_core or t_core in ex_core:
                    return original_t_name
        return np.nan

    # Check for case-insensitive column matching if needed
    matched_col = None
    for col in df.columns:
        if col.strip().lower() == args.id_col.strip().lower():
            matched_col = col
            break

    if matched_col:
        df['_tree_match_id'] = df[matched_col].apply(find_tree_match)
        df_matched = df.dropna(subset=['_tree_match_id']).drop_duplicates(subset=['_tree_match_id'])
        print(f"Success: Dynamically matched {len(df_matched)} rows to your phylogenetic tree nodes!")
        
        replicon_col = "replicon_type"
        replicon_map = {}
        if replicon_col in df_matched.columns:
            for _, row in df_matched.iterrows():
                val = str(row[replicon_col]).strip()
                replicon_map[row['_tree_match_id']] = val if val.lower() not in ['nan', 'none', ''] else "unknown"
        else:
            print(f"Warning: '{replicon_col}' not found in sheet.")
    else:
        print(f"CRITICAL ERROR: Column '{args.id_col}' not found in Excel sheet. Available columns: {list(df.columns)}")
        return

    print("\n--- Step 3: Generating iTOL Colorstrip Profiles ---")
    
    # --- 3A: REPLICON TYPE EXPORT ---
    replicon_tuples = []
    for taxon in tree_taxa:
        replicon_val = replicon_map.get(taxon, "unknown")
        replicon_tuples.append((taxon, replicon_val))
        
    unique_replicons = list(set([val for _, val in replicon_tuples]))
    rep_color_map = get_color_map(unique_replicons, is_replicon=True)
    
    rep_legend_vals = list(rep_color_map.keys())
    rep_shapes = ",".join(["1"] * len(rep_legend_vals))  
    rep_colors = ",".join([rep_color_map[v] for v in rep_legend_vals])
    rep_labels = ",".join(rep_legend_vals)
    
    with open(f"itol_{args.gene}_replicon_type.txt", "w") as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\n")
        f.write(f"DATASET_LABEL,{args.gene}_REPLICON\nCOLOR,#000000\n")
        f.write(f"LEGEND_TITLE,Replicon types\nLEGEND_SHAPES,{rep_shapes}\n")
        f.write(f"LEGEND_COLORS,{rep_colors}\nLEGEND_LABELS,{rep_labels}\nDATA\n")
        for taxon, val in replicon_tuples:
            f.write(f"{taxon},{rep_color_map[val]},{val}\n")
    print(f" => Created: itol_{args.gene}_replicon_type.txt")

    # --- 3B: VARIANT EXPORT ---
    variant_series = pd.Series(list(tree_variants.values()))
    value_counts = variant_series.value_counts()
    
    # Identify singletons (count == 1)
    singletons = value_counts[value_counts == 1].index.tolist()
    unique_variants_ordered = value_counts.index.tolist()
    
    full_var_color_map = get_color_map(unique_variants_ordered, singletons, is_replicon=False)
    
    # Clean up legend list
    legend_variants = [v for v in unique_variants_ordered if v not in singletons and v.lower() != 'unknown']
    
    if singletons:
        legend_variants.append("Singleton Variants")
        full_var_color_map["Singleton Variants"] = "#000000"
    if 'Unknown' in unique_variants_ordered or 'unknown' in unique_variants_ordered:
        legend_variants.append("Unknown")
        full_var_color_map["Unknown"] = "#000000"

    leg_shapes = ",".join(["1"] * len(legend_variants))
    leg_colors = ",".join([full_var_color_map.get(v, "#000000") for v in legend_variants])
    leg_labels = ",".join(legend_variants)
    
    with open(f"itol_{args.gene}_variant.txt", "w") as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR COMMA\n")
        f.write(f"DATASET_LABEL,{args.gene}_VARIANT\nCOLOR,#000000\n")
        f.write(f"LEGEND_TITLE,Variants\nLEGEND_SHAPES,{leg_shapes}\n")
        f.write(f"LEGEND_COLORS,{leg_colors}\nLEGEND_LABELS,{leg_labels}\nDATA\n")
        
        for taxon, val in tree_variants.items():
            if val in singletons:
                display_val = "Singleton Variants"
                color = "#000000"
            elif val.lower() == 'unknown':
                display_val = "Unknown"
                color = "#000000"
            else:
                display_val = val
                color = full_var_color_map.get(val, "#000000")
                
            f.write(f"{taxon},{color},{display_val}\n")
            
    print(f" => Created: itol_{args.gene}_variant.txt")
    print("\nAll done! Drag and drop both generated text files into iTOL.")

if __name__ == "__main__":
    main()