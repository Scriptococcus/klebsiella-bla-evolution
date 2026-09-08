import pandas as pd
import numpy as np
import random
import re
import argparse
import time
from Bio import Phylo

# Try importing geopy for city/state to country resolution
try:
    from geopy.geocoders import Nominatim
    from geopy.exc import GeocoderTimedOut, GeocoderUnavailable
    GEOPY_AVAILABLE = True
except ImportError:
    GEOPY_AVAILABLE = False
    print("WARNING: 'geopy' is not installed. Country resolution from cities/states will be limited.")
    print("Please run: pip install geopy")

# Pre-defined distinct colors for variants and generic properties.
# NOTE: Yellow (#ffe119) and Blue (#4363d8) have been EXPLICITLY REMOVED 
# from this list so they are never used for variants.
DISTINCT_COLORS = [
    "#e6194b", "#3cb44b", "#f58231", "#911eb4", "#46f0f0", 
    "#f032e6", "#bcf60c", "#fabebe", "#008080", "#e6beff", 
    "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", 
    "#ffd8b1", "#000075", "#808080", "#a9a9a9", "#ffffff"
]

# Dedicated colors for Replicon Types so they don't overlap with variants
REPLICON_COLORS = {
    "plasmid": "#ffe119",    # Yellow
    "chromosome": "#4363d8"  # Blue
}

GEO_CACHE = {}

def get_country_from_geopy(location_str):
    """Uses Geopy to resolve a city/state to its root country."""
    if not GEOPY_AVAILABLE:
        return location_str
        
    if location_str in GEO_CACHE:
        return GEO_CACHE[location_str]
        
    geolocator = Nominatim(user_agent="itol_tree_annotator_script")
    try:
        time.sleep(1)
        location = geolocator.geocode(location_str, language='en', timeout=10)
        if location:
            country = location.address.split(",")[-1].strip()
            GEO_CACHE[location_str] = country
            return country
    except (GeocoderTimedOut, GeocoderUnavailable):
        pass
        
    GEO_CACHE[location_str] = location_str 
    return location_str

def clean_year(val):
    val_str = str(val).strip()
    if val_str.lower() in ['unknown', 'nan', 'none', '', 'na', 'missing']:
        return "unknown"
    match = re.search(r'\b(19\d{2}|20\d{2})\b', val_str)
    if match:
        return match.group(1)
    return "unknown"

def clean_country(val):
    val_str = str(val).strip()
    if val_str.lower() in ['unknown', 'nan', 'none', '', 'na', 'missing']:
        return "unknown", False, val_str
    
    original = val_str
    if ":" in val_str:
        val_str = val_str.split(":")[0].strip()
    elif "," in val_str:
        val_str = val_str.split(",")[0].strip()
        
    resolved_country = get_country_from_geopy(val_str)
    was_modified = (original != resolved_country)
    return resolved_country, was_modified, val_str

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
        if val_lower == 'unknown':
            color_map[val] = "#000000" 
        elif val in singletons:
            color_map[val] = "#000000" 
        elif is_replicon and val_lower in REPLICON_COLORS:
            color_map[val] = REPLICON_COLORS[val_lower]
        else:
            if color_idx < len(DISTINCT_COLORS):
                color_map[val] = DISTINCT_COLORS[color_idx]
                color_idx += 1
            else:
                color_map[val] = "#{:06x}".format(random.randint(0, 0xFFFFFF))
    return color_map

def main():
    parser = argparse.ArgumentParser(description="Generate robust iTOL annotations directly from a single Excel sheet.")
    parser.add_argument("--tree", required=True, help="Path to the input .treefile")
    parser.add_argument("--metadata", required=True, help="Path to the metadata Excel file")
    parser.add_argument("--gene", required=True, help="The exact sheet name in Excel to read (e.g., IMP)")
    parser.add_argument("--id_col", default="fasta_id", help="The column name containing tree leaf IDs")
    args = parser.parse_args()

    print("--- Step 1: Parsing Tree Leaf Names ---")
    try:
        tree = Phylo.read(args.tree, "newick")
        tree_taxa = [leaf.name for leaf in tree.get_terminals()]
        print(f"Success: Found {len(tree_taxa)} taxa in the tree file.")
    except Exception as e:
        print(f"CRITICAL ERROR reading tree file: {e}")
        return

    print(f"\n--- Step 2: Reading Excel Sheet: {args.gene} ---")
    try:
        # Strictly reading ONLY the sheet provided by --gene
        df = pd.read_excel(args.metadata, sheet_name=args.gene)
        print(f"Success: Loaded sheet '{args.gene}' with {len(df)} rows.")
    except Exception as e:
        print(f"CRITICAL ERROR reading Excel sheet '{args.gene}': {e}")
        return
    
    if args.id_col not in df.columns:
        print(f"CRITICAL ERROR: ID column '{args.id_col}' not found in the sheet.")
        return

    df[args.id_col] = df[args.id_col].astype(str).str.strip()
    tree_name_mapping = {name.strip("'\" "): name for name in tree_taxa}
    
    def find_tree_match(excel_id):
        excel_id_clean = str(excel_id).strip("'\" ")
        if excel_id_clean in tree_name_mapping:
            return tree_name_mapping[excel_id_clean]
        if "|" in excel_id_clean:
            parts = excel_id_clean.split("|")
            if len(parts) >= 2:
                short_id = f"{parts[0]}|{parts[1]}"
                if short_id in tree_name_mapping:
                    return tree_name_mapping[short_id]
        return np.nan

    df['_tree_match_id'] = df[args.id_col].apply(find_tree_match)
    df = df.dropna(subset=['_tree_match_id'])
    df = df.drop_duplicates(subset=['_tree_match_id'])
    
    print(f"Success: Dynamically matched {len(df)} rows to your phylogenetic tree nodes!")

    if len(df) == 0:
        print("\nCRITICAL ERROR: No matching nodes found.")
        return
    
    print("\n--- Step 3: Determining Best Country Column & Cleaning Data ---")
    
    country_col = "source_country"
    if "source_geo_loc_name" in df.columns:
        geo_count = df["source_geo_loc_name"].dropna().astype(str).str.strip().loc[lambda s: s != ""].count()
        country_count = df["source_country"].dropna().astype(str).str.strip().loc[lambda s: s != ""].count() if "source_country" in df.columns else 0
        
        if geo_count >= country_count:
            country_col = "source_geo_loc_name"
            print(f" -> Auto-detected data layout: Using '{country_col}'.")
        else:
            print(f" -> Auto-detected data layout: Using standard '{country_col}' column.")
    
    FEATURE_MAPPING = {
        "variant": "variant",
        "replicon_type": "replicon_type",
        "collection_year": "source_collection_date",
        "collection_country": country_col
    }
    
    # Process Collection Year
    year_col = FEATURE_MAPPING["collection_year"]
    if year_col in df.columns:
        df[year_col] = df[year_col].apply(clean_year)
    
    # Process Country
    if country_col in df.columns:
        print("\n-> Geocoding and cleaning locations (this may take a moment to avoid rate limits)...")
        cleaned_countries = []
        for idx, row in df.iterrows():
            orig_val = row[country_col]
            pure_country, modified, _ = clean_country(orig_val)
            cleaned_countries.append(pure_country)
        df[country_col] = cleaned_countries

    print("\n--- Step 4: Generating iTOL Colorstrip Profiles ---")
    for display_name, actual_col in FEATURE_MAPPING.items():
        if actual_col not in df.columns:
            print(f"Warning: Column '{actual_col}' not found in Excel. Skipping.")
            continue
            
        df[actual_col] = df[actual_col].replace([np.nan, None, "", "nan", "NaN"], "unknown")
        df[actual_col] = df[actual_col].astype(str).str.replace(",", ";").str.strip()
        
        value_counts = df[actual_col].value_counts()
        unique_values_ordered = value_counts.index.tolist()
        
        is_replicon = (display_name == "replicon_type")
        
        singletons = []
        if display_name == "variant":
            singletons = value_counts[value_counts == 1].index.tolist()
            if singletons:
                print(f" -> Note: Found {len(singletons)} singleton variants. They will be colored black.")
        
        color_map = get_color_map(unique_values_ordered, singletons, is_replicon=is_replicon)
        
        legend_vals = [v for v in color_map.keys() if v not in singletons]
        if singletons:
            legend_vals.append("Singleton Variants") 
            color_map["Singleton Variants"] = "#000000"
            
        legend_shapes = ",".join(["1"] * len(legend_vals))  
        legend_colors = ",".join([color_map[val] for val in legend_vals])
        legend_labels = ",".join([val for val in legend_vals])
        
        if display_name == "variant":
            legend_title = "Variants"
        elif display_name == "replicon_type":
            legend_title = "Replicon types"
        else:
            legend_title = f"{args.gene} {display_name.replace('_', ' ').capitalize()}"
        
        output_filename = f"itol_{args.gene}_{display_name}.txt"
        with open(output_filename, "w") as f:
            f.write("DATASET_COLORSTRIP\n")
            f.write("SEPARATOR COMMA\n")
            f.write(f"DATASET_LABEL,{args.gene}_{display_name.upper()}\n")
            f.write("COLOR,#ff5722\n")
            
            f.write(f"LEGEND_TITLE,{legend_title}\n")
            f.write(f"LEGEND_SHAPES,{legend_shapes}\n")
            f.write(f"LEGEND_COLORS,{legend_colors}\n")
            f.write(f"LEGEND_LABELS,{legend_labels}\n")
            f.write("DATA\n")
            
            for _, row in df.iterrows():
                taxon = row['_tree_match_id']  
                val = row[actual_col]
                color = color_map[val] if val not in singletons else color_map["Singleton Variants"]
                f.write(f"{taxon},{color},{val}\n")
                
        print(f" => Created Annotation File: {output_filename}")

    print("\nAll done! Drag and drop your newly generated text files directly onto your iTOL tree interface.")

if __name__ == "__main__":
    main()