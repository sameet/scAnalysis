import cellxgene_census
import pandas as pd
import argparse
import os
import sys

def query_metadata(species, tissue, output_path, tissue_column="tissue"):
    """
    Queries CellxGene Census for cell metadata based on species and tissue.
    Saves the result to a Parquet file.
    """
    print(f"Connecting to CellxGene Census for species: {species}...")
    
    try:
        with cellxgene_census.open_soma() as census:
            print(f"Querying metadata where {tissue_column} == '{tissue}'...")
            
            # get_obs returns a pandas DataFrame
            # We filter for is_primary_data == True to avoid duplicates
            cell_metadata = cellxgene_census.get_obs(
                census,
                organism=species,
                value_filter=f"{tissue_column} == '{tissue}' and is_primary_data == True"
            )
            
            if cell_metadata.empty:
                print(f"No metadata found for {species} in {tissue} ({tissue_column}).")
                return False

            print(f"Found {len(cell_metadata)} cells.")
            print(f"Saving to {output_path}...")
            
            # Ensure the output directory exists
            os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
            
            cell_metadata.to_parquet(output_path, index=False)
            print("Successfully saved metadata.")
            return True
            
    except Exception as e:
        print(f"Error during query: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Query CellxGene Atlas for cell metadata.")
    parser.add_argument("--species", type=str, required=True, help="Species name (e.g., 'homo_sapiens' or 'mus_musculus')")
    parser.add_argument("--tissue", type=str, required=True, help="Tissue name (e.g., 'lung')")
    parser.add_argument("--output", type=str, default="cell_metadata.parquet", help="Output Parquet file path")
    parser.add_argument("--use-general-tissue", action="store_true", help="Use 'tissue_general' instead of 'tissue' for filtering")

    args = parser.parse_args()

    tissue_col = "tissue_general" if args.use_general_tissue else "tissue"
    
    success = query_metadata(args.species, args.tissue, args.output, tissue_col)
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
