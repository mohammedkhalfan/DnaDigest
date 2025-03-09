import argparse
import pandas as pd
import pickle

def parse_gff(gff_file, feature_preference=['gene', 'mRNA'], source_preference=None):
    """Parses a GFF3 file using pandas and returns a Pandas DataFrame."""

    try:
        # Step 1: Read the GFF file into a DataFrame.
        gff_data = pd.read_csv(gff_file, sep='\t', comment='#', header=None,
                               names=['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

        # Step 2: Filter for preferred features and sources.
        feature_data = gff_data[gff_data['type'].isin(feature_preference)].copy()

        # Step 3: Optionally filter by source preference.
        if source_preference:
            feature_data = feature_data[feature_data['source'].isin(source_preference)]

        # Step 4: Extract feature information.
        feature_data['id'] = feature_data['attributes'].str.extract(r'ID=([^;]+)')
        feature_data['name'] = feature_data['attributes'].str.extract(r'Name=([^;]+)')
        feature_data['parent_id'] = feature_data['attributes'].str.extract(r'Parent=([^;]+)')

        # Step 5: Calculate feature length.
        feature_data['length'] = abs(feature_data['end'] - feature_data['start']) + 1

         # Step 7: Select and order the desired columns.
        feature_data = feature_data[['id', 'name', 'parent_id', 'seqid', 'start', 'end', 'length', 'strand', 'source', 'type']]
        feature_data = feature_data.reset_index(drop=True)

        return feature_data

    except FileNotFoundError:
        print(f"Error: GFF file '{gff_file}' not found.")
        return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gene information from a GFF3 file.")
    parser.add_argument("gff_file", help="Path to the GFF3 file.")
    parser.add_argument("--source", nargs="+", help="Preferred annotation sources.")
    parser.add_argument("-o", "--output", default="feature_data.pkl", help="Output pickle file name (default: feature_data.pkl).")
    args = parser.parse_args()

    feature_df = parse_gff(args.gff_file, args.source)

    if feature_df is None:
        exit(1)

    if not feature_df.empty:
        try:
            with open(args.output, 'wb') as f:
                pickle.dump(feature_df, f)
            print(f"Feature data saved to {args.output}")
        except Exception as e:
            print(f"Error saving pickle file: {e}")
    else:
        print("No genes found in the GFF file (or none matching the specified source).")
