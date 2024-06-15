import pandas as pd
import argparse

# Function to load the BED file
def load_bed_file(bed_file_path):
    bed_df = pd.read_csv(bed_file_path, sep='\t', header=None)
    return bed_df

# Function to load the name conversion file
def load_conversion_file(conversion_file_path):
    conversion_df = pd.read_csv(conversion_file_path, sep='\t', header=None)
    return conversion_df

# Function to process the BED file using the conversion data
def process_bed_file(bed_df, conversion_dict):
    def replace_values(value):
        values = value.split(',')
        replaced_values = [conversion_dict.get(v, v) for v in values]
        return ','.join(replaced_values)

    for col in bed_df.columns:
        bed_df[col] = bed_df[col].astype(str).apply(replace_values)
    
    return bed_df

# Main function to run the script
def main():
    parser = argparse.ArgumentParser(description="Process a BED file using a name conversion file.")
    parser.add_argument("bed_file", help="Path to the BED file")
    parser.add_argument("conversion_file", help="Path to the name conversion file")
    parser.add_argument("-o", "--output", default="processed_bed_file.bed", help="Path to save the processed BED file (default: processed_bed_file.bed)")
    
    args = parser.parse_args()
    
    bed_df = load_bed_file(args.bed_file)
    conversion_df = load_conversion_file(args.conversion_file)
    
    conversion_dict = dict(zip(conversion_df[0], conversion_df[5]))
    
    processed_bed_df = process_bed_file(bed_df, conversion_dict)
    
    processed_bed_df.to_csv(args.output, sep='\t', header=False, index=False)
    print(f"Processed BED file saved to {args.output}")

if __name__ == "__main__":
    main()
