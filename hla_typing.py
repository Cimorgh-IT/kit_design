import os
import argparse
import subprocess

# List of HLA gene names
hla_genes = ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", 
             "HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", 
             "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
             "HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPB2",
             "HLA-DQB2","HLA-DQA2","HLA-DRA","HLA-DRB6","HLA-F-AS1",
             "HLA-G","HLA-H","HLA-J","HLA-L","HLA-V"]

# Define function to convert Ensembl GTF to BED format with padding
def gtf_to_bed(input_file, exon_output_file, padding):
    with open(input_file, 'r') as f_in, open(exon_output_file, 'w') as f_out:
        for line in f_in:
            # Skip comments
            if line.startswith('#'):
                continue
            
            # Split line into fields
            fields = line.strip().split('\t')
            
            # Extract relevant fields
            chrom = fields[0].replace('chr', '')  # Remove 'chr' prefix if present
            feature = fields[2]
            start = int(fields[3]) - 1  # Convert to 0-based
            end = int(fields[4])
            strand = fields[6]
            
            # Extract gene_id, exon_number, and gene_name from the attributes field
            attributes = fields[8]
            gene_id = None
            exon_number = None
            gene_name = None
            transcript_id = None
            for attribute in attributes.split(';'):
                attribute = attribute.strip()
                if attribute.startswith('gene_id'):
                    gene_id = attribute.split(' ')[1].replace('"', '')
                elif attribute.startswith('exon_number'):
                    exon_number = attribute.split(' ')[1].replace('"', '')
                elif attribute.startswith('gene_name'):
                    gene_name = attribute.split(' ')[1].replace('"', '')
                elif attribute.startswith('transcript_id'):
                    transcript_id = attribute.split(' ')[1].replace('"','')
            
            # Check if the feature is an exon
            if gene_name in hla_genes and gene_id and exon_number and '_' not in chrom:
                # Adjust start and end positions by padding
                padded_start = max(start - padding, 0)  # Ensure start is not negative
                padded_end = end + padding
                
                # Write BED-formatted line with strand, gene_id, and exon_number
                bed_line = f"{chrom}\t{padded_start}\t{padded_end}\t{feature}\t{strand}\t{gene_name}\t{transcript_id}\t{exon_number}\n"
                f_out.write(bed_line)

def add_exon_ids_and_sort(exon_output_file, sorted_exon_output_file):
    # Sort the exon-level BED file using the sort command
    subprocess.run(f"sort -k6,6V -k1,1V -k2,2n -k3,3n {exon_output_file} | uniq > {sorted_exon_output_file}", shell=True, check=True)

def generate_gene_level_bed(sorted_exon_output_file, gene_output_file):
    gene_dict = {}

    # Read the sorted exon-level BED file
    with open(sorted_exon_output_file, 'r') as f_in:
        for line in f_in:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[4]
            gene_name = fields[5]

            # Update gene start and end positions
            if gene_name not in gene_dict:
                gene_dict[gene_name] = [chrom, start, end, strand]
            else:
                gene_dict[gene_name][1] = min(gene_dict[gene_name][1], start)
                gene_dict[gene_name][2] = max(gene_dict[gene_name][2], end)
    
    # Write the gene-level BED file
    with open(gene_output_file, 'w') as f_out:
        for gene_name, values in gene_dict.items():
            chrom, start, end, strand = values
            bed_line = f"{chrom}\t{start}\t{end}\t{gene_name}\t{strand}\n"
            f_out.write(bed_line)

def sort_gene_level_bed(gene_output_file, sorted_gene_output_file):
    # Sort the gene-level BED file using the sort command
    subprocess.run(f"sort -k4,4V -k1,1V -k2,2n -k3,3n {gene_output_file} > {sorted_gene_output_file}", shell=True, check=True)
    os.remove(gene_output_file)

def generate_filtered_sorted_exon_bed(sorted_exon_output_file, filtered_exon_output_file):
    with open(sorted_exon_output_file, 'r') as f_in, open(filtered_exon_output_file, 'w') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            feature = fields[3]
            transcript_id = fields[6]
            
            if feature == "CDS" and transcript_id.startswith("NM_"):
                f_out.write(line)
    
    # Sort the filtered exon-level BED file using the sort command
    sorted_filtered_exon_output_file = filtered_exon_output_file.replace('.bed', '_sorted.bed')
    subprocess.run(f"sort -k6,6V -k1,1V -k2,2n -k3,3n {filtered_exon_output_file} > {sorted_filtered_exon_output_file}", shell=True, check=True)
    os.remove(filtered_exon_output_file)
    return sorted_filtered_exon_output_file

def main():
    parser = argparse.ArgumentParser(description='Convert Ensembl GTF to BED format with exon IDs and padding.')
    parser.add_argument('input_file', type=str, help='Path to the input Ensembl GTF file')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output BED files')
    parser.add_argument('padding', type=int, help='Padding distance to add to exon start and end positions')
    
    args = parser.parse_args()

    exon_output_file = f"{args.output_prefix}_exon_level.bed"
    sorted_exon_output_file = f"{args.output_prefix}_sorted_exon_level.bed"
    gene_output_file = f"{args.output_prefix}_gene_level.bed"
    sorted_gene_output_file = f"{args.output_prefix}_sorted_gene_level.bed"
    filtered_exon_output_file = f"{args.output_prefix}_filtered_exon_level.bed"

    # Convert Ensembl GTF to exon-level BED format with padding
    gtf_to_bed(args.input_file, exon_output_file, args.padding)
    
    # Sort the exon-level BED file and generate sorted exon-level BED file
    add_exon_ids_and_sort(exon_output_file, sorted_exon_output_file)
    
    # Generate gene-level BED file
    generate_gene_level_bed(sorted_exon_output_file, gene_output_file)
    
    # Sort the gene-level BED file and generate sorted gene-level BED file
    sort_gene_level_bed(gene_output_file, sorted_gene_output_file)

    # Generate filtered and sorted exon-level BED file
    sorted_filtered_exon_output_file = generate_filtered_sorted_exon_bed(sorted_exon_output_file, filtered_exon_output_file)

    # Remove temporary files
    os.remove(exon_output_file)
    os.remove(sorted_exon_output_file)

    print("Sorted exon-level BED file created successfully.")
    print("Sorted gene-level BED file created successfully.")
    print(f"Filtered and sorted exon-level BED file created successfully: {sorted_filtered_exon_output_file}")

if __name__ == '__main__':
    main()
