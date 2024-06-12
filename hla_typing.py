import os
import argparse
import subprocess
import pandas as pd

def read_hla_genes(genes_list_file):
    hla_genes = []
    with open(genes_list_file, 'r') as f:
        for line in f:
            gene_name = line.strip().split('\t')[0]
            hla_genes.append(gene_name)
    return hla_genes

def gtf_to_bed(input_file, exon_output_file, padding, hla_genes):
    with open(input_file, 'r') as f_in, open(exon_output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom = fields[0].replace('chr', '')
            feature = fields[2]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            gene_id, exon_number, gene_name, transcript_id = None, None, None, None
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
            if gene_name in hla_genes and gene_id and exon_number and '_' not in chrom and (transcript_id.startswith("NM_") or transcript_id.startswith("ENST")):
                padded_start = max(start - padding, 0)
                padded_end = end + padding
                bed_line = f"{chrom}\t{padded_start}\t{padded_end}\t{feature}\t{strand}\t{gene_name}\t{transcript_id}\t{exon_number}\n"
                f_out.write(bed_line)

def add_exon_ids_and_sort(exon_output_file, sorted_exon_output_file):
    subprocess.run(f"sort -k1,1V -k2,2n -k3,3n {exon_output_file} | uniq > {sorted_exon_output_file}", shell=True, check=True)

def generate_gene_level_bed(sorted_exon_output_file, gene_output_file):
    gene_dict = {}
    with open(sorted_exon_output_file, 'r') as f_in:
        for line in f_in:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            strand = fields[4]
            gene_name = fields[5]
            if gene_name not in gene_dict:
                gene_dict[gene_name] = [chrom, start, end, strand]
            else:
                gene_dict[gene_name][1] = min(gene_dict[gene_name][1], start)
                gene_dict[gene_name][2] = max(gene_dict[gene_name][2], end)
    with open(gene_output_file, 'w') as f_out:
        for gene_name, values in gene_dict.items():
            chrom, start, end, strand = values
            bed_line = f"{chrom}\t{start}\t{end}\t{gene_name}\t{strand}\n"
            f_out.write(bed_line)

def sort_gene_level_bed(gene_output_file, sorted_gene_output_file):
    subprocess.run(f"sort -k4,4V -k1,1V -k2,2n -k3,3n {gene_output_file} > {sorted_gene_output_file}", shell=True, check=True)
    os.remove(gene_output_file)

def generate_filtered_sorted_exon_bed(sorted_exon_output_file, filtered_exon_output_file):
    with open(sorted_exon_output_file, 'r') as f_in, open(filtered_exon_output_file, 'w') as f_out:
        for line in f_in:
            fields = line.strip().split('\t')
            feature = fields[3]
            transcript_id = fields[6]
            if feature == "exon" and (transcript_id.startswith("NM_") or transcript_id.startswith("ENST")):
                f_out.write(line)
    sorted_filtered_exon_output_file = filtered_exon_output_file.replace('.bed', '_sorted.bed')
    subprocess.run(f"sort -k1,1V -k2,2n -k3,3n {filtered_exon_output_file} > {sorted_filtered_exon_output_file}", shell=True, check=True)
    os.remove(filtered_exon_output_file)
    return sorted_filtered_exon_output_file

def main():
    parser = argparse.ArgumentParser(description='Convert Ensembl GTF to BED format with exon IDs and padding.')
    parser.add_argument('input_file', type=str, help='Path to the input Ensembl GTF file')
    parser.add_argument('output_prefix', type=str, help='Prefix for the output BED files')
    parser.add_argument('padding', type=int, help='Padding distance to add to exon start and end positions')
    parser.add_argument('genes_list_file', type=str, help='Path to the panel genes list file')
    args = parser.parse_args()

    exon_output_file = f"{args.output_prefix}_exon_level.bed"
    sorted_exon_output_file = f"{args.output_prefix}_sorted_exon_level.bed"
    gene_output_file = f"{args.output_prefix}_gene_level.bed"
    sorted_gene_output_file = f"{args.output_prefix}_sorted_gene_level.bed"
    filtered_exon_output_file = f"{args.output_prefix}_filtered_exon_level.bed"

    hla_genes = read_hla_genes(args.genes_list_file)

    gtf_to_bed(args.input_file, exon_output_file, args.padding, hla_genes)
    add_exon_ids_and_sort(exon_output_file, sorted_exon_output_file)
    generate_gene_level_bed(sorted_exon_output_file, gene_output_file)
    sort_gene_level_bed(gene_output_file, sorted_gene_output_file)
    sorted_filtered_exon_output_file = generate_filtered_sorted_exon_bed(sorted_exon_output_file, filtered_exon_output_file)

    df = pd.read_csv(sorted_filtered_exon_output_file, sep='\t', header=None)
    df_unique = df.drop_duplicates(subset=[0, 1, 2, 3, 4, 5])
    unique_exon_level_output = f"{args.output_prefix}_unique_exon_level.bed"
    df_unique.to_csv(unique_exon_level_output, sep='\t', index=False, header=False)

    # Sort the final unique exon-level BED file
    final_sorted_exon_output_file = f"{args.output_prefix}_final_sorted_exon_level.bed"
    subprocess.run(f"sort -k6,6V -k1,1V -k2,2n -k3,3n {unique_exon_level_output} > {final_sorted_exon_output_file}", shell=True, check=True)

    os.remove(exon_output_file)
    os.remove(sorted_exon_output_file)
    os.remove(sorted_filtered_exon_output_file)
    os.remove(unique_exon_level_output)

    print("Sorted exon-level BED file created successfully.")
    print("Sorted gene-level BED file created successfully.")
    print(f"Filtered, unique, and sorted exon-level BED file created successfully: {final_sorted_exon_output_file}")

if __name__ == '__main__':
    main()
