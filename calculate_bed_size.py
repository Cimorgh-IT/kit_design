import argparse

def calculate_bed_size_in_kb(bed_file):
    total_size = 0

    with open(bed_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            start = int(fields[1])
            end = int(fields[2])
            size = end - start
            total_size += size

    size_in_kb = total_size / 1000
    return size_in_kb

def main():
    parser = argparse.ArgumentParser(description='Calculate the total size of BED file regions in kilobases (kb).')
    parser.add_argument('bed_file', type=str, help='Path to the input BED file')
    
    args = parser.parse_args()
    
    size_in_kb = calculate_bed_size_in_kb(args.bed_file)
    print(f"Total size of regions in {args.bed_file}: {size_in_kb:.2f} kb")

if __name__ == '__main__':
    main()
