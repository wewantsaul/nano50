#!/usr/bin/env python3

import os
import argparse
import gzip
import csv
from pathlib import Path
from Bio import SeqIO


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Calculate N50 for nanopore reads.')
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-f', '--file', type=str, help='Input FASTQ file (can be .gz compressed).')
    input_group.add_argument('-d', '--directory', type=str, help='Directory containing FASTQ files.')
    parser.add_argument('-o', '--output', type=str, default='n50_results.csv', 
                        help='Output CSV file (default: n50_results.csv).')
    return parser.parse_args()


def is_fastq_file(filename):
    """Check if a file is a FASTQ file based on its extension."""
    extensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
    return filename.lower().endswith(extensions)


def get_fastq_files(directory):
    """Get all FASTQ files in a directory."""
    fastq_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if is_fastq_file(file):
                fastq_files.append(os.path.join(root, file))
    return fastq_files


def calculate_n50(read_lengths):
    """Calculate N50 from a list of read lengths."""
    if not read_lengths:
        return 0
    
    # Sort read lengths in descending order
    sorted_lengths = sorted(read_lengths, reverse=True)
    
    # Calculate total sum of lengths
    total_length = sum(sorted_lengths)
    
    # Find N50
    current_sum = 0
    for length in sorted_lengths:
        current_sum += length
        if current_sum >= total_length / 2:
            return length
    
    return 0


def process_fastq(filename):
    """Process a FASTQ file and return all read lengths."""
    read_lengths = []
    
    # Determine if the file is gzipped
    is_gzipped = filename.endswith('.gz')
    
    # Open file accordingly
    if is_gzipped:
        with gzip.open(filename, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                read_lengths.append(len(record.seq))
    else:
        with open(filename, 'r') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                read_lengths.append(len(record.seq))
    
    return read_lengths


def main():
    args = parse_arguments()
    
    # Collect files to process
    files_to_process = []
    if args.file:
        files_to_process.append(args.file)
    else:  # args.directory
        files_to_process = get_fastq_files(args.directory)
    
    if not files_to_process:
        print("No FASTQ files found to process.")
        return
    
    # Calculate N50 for each file
    results = []
    for file in files_to_process:
        print(f"Processing {file}...")
        try:
            read_lengths = process_fastq(file)
            
            if read_lengths:
                n50 = calculate_n50(read_lengths)
                total_reads = len(read_lengths)
                total_bases = sum(read_lengths)
                mean_length = total_bases / total_reads if total_reads > 0 else 0
                max_length = max(read_lengths) if read_lengths else 0
                min_length = min(read_lengths) if read_lengths else 0
                
                results.append({
                    'file': os.path.basename(file),
                    'n50': n50,
                    'total_reads': total_reads,
                    'total_bases': total_bases,
                    'mean_length': mean_length,
                    'max_length': max_length,
                    'min_length': min_length
                })
            else:
                print(f"No reads found in {file}")
        except Exception as e:
            print(f"Error processing {file}: {e}")
    
    # Write results to CSV
    if results:
        with open(args.output, 'w', newline='') as csvfile:
            fieldnames = ['file', 'n50', 'total_reads', 'total_bases', 'mean_length', 'max_length', 'min_length']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for result in results:
                writer.writerow(result)
            
        print(f"Results written to {args.output}")
    else:
        print("No results to write.")


if __name__ == '__main__':
    main()
