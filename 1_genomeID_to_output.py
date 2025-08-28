#!/usr/bin/env python3
"""
Script to update DefenseFinder and PADLOC output files to include genome ID prefix in locus tags
Author: Emmet Campbell
Date: 2025-08-21

Usage:
    # Process a single file
    python3 1_genomeID_to_output.py 123456_padloc.csv
    python3 1_genomeID_to_output.py 123456_defense_finder_genes.tsv
    
    # Process entire directory structure (original functionality)
    python3 1_genomeID_to_output.py --directory AllGenomeOutputDirs
          where "AllGenomeOutputDirs" is a directory containing an output directory for each genome with the suff  "_TOOLOUTPUTS"
          AllGenomeOutputDirs/123456_TOOLOUTPUTS/123456_padloc.csv (and 123456_defense_finder_genes.tsv) 
"""

import os
import sys
import shutil
import pandas as pd
import argparse
from pathlib import Path

def extract_genome_id_from_filename(filename):
    """Extract genome ID from filename by removing known suffixes"""
    base_name = Path(filename).stem
    
    # Remove common suffixes to get genome ID
    suffixes_to_remove = ['_padloc', '_defense_finder_genes']
    
    for suffix in suffixes_to_remove:
        if base_name.endswith(suffix):
            return base_name.replace(suffix, '')
    
    # If no known suffix found, assume everything before first underscore is genome ID
    return base_name.split('_')[0]

def update_padloc_csv(csv_file, genome_id):
    """Update PADLOC CSV file to add genome ID prefix to target.name column"""
    print(f"Processing PADLOC file: {csv_file}")

    # Create backup
    backup_file = str(csv_file) + ".original"
    shutil.copy2(csv_file, backup_file)
    print(f"  Created backup: {backup_file}")

    # Read CSV
    df = pd.read_csv(csv_file)

    # Check if target.name column exists
    if 'target.name' not in df.columns:
        print(f"  WARNING: 'target.name' column not found in {csv_file}")
        return False

    # Count original entries
    original_count = len(df)

    # Update target.name column to add genome ID prefix
    df['target.name'] = genome_id + '@' + df['target.name'].astype(str)

    # Save updated file
    df.to_csv(csv_file, index=False)
    print(f"  Updated {original_count} entries in target.name column")
    return True

def update_defensefinder_tsv(tsv_file, genome_id):
    """Update DefenseFinder TSV file to add genome ID prefix to hit_id column"""
    print(f"Processing DefenseFinder file: {tsv_file}")

    # Create backup
    backup_file = str(tsv_file) + ".original"
    shutil.copy2(tsv_file, backup_file)
    print(f"  Created backup: {backup_file}")

    # Read TSV
    df = pd.read_csv(tsv_file, sep='\t')

    # Check if hit_id column exists
    if 'hit_id' not in df.columns:
        print(f"  WARNING: 'hit_id' column not found in {tsv_file}")
        return False

    # Count original entries
    original_count = len(df)

    # Update hit_id column to add genome ID prefix
    df['hit_id'] = genome_id + '@' + df['hit_id'].astype(str)

    # Save updated file
    df.to_csv(tsv_file, sep='\t', index=False)
    print(f"  Updated {original_count} entries in hit_id column")
    return True

def process_single_file(file_path):
    """Process a single PADLOC or DefenseFinder file"""
    file_path = Path(file_path)
    
    if not file_path.exists():
        print(f"ERROR: File {file_path} does not exist!")
        return False
    
    # Extract genome ID from filename
    genome_id = extract_genome_id_from_filename(file_path.name)
    print(f"Extracted genome ID: {genome_id}")
    print()
    
    success = False
    
    # Determine file type and process accordingly
    if file_path.suffix.lower() == '.csv' and 'padloc' in file_path.name.lower():
        success = update_padloc_csv(file_path, genome_id)
    elif file_path.suffix.lower() == '.tsv' and 'defense_finder_genes' in file_path.name.lower():
        success = update_defensefinder_tsv(file_path, genome_id)
    else:
        print(f"ERROR: File type not recognised. Expected:")
        print(f"  - *padloc*.csv for PADLOC files")
        print(f"  - *defense_finder_genes*.tsv for DefenseFinder files")
        return False
    
    if success:
        print(f"Successfully processed {file_path}")
    else:
        print(f"Failed to process {file_path}")
    
    return success

def process_directory(base_dir):
    """Process all defense output files in directory structure (original functionality)"""
    base_dir = Path(base_dir)
    
    if not base_dir.exists():
        print(f"ERROR: Directory {base_dir} does not exist!")
        return

    # Find all genome directories
    genome_dirs = [d for d in base_dir.iterdir() if d.is_dir() and d.name.endswith("_TOOLOUTPUTS")]

    print(f"Found {len(genome_dirs)} genome directories to process")
    print()

    padloc_success = 0
    defensefinder_success = 0
    padloc_total = 0
    defensefinder_total = 0

    for genome_dir in sorted(genome_dirs):
        # Extract genome ID from directory name (remove "_Lucy_DSs" suffix)
        genome_id = genome_dir.name.replace("_TOOLOUTPUTS", "")
        print(f"Processing genome: {genome_id}")

        # Process PADLOC CSV file
        padloc_file = genome_dir / f"{genome_id}_padloc.csv"
        if padloc_file.exists():
            padloc_total += 1
            if update_padloc_csv(padloc_file, genome_id):
                padloc_success += 1
        else:
            print(f"  PADLOC file not found: {padloc_file}")

        # Process DefenseFinder TSV file
        defensefinder_file = genome_dir / f"{genome_id}_defense_finder_genes.tsv"
        if defensefinder_file.exists():
            defensefinder_total += 1
            if update_defensefinder_tsv(defensefinder_file, genome_id):
                defensefinder_success += 1
        else:
            print(f"  DefenseFinder file not found: {defensefinder_file}")

        print()

    print("=== SUMMARY ===")
    print(f"PADLOC files: {padloc_success}/{padloc_total} successfully updated")
    print(f"DefenseFinder files: {defensefinder_success}/{defensefinder_total} successfully updated")
    print()
    print("All original files have been backed up with .original extension")

def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description="Update DefenseFinder and PADLOC output files to include genome ID prefix",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single PADLOC file
  python3 1_genomeID_to_output.py 123456_padloc.csv
  
  # Process a single DefenseFinder file  
  python3 1_genomeID_to_output.py 123456_defense_finder_genes.tsv
  
  # Process entire directory structure
  python3 1_genomeID_to_output.py --directory AllGenomeOutputDirs
        """
    )
    
    parser.add_argument('input', nargs='?', help='Input file to process')
    parser.add_argument('--directory', '-d', help='Process entire directory structure')
    
    args = parser.parse_args()
    
    # Check arguments
    if args.directory:
        # Directory mode (original functionality)
        process_directory(args.directory)
    elif args.input:
        # Single file mode
        process_single_file(args.input)
    else:
        # No arguments provided
        print("ERROR: Please provide either a file to process or use --directory flag")
        print()
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
