#!/usr/bin/env python3
"""
Script to add genome ID prefix to protein sequence headers in FAA files
Author: Emmet Campbell
Date: 2025-08-28

Usage:
    # Process a single FAA file
    python3 add_genomeID_to_faa.py 123456.faa
    
    # Process multiple FAA files in a directory
    python3 add_genomeID_to_faa.py --directory /path/to/faa_files/
"""

import os
import sys
import shutil
import argparse
from pathlib import Path

def extract_genome_id_from_filename(filename):
    """Extract genome ID from filename by removing .faa extension"""
    return Path(filename).stem

def process_faa_file(faa_file):
    """Process a single FAA file to add genome ID prefix to headers"""
    faa_path = Path(faa_file)
    
    if not faa_path.exists():
        print(f"ERROR: File {faa_path} does not exist!")
        return False
    
    if faa_path.suffix.lower() != '.faa':
        print(f"ERROR: File {faa_path} is not a .faa file!")
        return False
    
    # Extract genome ID from filename
    genome_id = extract_genome_id_from_filename(faa_path.name)
    print(f"Processing FAA file: {faa_path}")
    print(f"Genome ID: {genome_id}")
    
    # Create backup
    backup_file = str(faa_path) + ".original"
    shutil.copy2(faa_path, backup_file)
    print(f"Created backup: {backup_file}")
    
    # Read and process the file
    updated_headers = 0
    total_sequences = 0
    
    try:
        with open(faa_path, 'r') as infile:
            lines = infile.readlines()
        
        with open(faa_path, 'w') as outfile:
            for line in lines:
                if line.startswith('>'):
                    total_sequences += 1
                    # Remove the '>' and strip whitespace
                    header_content = line[1:].strip()
                    
                    # Check if genome ID is already present
                    if not header_content.startswith(f"{genome_id}@"):
                        # Add genome ID prefix
                        new_header = f">{genome_id}@{header_content}\n"
                        updated_headers += 1
                    else:
                        # Keep original if already has prefix
                        new_header = line
                    
                    outfile.write(new_header)
                else:
                    # Write sequence lines unchanged
                    outfile.write(line)
        
        print(f"Updated {updated_headers}/{total_sequences} sequence headers")
        
        if updated_headers == 0 and total_sequences > 0:
            print("  All headers already had genome ID prefix")
        
        return True
        
    except Exception as e:
        print(f"ERROR processing {faa_path}: {e}")
        # Restore backup if something went wrong
        shutil.copy2(backup_file, faa_path)
        return False

def process_directory(directory):
    """Process all FAA files in a directory"""
    dir_path = Path(directory)
    
    if not dir_path.exists():
        print(f"ERROR: Directory {dir_path} does not exist!")
        return
    
    if not dir_path.is_dir():
        print(f"ERROR: {dir_path} is not a directory!")
        return
    
    # Find all FAA files
    faa_files = list(dir_path.glob("*.faa"))
    
    if not faa_files:
        print(f"No .faa files found in {dir_path}")
        return
    
    print(f"Found {len(faa_files)} FAA files to process")
    print()
    
    success_count = 0
    
    for faa_file in sorted(faa_files):
        if process_faa_file(faa_file):
            success_count += 1
        print()
    
    print("=== SUMMARY ===")
    print(f"Successfully processed {success_count}/{len(faa_files)} FAA files")
    print("All original files have been backed up with .original extension")

def main():
    """Main function with command line argument parsing"""
    parser = argparse.ArgumentParser(
        description="Add genome ID prefix to protein sequence headers in FAA files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process a single FAA file
  python3 add_genomeID_to_faa.py 123456.faa
  
  # Process all FAA files in a directory
  python3 add_genomeID_to_faa.py --directory /path/to/faa_files/
  
  # Process FAA files in current directory
  python3 add_genomeID_to_faa.py --directory .

Header transformation:
  Before: >locus_001 hypothetical protein
  After:  >123456@locus_001 hypothetical protein
        """
    )
    
    parser.add_argument('input', nargs='?', help='Input FAA file to process')
    parser.add_argument('--directory', '-d', help='Process all FAA files in directory')
    
    args = parser.parse_args()
    
    # Check arguments
    if args.directory:
        # Directory mode
        process_directory(args.directory)
    elif args.input:
        # Single file mode
        process_faa_file(args.input)
    else:
        # No arguments provided
        print("ERROR: Please provide either a FAA file to process or use --directory flag")
        print()
        parser.print_help()
        sys.exit(1)

if __name__ == "__main__":
    main()
