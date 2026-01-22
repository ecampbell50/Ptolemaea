#!/usr/bin/env python3
"""
extract_unresolved_patterns.py
Extract unique patterns of unresolved defence gene annotations for manual curation

PURPOSE:
Identifies all problematic defence gene annotations (MAPPING and CONFLICT status)
and groups them by unique pattern for efficient manual curation.

Usage:
    python3 extract_unresolved_patterns.py \
        --consensus-dir 05_consensus/ \
        --output unresolved_patterns.csv
"""

import pandas as pd
import argparse
import sys
from pathlib import Path
from collections import defaultdict
import re

def extract_system_name_from_blast(blast_str):
    """
    Extract just the system name from BLAST hit string.
    Format: System(X.X%, E=Y.Ye-ZZ, ...) -> System

    Args:
        blast_str (str): BLAST hit string

    Returns:
        str: System name only
    """
    if not blast_str or blast_str == 'No_hit':
        return 'No_hit'

    # Extract everything before the first parenthesis
    match = re.match(r'^([^(]+)', blast_str)
    if match:
        return match.group(1).strip()

    return blast_str

def extract_genome_id(protein_id):
    """
    Extract genome ID from protein ID.
    Format: genome_id@locus_tag
    """
    if '@' in protein_id:
        return protein_id.split('@')[0]
    return protein_id

def load_consensus_files(consensus_dir):
    """
    Load all defence profile CSV files and extract problematic entries.

    Args:
        consensus_dir (Path): Directory containing defence profile files

    Returns:
        pd.DataFrame: All problematic entries
    """
    consensus_files = list(Path(consensus_dir).glob('*_defenceprofile.csv'))

    if not consensus_files:
        print(f"ERROR: No defence profile files found in {consensus_dir}")
        print(f"Looking for files matching: *_defenceprofile.csv")
        sys.exit(1)

    print(f"Found {len(consensus_files)} defence profile files")

    all_problematic = []

    for file_path in consensus_files:
        try:
            df = pd.read_csv(file_path)

            # Filter for problematic statuses
            problematic = df[df['status'].isin(['MAPPING', 'CONFLICT'])].copy()

            if not problematic.empty:
                # Add genome ID
                problematic['genome_id'] = problematic['protein_id'].apply(extract_genome_id)
                all_problematic.append(problematic)

        except Exception as e:
            print(f"WARNING: Error processing {file_path}: {e}")
            continue

    if not all_problematic:
        print("\nNo problematic genes found! All consensus files are clean.")
        return pd.DataFrame()

    return pd.concat(all_problematic, ignore_index=True)

def group_by_pattern(problematic_df):
    """
    Group problematic entries by unique annotation pattern.
    Uses ORIGINAL tool outputs to see what was actually detected.

    Args:
        problematic_df (pd.DataFrame): Problematic genes

    Returns:
        dict: Pattern -> list of protein IDs
    """
    patterns = defaultdict(list)

    for _, row in problematic_df.iterrows():
        # Use _original columns to see what tools actually found
        padloc = row.get('padloc_original', 'No_hit')
        if pd.isna(padloc) or padloc == '':
            padloc = 'No_hit'

        defensefinder = row.get('deffind_original', 'No_hit')
        if pd.isna(defensefinder) or defensefinder == '':
            defensefinder = 'No_hit'

        fwd_blast = row.get('fwd_blast', 'No_hit')
        if pd.isna(fwd_blast) or fwd_blast == '':
            fwd_blast = 'No_hit'
        else:
            fwd_blast = extract_system_name_from_blast(fwd_blast)

        rev_blast = row.get('rev_blast', 'No_hit')
        if pd.isna(rev_blast) or rev_blast == '':
            rev_blast = 'No_hit'
        else:
            rev_blast = extract_system_name_from_blast(rev_blast)

        # Create pattern tuple
        pattern = (padloc, defensefinder, fwd_blast, rev_blast)
        patterns[pattern].append(row['protein_id'])

    return patterns

def create_unresolved_patterns_csv(patterns):
    """
    Create DataFrame for manual curation.

    Args:
        patterns (dict): Pattern -> list of protein IDs

    Returns:
        pd.DataFrame: Unresolved patterns for curation
    """
    rows = []

    for pattern, protein_ids in sorted(patterns.items(), key=lambda x: len(x[1]), reverse=True):
        padloc, defensefinder, fwd_blast, rev_blast = pattern

        # Get example protein IDs (up to 5)
        example_proteins = ', '.join(protein_ids[:5])
        if len(protein_ids) > 5:
            example_proteins += f', ... ({len(protein_ids)-5} more)'

        rows.append({
            'PADLOC': padloc,
            'DefenseFinder': defensefinder,
            'BLAST_fwd': fwd_blast,
            'BLAST_rev': rev_blast,
            'protein_count': len(protein_ids),
            'example_proteins': example_proteins,
            'TYPE': '',
            'SUBTYPE': '',
            'OUTCOME': ''
        })

    return pd.DataFrame(rows)

def main():
    parser = argparse.ArgumentParser(
        description="Extract unique patterns of unresolved defence gene annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python3 extract_unresolved_patterns.py \\
        --consensus-dir 05_consensus/ \\
        --output unresolved_patterns.csv

After running:
    1. Open unresolved_patterns.csv in Excel
    2. Fill in TYPE, SUBTYPE, and OUTCOME columns for each pattern
    3. Use 'type_unresolved', 'subtype_unresolved', or 'outcome_unresolved'
       when you cannot determine a value with confidence
    4. Save the file (e.g., as unresolved_patterns_CURATED.csv)
    5. Run create_final_defence_matrix.py with your curated file
        """
    )

    parser.add_argument('--consensus-dir', required=True,
                       help='Directory containing defence profile CSV files (*_defenceprofile.csv)')
    parser.add_argument('--output', default='unresolved_patterns.csv',
                       help='Output CSV file for manual curation (default: unresolved_patterns.csv)')

    args = parser.parse_args()

    consensus_dir = Path(args.consensus_dir)
    output_file = Path(args.output)

    print("=" * 70)
    print("EXTRACTING UNRESOLVED DEFENCE GENE PATTERNS")
    print("=" * 70)
    print()

    # Check directory exists
    if not consensus_dir.exists():
        print(f"ERROR: Directory not found: {consensus_dir}")
        sys.exit(1)

    # Load all problematic entries
    print("Step 1: Loading consensus files...")
    problematic_df = load_consensus_files(consensus_dir)

    if problematic_df.empty:
        print("\n" + "=" * 70)
        print("SUCCESS: No problematic genes found!")
        print("All defence profile files are clean and ready for matrix creation.")
        print("=" * 70)
        sys.exit(0)

    print(f"  Total problematic genes: {len(problematic_df)}")
    print(f"  Status breakdown:")
    for status, count in problematic_df['status'].value_counts().items():
        print(f"    {status}: {count}")

    # Group by pattern
    print("\nStep 2: Grouping by unique annotation pattern...")
    patterns = group_by_pattern(problematic_df)
    print(f"  Unique patterns found: {len(patterns)}")

    # Create output DataFrame
    print("\nStep 3: Creating curation template...")
    unresolved_df = create_unresolved_patterns_csv(patterns)

    # Save to CSV
    unresolved_df.to_csv(output_file, index=False)

    print(f"  Saved to: {output_file}")

    # Print summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Total problematic proteins: {len(problematic_df)}")
    print(f"Unique patterns to review: {len(patterns)}")
    print(f"\nTop 10 most common patterns:")
    print("-" * 70)

    for i, (pattern, protein_ids) in enumerate(sorted(patterns.items(),
                                                       key=lambda x: len(x[1]),
                                                       reverse=True)[:10], 1):
        padloc, defensefinder, fwd_blast, rev_blast = pattern
        print(f"{i:2d}. n={len(protein_ids):4d}  PADLOC:{padloc:15s} DF:{defensefinder:15s} "
              f"Fwd:{fwd_blast:15s} Rev:{rev_blast:15s}")

    # Print instructions
    print("\n" + "=" * 70)
    print("NEXT STEPS")
    print("=" * 70)
    print(f"1. Open {output_file} in Excel or similar")
    print("2. Review each pattern and fill in TYPE, SUBTYPE, and OUTCOME columns")
    print()
    print("IMPORTANT GUIDELINES:")
    print("  • Use 'type_unresolved' when type cannot be determined with confidence")
    print("  • Use 'subtype_unresolved' when subtype cannot be determined")
    print("  • Use 'outcome_unresolved' when outcome cannot be determined")
    print("  • These are DIFFERENT from 'Unknown' values in the original data")
    print("  • 'Unknown' = tool found it but classified as unknown")
    print("  • '_unresolved' = manual curator could not resolve conflicting annotations")
    print()
    print("EXAMPLE RESOLUTIONS:")
    print("  Pattern: PADLOC=CBASS_other, DF=CBASS_IIs, BLAST_fwd=CBASS_II, BLAST_rev=CBASS_II")
    print("  → TYPE=CBASS, SUBTYPE=CBASS_IIs, OUTCOME=Abi")
    print("    (BLAST and DF agree on specific subtype)")
    print()
    print("  Pattern: PADLOC=Dynamins, DF=Eleos, BLAST_fwd=No_hit, BLAST_rev=No_hit")
    print("  → TYPE=Dynamins, SUBTYPE=type_unresolved, OUTCOME=outcome_unresolved")
    print("    (Tools disagree, no BLAST to help resolve)")
    print()
    print("3. Save your curated file (e.g., unresolved_patterns_CURATED.csv)")
    print("4. Run: python3 create_final_defence_matrix.py \\")
    print(f"         --consensus-dir {consensus_dir} \\")
    print("         --resolutions unresolved_patterns_CURATED.csv \\")
    print("         --output-prefix species_defence")
    print()
    print("=" * 70)

if __name__ == "__main__":
    main()
