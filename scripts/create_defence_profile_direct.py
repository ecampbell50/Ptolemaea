#!/usr/bin/env python3
"""
create_defence_profile_direct.py
Create final defence profile directly from raw tool outputs with voting logic

PURPOSE:
Takes raw tool outputs and applies voting logic to create
a final defence profile with maximum interpretability and traceability.

Usage:
    python3 create_defence_profile_direct.py \
        --padloc 1004153.3_padloc.csv \
        --defensefinder 1004153.3_defense_finder_genes.tsv \
        --forward-blast 1004153.3_vs_bcereus_forward_best.txt \
        --reverse-blast 1004153.3_vs_bcereus_reverse_best.txt \
        --master-key MASTER_ToolKey.tsv \
        --output 1004153.3_defenceprofile.csv
"""

import pandas as pd
import argparse
import sys
import re
from pathlib import Path

def load_master_key_mappings(master_key_file):
    """
    Load master key file and create mapping dictionaries.
    Returns: (padloc_mapping, defensefinder_mapping, classification_lookup)
    """
    print(f"Loading master key: {master_key_file}")

    try:
        master_df = pd.read_csv(master_key_file, sep='\t')
        print(f"   Master key entries: {len(master_df)}")

        # Create mapping dictionaries
        # PADLOC: PADLOC_systems -> Novel_subtypes
        padloc_mapping = {}
        for _, row in master_df.iterrows():
            padloc_name = row['PADLOC_systems']
            consensus_name = row['Novel_subtypes']
            if pd.notna(padloc_name) and padloc_name != '/':
                padloc_mapping[padloc_name] = consensus_name

        # DefenseFinder: DefenseFinder_subtypes -> Novel_subtypes
        defensefinder_mapping = {}
        for _, row in master_df.iterrows():
            df_name = row['DefenseFinder_subtypes']
            consensus_name = row['Novel_subtypes']
            if pd.notna(df_name) and df_name != '/':
                defensefinder_mapping[df_name] = consensus_name

        # Classification lookup: Novel_subtypes -> (Novel_types, Defense_outcome)
        classification_lookup = {}
        for _, row in master_df.iterrows():
            subtype = row['Novel_subtypes']
            system_type = row['Novel_types']
            outcome = row['Defense_outcome']

            if pd.notna(subtype) and subtype != '/':
                classification_lookup[subtype] = {
                    'type': system_type if pd.notna(system_type) else 'UNMAPPED_TYPE',
                    'outcome': outcome if pd.notna(outcome) else 'UNMAPPED_OUTCOME'
                }

        print(f"   PADLOC mappings: {len(padloc_mapping)}")
        print(f"   DefenseFinder mappings: {len(defensefinder_mapping)}")
        print(f"   Classification entries: {len(classification_lookup)}")

        return padloc_mapping, defensefinder_mapping, classification_lookup

    except Exception as e:
        print(f"ERROR loading master key file: {e}")
        sys.exit(1)

def clean_trailing_underscore_one(name):
    """
    Remove trailing '_1' from defense system names, except for specific cases.
    """
    # COMMENT OUT THE NEXT LINE TO DISABLE THIS CLEANING FUNCTION
    # return name

    # Define exceptions - systems that should legitimately have '_1'
    exceptions = [
        'DISARM_1',
        'PD-T7-5_1',
        'GAO_19'  # This has '_19' not '_1', so won't be affected anyway
    ]

    # If the name is in our exceptions list, don't clean it
    if name in exceptions:
        return name

    # If it ends with '_1' and it's not an exception, remove it
    if name.endswith('_1'):
        return name[:-2]  # Remove last 2 characters ('_1')

    return name

def clean_defense_system_name(name):
    """
    Remove trailing '_1' from defense system names that appear in BLAST results.
    """
    if name.endswith('_1'):
        return name[:-2]  # Remove last 2 characters ('_1')
    return name

def extract_defense_system_from_blast_id(blast_id):
    """
    Extract defense system name from BLAST subject/query ID.
    Expected format: locustag#defensesystem_subtype_1
    """
    if '#' in blast_id:
        # Split by '#' and take the second part (defense system)
        defense_part = blast_id.split('#')[1]
        # Remove the trailing '_1' if present
        return clean_defense_system_name(defense_part)
    return blast_id

def extract_defense_name_from_blast(blast_result):
    """
    Extract defense system name from BLAST result string.
    """
    if blast_result == "No_hit":
        return None

    match = re.match(r'^([^(]+)', blast_result)
    if match:
        return match.group(1).strip()
    return blast_result

def extract_blast_metrics(blast_result):
    """
    Extract metrics from forward BLAST result.
    Returns: dict with pident, evalue, length, qlen, slen
    """
    if blast_result is None or blast_result == "No_hit":
        return None

    # Pattern: "name(pident%, E=evalue, L=length, Q=qlen, S=slen)"
    pattern = r'([^(]+)\(([0-9.]+)%, E=([0-9.e+-]+), L=([0-9]+), Q=([0-9]+), S=([0-9]+)\)'
    match = re.match(pattern, blast_result)

    if match:
        name, pident, evalue, length, qlen, slen = match.groups()
        return {
            'name': name.strip(),
            'pident': float(pident),
            'evalue': float(evalue),
            'length': int(length),
            'qlen': int(qlen),
            'slen': int(slen)
        }
    return None

def passes_blast_filtering(forward_metrics):
    """
    Check if BLAST-only hit passes filtering criteria.
    """
    if not forward_metrics:
        return False, "No metrics"

    Q = forward_metrics['qlen']
    S = forward_metrics['slen']
    L = forward_metrics['length']

    # Q/S ratio check
    qs_ratio = Q / S if S > 0 else 0
    if not (0.8 <= qs_ratio <= 1.25):
        return False, f"Q/S ratio {qs_ratio:.3f} outside 0.8-1.25"

    # Alignment coverage check
    avg_length = (Q + S) / 2
    coverage = L / avg_length if avg_length > 0 else 0
    if not (0.8 <= coverage <= 1.25):
        return False, f"Coverage {coverage:.3f} outside 0.8-1.25"

    return True, "Passed filtering"

def lookup_final_classification(final_consensus, classification_lookup):
    """
    Look up final system type and outcome from master key.
    Handles complex consensus names with special formatting.
    """
    if not classification_lookup:
        return 'UNMAPPED_TYPE', 'UNMAPPED_OUTCOME'

    # Handle special formatting cases
    if final_consensus.startswith('(p::') and '|d::' in final_consensus:
        # Extract individual names from MAPPING or CONFLICT cases
        pattern = r'\(p::([^|]*)\|d::([^)]*)\)'
        match = re.match(pattern, final_consensus)

        if match:
            padloc_name = match.group(1).strip() if match.group(1).strip() else None
            df_name = match.group(2).strip() if match.group(2).strip() else None

            # Look up each name
            padloc_type = padloc_outcome = None
            df_type = df_outcome = None

            if padloc_name and padloc_name in classification_lookup:
                padloc_type = classification_lookup[padloc_name]['type']
                padloc_outcome = classification_lookup[padloc_name]['outcome']

            if df_name and df_name in classification_lookup:
                df_type = classification_lookup[df_name]['type']
                df_outcome = classification_lookup[df_name]['outcome']

            # Format the composite result
            type_parts = []
            outcome_parts = []

            if padloc_type:
                type_parts.append(f"p::{padloc_type}")
            else:
                type_parts.append("p::UNMAPPED")

            if df_type:
                type_parts.append(f"d::{df_type}")
            else:
                type_parts.append("d::UNMAPPED")

            if padloc_outcome:
                outcome_parts.append(f"p::{padloc_outcome}")
            else:
                outcome_parts.append("p::UNMAPPED")

            if df_outcome:
                outcome_parts.append(f"d::{df_outcome}")
            else:
                outcome_parts.append("d::UNMAPPED")

            composite_type = f"({type_parts[0]}|{type_parts[1]})"
            composite_outcome = f"({outcome_parts[0]}|{outcome_parts[1]})"

            return composite_type, composite_outcome

        return 'UNMAPPED_TYPE', 'UNMAPPED_OUTCOME'

    # Simple consensus name - direct lookup
    if final_consensus in classification_lookup:
        result = classification_lookup[final_consensus]
        return result['type'], result['outcome']

    # No mapping found
    return 'UNMAPPED_TYPE', 'UNMAPPED_OUTCOME'

def process_raw_tool_outputs(padloc_file, defensefinder_file, forward_blast_file, reverse_blast_file, padloc_mapping, defensefinder_mapping):
    """
    Process raw tool outputs and create protein results dictionary.
    Enhanced with comprehensive validation and error handling.
    """
    print("\n" + "="*60)
    print("PROCESSING RAW TOOL OUTPUTS")
    print("="*60)

    protein_results = {}

    # Process PADLOC results
    print("\n1. Processing PADLOC results...")
    if not padloc_file.exists():
        print("   PADLOC file does not exist")
    elif padloc_file.stat().st_size == 0:
        print("   PADLOC file is empty (0 bytes) - no defence systems found")
    else:
        try:
            padloc_df = pd.read_csv(padloc_file)

            if len(padloc_df) == 0:
                print("   PADLOC file has headers but no data rows - no defence systems found")
            elif 'target.name' not in padloc_df.columns or 'system' not in padloc_df.columns:
                print("   Warning: PADLOC file missing required columns")
            else:
                print(f"   PADLOC entries found: {len(padloc_df)}")
                processed_count = 0

                for _, row in padloc_df.iterrows():
                    protein_id = row['target.name']
                    system = row['system']

                    # Validate protein ID
                    if pd.isna(protein_id) or not isinstance(protein_id, str) or protein_id.strip() == '':
                        print(f"   Warning: Skipping PADLOC row with invalid protein_id: {protein_id}")
                        continue

                    # Validate system
                    if pd.isna(system) or (isinstance(system, str) and system.strip() == ''):
                        print(f"   Warning: Skipping PADLOC row with invalid system for {protein_id}")
                        continue

                    # Clean trailing '_1' if inappropriate
                    system_cleaned = clean_trailing_underscore_one(str(system))

                    # Get consensus mapping
                    consensus_name = padloc_mapping.get(system_cleaned, 'No_mapping')

                    # Initialize protein entry
                    if protein_id not in protein_results:
                        protein_results[protein_id] = {
                            'padloc_orig': None, 'padloc_mapped': None,
                            'df_orig': None, 'df_mapped': None,
                            'forward_blast': None, 'reverse_blast': None
                        }

                    protein_results[protein_id]['padloc_orig'] = system_cleaned
                    protein_results[protein_id]['padloc_mapped'] = consensus_name
                    processed_count += 1

                print(f"   Processed {processed_count} valid PADLOC hits")

        except pd.errors.EmptyDataError:
            print("   PADLOC file is empty (no data) - no defence systems found")
        except pd.errors.ParserError as e:
            print(f"   Error parsing PADLOC file - possibly corrupted: {e}")
        except Exception as e:
            print(f"   Error processing PADLOC file: {e}")

    # Process DefenseFinder results
    print("\n2. Processing DefenseFinder results...")
    if not defensefinder_file.exists():
        print("   DefenseFinder file does not exist")
    elif defensefinder_file.stat().st_size == 0:
        print("   DefenseFinder file is empty (0 bytes) - no defence systems found")
    else:
        try:
            defensefinder_df = pd.read_csv(defensefinder_file, sep='\t')

            if len(defensefinder_df) == 0:
                print("   DefenseFinder file has headers but no data rows - no defence systems found")
            elif 'hit_id' not in defensefinder_df.columns or 'subtype' not in defensefinder_df.columns:
                print("   Warning: DefenseFinder file missing required columns")
            else:
                print(f"   DefenseFinder entries found: {len(defensefinder_df)}")
                processed_count = 0

                for _, row in defensefinder_df.iterrows():
                    protein_id = row['hit_id']
                    subtype = row['subtype']

                    # Validate protein ID
                    if pd.isna(protein_id) or not isinstance(protein_id, str) or protein_id.strip() == '':
                        print(f"   Warning: Skipping DefenseFinder row with invalid protein_id: {protein_id}")
                        continue

                    # Validate subtype
                    if pd.isna(subtype) or (isinstance(subtype, str) and subtype.strip() == ''):
                        print(f"   Warning: Skipping DefenseFinder row with invalid subtype for {protein_id}")
                        continue

                    # Clean trailing '_1' if inappropriate
                    subtype_cleaned = clean_trailing_underscore_one(str(subtype))

                    # Get consensus mapping
                    consensus_name = defensefinder_mapping.get(subtype_cleaned, 'No_mapping')

                    # Initialize protein entry
                    if protein_id not in protein_results:
                        protein_results[protein_id] = {
                            'padloc_orig': None, 'padloc_mapped': None,
                            'df_orig': None, 'df_mapped': None,
                            'forward_blast': None, 'reverse_blast': None
                        }

                    protein_results[protein_id]['df_orig'] = subtype_cleaned
                    protein_results[protein_id]['df_mapped'] = consensus_name
                    processed_count += 1

                print(f"   Processed {processed_count} valid DefenseFinder hits")

        except pd.errors.EmptyDataError:
            print("   DefenseFinder file is empty (no data) - no defence systems found")
        except pd.errors.ParserError as e:
            print(f"   Error parsing DefenseFinder file - possibly corrupted: {e}")
        except Exception as e:
            print(f"   Error processing DefenseFinder file: {e}")

    # Process Forward BLAST results
    print("\n3. Processing Forward BLAST results...")
    if not forward_blast_file.exists():
        print("   Forward BLAST file does not exist")
    elif forward_blast_file.stat().st_size == 0:
        print("   Forward BLAST file is empty - no matches found")
    else:
        try:
            blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
                            'qcovs', 'qlen', 'slen']

            forward_blast_df = pd.read_csv(forward_blast_file, sep='\t', header=None, names=blast_columns)

            if len(forward_blast_df) == 0:
                print("   Forward BLAST file has no data - no matches found")
            else:
                print(f"   Forward BLAST entries found: {len(forward_blast_df)}")
                processed_count = 0

                for _, row in forward_blast_df.iterrows():
                    protein_id = row['qseqid']
                    subject_id = row['sseqid']

                    # Validate protein ID
                    if pd.isna(protein_id) or not isinstance(protein_id, str) or protein_id.strip() == '':
                        print(f"   Warning: Skipping Forward BLAST row with invalid protein_id: {protein_id}")
                        continue

                    # Validate required numeric fields
                    try:
                        pident = float(row['pident'])
                        evalue = float(row['evalue'])
                        length = int(row['length'])
                        qlen = int(row['qlen'])
                        slen = int(row['slen'])
                    except (ValueError, TypeError):
                        print(f"   Warning: Skipping Forward BLAST row with invalid numeric data for {protein_id}")
                        continue

                    # Extract defense system name
                    cleaned_subject_id = extract_defense_system_from_blast_id(str(subject_id))

                    # Create comprehensive summary string
                    forward_result = f"{cleaned_subject_id}({pident:.1f}%, E={evalue:.1e}, L={length}, Q={qlen}, S={slen})"

                    # Initialize protein entry
                    if protein_id not in protein_results:
                        protein_results[protein_id] = {
                            'padloc_orig': None, 'padloc_mapped': None,
                            'df_orig': None, 'df_mapped': None,
                            'forward_blast': None, 'reverse_blast': None
                        }

                    protein_results[protein_id]['forward_blast'] = forward_result
                    processed_count += 1

                print(f"   Processed {processed_count} valid Forward BLAST hits")

        except Exception as e:
            print(f"   Error processing Forward BLAST file: {e}")

    # Process Reverse BLAST results
    print("\n4. Processing Reverse BLAST results...")
    if not reverse_blast_file.exists():
        print("   Reverse BLAST file does not exist")
    elif reverse_blast_file.stat().st_size == 0:
        print("   Reverse BLAST file is empty - no matches found")
    else:
        try:
            blast_columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                            'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
                            'qcovs', 'qlen', 'slen']

            reverse_blast_df = pd.read_csv(reverse_blast_file, sep='\t', header=None, names=blast_columns)

            if len(reverse_blast_df) == 0:
                print("   Reverse BLAST file has no data - no matches found")
            else:
                print(f"   Reverse BLAST entries found: {len(reverse_blast_df)}")
                processed_count = 0

                for _, row in reverse_blast_df.iterrows():
                    protein_id = row['sseqid']  # Subject is the genome protein in reverse BLAST
                    query_id = row['qseqid']

                    # Validate protein ID
                    if pd.isna(protein_id) or not isinstance(protein_id, str) or protein_id.strip() == '':
                        print(f"   Warning: Skipping Reverse BLAST row with invalid protein_id: {protein_id}")
                        continue

                    # Validate required numeric fields
                    try:
                        pident = float(row['pident'])
                        evalue = float(row['evalue'])
                    except (ValueError, TypeError):
                        print(f"   Warning: Skipping Reverse BLAST row with invalid numeric data for {protein_id}")
                        continue

                    # Extract defense system name
                    cleaned_query_id = extract_defense_system_from_blast_id(str(query_id))

                    # Create summary string
                    reverse_result = f"{cleaned_query_id}({pident:.1f}%, E={evalue:.1e})"

                    # Initialize protein entry
                    if protein_id not in protein_results:
                        protein_results[protein_id] = {
                            'padloc_orig': None, 'padloc_mapped': None,
                            'df_orig': None, 'df_mapped': None,
                            'forward_blast': None, 'reverse_blast': None
                        }

                    protein_results[protein_id]['reverse_blast'] = reverse_result
                    processed_count += 1

                print(f"   Processed {processed_count} valid Reverse BLAST hits")

        except Exception as e:
            print(f"   Error processing Reverse BLAST file: {e}")

    return protein_results

def determine_consensus(protein_data):
    """
    Apply voting logic to determine final consensus for a protein.
    """
    padloc_orig = protein_data['padloc_orig']
    padloc_mapped = protein_data['padloc_mapped']
    df_orig = protein_data['df_orig']
    df_mapped = protein_data['df_mapped']
    forward_blast = protein_data['forward_blast']
    reverse_blast = protein_data['reverse_blast']

    # Extract BLAST names
    forward_name = extract_defense_name_from_blast(forward_blast) if forward_blast else None
    reverse_name = extract_defense_name_from_blast(reverse_blast) if reverse_blast else None

    has_padloc = padloc_orig is not None
    has_df = df_orig is not None
    has_forward = forward_name is not None
    has_reverse = reverse_name is not None

    padloc_has_mapping = padloc_mapped not in [None, "No_mapping"]
    df_has_mapping = df_mapped not in [None, "No_mapping"]

    # Case 1: BLAST-only hits
    if not has_padloc and not has_df:
        if has_forward and has_reverse:
            if forward_name == reverse_name:
                # Check BLAST filtering criteria
                forward_metrics = extract_blast_metrics(forward_blast)
                if forward_metrics:
                    passed, reason = passes_blast_filtering(forward_metrics)
                    if passed:
                        return forward_name, "BLAST", f"BLAST-only hit passed filtering: {reason}"
                    else:
                        return None, "FILTERED", f"BLAST-only hit failed filtering: {reason}"
                else:
                    return None, "FILTERED", "Could not parse forward BLAST metrics"
            else:
                return None, "FILTERED", f"BLAST names disagree: {forward_name} vs {reverse_name}"
        else:
            return None, "FILTERED", "Insufficient BLAST evidence (need both forward and reverse)"

    # Case 2: Single tool hits
    elif has_padloc and not has_df:
        if padloc_has_mapping:
            return padloc_mapped, "SINGLE", "PADLOC only with mapping"
        else:
            return f"(p::{padloc_orig}|d::)", "MAPPING", "PADLOC only without mapping"

    elif not has_padloc and has_df:
        if df_has_mapping:
            return df_mapped, "SINGLE", "DefenseFinder only with mapping"
        else:
            return f"(p::|d::{df_orig})", "MAPPING", "DefenseFinder only without mapping"

    # Case 3: Both tools have hits
    elif has_padloc and has_df:
        # Check if we can do meaningful voting
        if padloc_has_mapping or df_has_mapping:
            # At least one tool has mapping - we can potentially resolve with BLAST
            padloc_consensus = padloc_mapped if padloc_has_mapping else None
            df_consensus = df_mapped if df_has_mapping else None

            # If both have mappings and agree
            if padloc_has_mapping and df_has_mapping and padloc_mapped == df_mapped:
                return padloc_mapped, "AGREE", "Both tools agree on consensus"

            # Use BLAST voting to resolve
            padloc_votes = 1  # Tool itself
            df_votes = 1      # Tool itself

            blast_evidence = []

            if has_forward:
                if padloc_consensus and forward_name == padloc_consensus:
                    padloc_votes += 1
                    blast_evidence.append(f"Forward supports PADLOC ({forward_name})")
                elif df_consensus and forward_name == df_consensus:
                    df_votes += 1
                    blast_evidence.append(f"Forward supports DefenseFinder ({forward_name})")
                else:
                    blast_evidence.append(f"Forward supports neither ({forward_name})")

            if has_reverse:
                if padloc_consensus and reverse_name == padloc_consensus:
                    padloc_votes += 1
                    blast_evidence.append(f"Reverse supports PADLOC ({reverse_name})")
                elif df_consensus and reverse_name == df_consensus:
                    df_votes += 1
                    blast_evidence.append(f"Reverse supports DefenseFinder ({reverse_name})")
                else:
                    blast_evidence.append(f"Reverse supports neither ({reverse_name})")

            # Determine winner
            if padloc_votes > df_votes:
                if padloc_has_mapping:
                    return padloc_mapped, "RESOLVED", f"PADLOC wins voting {padloc_votes}vs{df_votes}: {'; '.join(blast_evidence)}"
                else:
                    # PADLOC wins but has no mapping - use original name
                    return f"(p::{padloc_orig}|d::{df_orig})", "MAPPING", f"PADLOC wins voting but no mapping: {'; '.join(blast_evidence)}"
            elif df_votes > padloc_votes:
                if df_has_mapping:
                    return df_mapped, "RESOLVED", f"DefenseFinder wins voting {df_votes}vs{padloc_votes}: {'; '.join(blast_evidence)}"
                else:
                    # DefenseFinder wins but has no mapping - use original name
                    return f"(p::{padloc_orig}|d::{df_orig})", "MAPPING", f"DefenseFinder wins voting but no mapping: {'; '.join(blast_evidence)}"
            else:
                # Tied votes - show what we can
                if padloc_has_mapping and df_has_mapping:
                    return f"(p::{padloc_mapped}|d::{df_mapped})", "CONFLICT", f"Tied votes {padloc_votes}vs{df_votes}: {'; '.join(blast_evidence)}"
                else:
                    return f"(p::{padloc_orig}|d::{df_orig})", "MAPPING", f"Tied votes with mapping issues: {'; '.join(blast_evidence)}"

        else:
            # Neither tool has mapping
            return f"(p::{padloc_orig}|d::{df_orig})", "MAPPING", "Both tools without mapping"

    return None, "ERROR", "Unexpected case in consensus logic"

def main():
    parser = argparse.ArgumentParser(
        description="Create defence profile directly from raw tool outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python3 create_defence_profile_direct.py \\
        --padloc 1004153.3_padloc.csv \\
        --defensefinder 1004153.3_defense_finder_genes.tsv \\
        --forward-blast 1004153.3_vs_bcereus_forward_best.txt \\
        --reverse-blast 1004153.3_vs_bcereus_reverse_best.txt \\
        --master-key MASTER_ToolKey.tsv \\
        --output 1004153.3_defenceprofile.csv
        """
    )

    parser.add_argument('--padloc', required=True, help='PADLOC CSV output file')
    parser.add_argument('--defensefinder', required=True, help='DefenseFinder TSV output file')
    parser.add_argument('--forward-blast', required=True, help='Forward BLAST output file')
    parser.add_argument('--reverse-blast', required=True, help='Reverse BLAST output file')
    parser.add_argument('--master-key', required=True, help='Master tool key TSV file')
    parser.add_argument('--output', help='Output defence profile CSV file (auto-generated if not specified)')

    args = parser.parse_args()

    # Convert to Path objects
    padloc_file = Path(args.padloc)
    defensefinder_file = Path(args.defensefinder)
    forward_blast_file = Path(args.forward_blast)
    reverse_blast_file = Path(args.reverse_blast)
    master_key_file = Path(args.master_key)

    # Auto-generate output filename if not provided
    if args.output:
        output_file = Path(args.output)
    else:
        # Extract genome ID from padloc file name
        genome_id = padloc_file.stem.replace('_padloc', '')
        output_file = Path(f"{genome_id}_defenceprofile.csv")

    print(f"=== Creating Defence Profile (Direct from Raw Outputs) ===")
    print(f"Genome ID: {output_file.stem.replace('_defenceprofile', '')}")
    print(f"Output file: {output_file}")

    # Check all files exist
    required_files = [
        (padloc_file, "PADLOC"),
        (defensefinder_file, "DefenseFinder"),
        (forward_blast_file, "Forward BLAST"),
        (reverse_blast_file, "Reverse BLAST"),
        (master_key_file, "Master Tool Key")
    ]

    for file_path, file_type in required_files:
        if not file_path.exists():
            print(f"ERROR: {file_type} file not found: {file_path}")
            sys.exit(1)
        print(f"Found {file_type}: {file_path}")

    # Load master key mappings
    padloc_mapping, defensefinder_mapping, classification_lookup = load_master_key_mappings(master_key_file)

    # Process raw tool outputs
    protein_results = process_raw_tool_outputs(
        padloc_file, defensefinder_file, forward_blast_file, reverse_blast_file,
        padloc_mapping, defensefinder_mapping
    )

    print(f"\n" + "="*60)
    print("DETERMINING FINAL CONSENSUS")
    print("="*60)

    # Process each protein and determine consensus
    final_results = []
    stats = {
        'AGREE': 0, 'RESOLVED': 0, 'CONFLICT': 0, 'MAPPING': 0,
        'SINGLE': 0, 'BLAST': 0, 'FILTERED': 0, 'ERROR': 0
    }

    for protein_id in sorted([k for k in protein_results.keys() if isinstance(k, str) and pd.notna(k) and k.strip() != '']):
        protein_data = protein_results[protein_id]
        final_name, status, explanation = determine_consensus(protein_data)

        if final_name is not None:
            # Look up final classifications
            final_type, final_outcome = lookup_final_classification(final_name, classification_lookup)

            final_results.append({
                'protein_id': protein_id,
                'padloc_original': protein_data['padloc_orig'] or 'No_hit',
                'padloc_final': protein_data['padloc_mapped'] if protein_data['padloc_mapped'] not in [None, 'No_mapping'] else 'No_hit',
                'deffind_original': protein_data['df_orig'] or 'No_hit',
                'deffind_final': protein_data['df_mapped'] if protein_data['df_mapped'] not in [None, 'No_mapping'] else 'No_hit',
                'fwd_blast': extract_defense_name_from_blast(protein_data['forward_blast']) if protein_data['forward_blast'] is not None else 'No_hit',
                'rev_blast': extract_defense_name_from_blast(protein_data['reverse_blast']) if protein_data['reverse_blast'] is not None else 'No_hit',
                'status': status,
                'final_consensus': final_name,
                'explanation': explanation,
                'final_system_type': final_type,
                'final_system_subtype': final_name,
                'final_system_outcome': final_outcome
            })

        stats[status] += 1

        # Print progress for complex cases
        if status in ['RESOLVED', 'CONFLICT', 'BLAST']:
            print(f"  {protein_id}: {status} -> {final_name}")

    # Create and save output
    profile_df = pd.DataFrame(final_results)
    profile_df.to_csv(output_file, index=False)

    print(f"\n{'='*80}")
    print("DEFENCE PROFILE SUMMARY:")
    print(f"Total proteins processed: {len(protein_results)}")
    print(f"Proteins in final profile: {len(profile_df)}")
    print(f"Proteins filtered out: {stats['FILTERED']}")
    print()
    print("STATUS BREAKDOWN:")
    for status, count in stats.items():
        if count > 0:
            print(f"  {status}: {count} ({count/len(protein_results)*100:.1f}%)")

    # Count final classifications
    if len(profile_df) > 0:
        mapped_count = len(profile_df[profile_df['final_system_type'] != 'UNMAPPED_TYPE'])
        unmapped_count = len(profile_df[profile_df['final_system_type'] == 'UNMAPPED_TYPE'])
        print(f"\nFINAL CLASSIFICATIONS:")
        print(f"  Mapped to master key: {mapped_count} ({mapped_count/len(profile_df)*100:.1f}%)")
        print(f"  Unmapped (novel/species-specific): {unmapped_count} ({unmapped_count/len(profile_df)*100:.1f}%)")

        if mapped_count > 0:
            mapped_df = profile_df[profile_df['final_system_type'] != 'UNMAPPED_TYPE']
            outcome_counts = mapped_df['final_system_outcome'].value_counts()
            print(f"\nDEFENCE OUTCOME BREAKDOWN (mapped entries):")
            for outcome, count in outcome_counts.items():
                print(f"  {outcome}: {count}")

    print(f"\nFinal defence profile saved: {output_file}")

    # Show sample results
    print(f"\nSAMPLE RESULTS (first 5 entries):")
    if len(profile_df) > 0:
        sample_cols = ['protein_id', 'status', 'final_consensus', 'padloc_original', 'deffind_original']
        print(profile_df[sample_cols].head().to_string(index=False))

if __name__ == "__main__":
    main()
