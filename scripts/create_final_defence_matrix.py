#!/usr/bin/env python3
"""
create_final_defence_matrix.py
Apply manual curation and aggregate per-genome defence profiles into final outputs.

PURPOSE:
This is the final post-processing step of the Ptolemaea pipeline. It:
  1. Reads every per-genome consensus profile (`*_defenceprofile.csv`) produced by
     create_defence_profile_direct.py.
  2. Optionally applies the curated resolutions produced by editing the template
     from extract_unresolved_patterns.py, substituting user-supplied TYPE/SUBTYPE/
     OUTCOME values back into every MAPPING and CONFLICT protein.
  3. Writes three outputs across all genomes:
       <prefix>_annotations.csv  - tidy/long table, one row per defence gene
       <prefix>_matrix.csv       - wide genome x system matrix (gene counts);
                                   columns are <type>#<subtype>#<outcome>
       <prefix>_summary.tsv      - per-genome summary statistics

NOTE ON UNITS:
The Ptolemaea consensus profile is gene-centric: one row per predicted protein that
received a defence annotation. The matrix therefore reports the NUMBER OF DEFENCE
GENES per (genome, system), not the number of collapsed multi-gene system instances.
Use --binary to convert to simple presence/absence (0/1).

Usage:
    python3 create_final_defence_matrix.py \
        --consensus-dir output/05_consensus/ \
        --resolutions unresolved_patterns_CURATED.csv \
        --output-prefix species_defence

    # Without curation (MAPPING/CONFLICT genes become *_unresolved):
    python3 create_final_defence_matrix.py \
        --consensus-dir output/05_consensus/ \
        --output-prefix species_defence
"""

import pandas as pd
import argparse
import sys
import re
from pathlib import Path

# Statuses that were resolved automatically by the consensus script.
AUTO_STATUSES = {'AGREE', 'RESOLVED', 'SINGLE', 'BLAST'}
# Statuses that require manual curation.
CURATION_STATUSES = {'MAPPING', 'CONFLICT'}

# Sentinel values for defence genes that could not be resolved.
TYPE_UNRESOLVED = 'type_unresolved'
SUBTYPE_UNRESOLVED = 'subtype_unresolved'
OUTCOME_UNRESOLVED = 'outcome_unresolved'


def extract_system_name_from_blast(blast_str):
    """
    Extract just the system name from a BLAST hit string.
    Format: System(X.X%, E=Y.Ye-ZZ, ...) -> System
    Idempotent: a bare name with no parenthesis is returned unchanged.
    Mirrors the function of the same name in extract_unresolved_patterns.py so that
    pattern keys built here match the curated template exactly.
    """
    if blast_str is None or (isinstance(blast_str, float) and pd.isna(blast_str)):
        return 'No_hit'
    blast_str = str(blast_str)
    if blast_str == '' or blast_str == 'No_hit':
        return 'No_hit'
    match = re.match(r'^([^(]+)', blast_str)
    if match:
        return match.group(1).strip()
    return blast_str


def normalise_cell(value, is_blast=False):
    """Normalise a profile cell to match extract_unresolved_patterns.py pattern keys."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return 'No_hit'
    value = str(value).strip()
    if value == '' or value.lower() == 'nan':
        return 'No_hit'
    if is_blast:
        return extract_system_name_from_blast(value)
    return value


def build_pattern_key(row):
    """
    Reconstruct the (PADLOC, DefenseFinder, BLAST_fwd, BLAST_rev) pattern key for a
    profile row, using the SAME normalisation as extract_unresolved_patterns.py.
    """
    padloc = normalise_cell(row.get('padloc_original', 'No_hit'))
    deffind = normalise_cell(row.get('deffind_original', 'No_hit'))
    fwd = normalise_cell(row.get('fwd_blast', 'No_hit'), is_blast=True)
    rev = normalise_cell(row.get('rev_blast', 'No_hit'), is_blast=True)
    return (padloc, deffind, fwd, rev)


def genome_id_from_filename(file_path):
    """Derive the genome ID from a profile filename: <genome_id>_defenceprofile.csv."""
    return file_path.name.replace('_defenceprofile.csv', '')


def genome_id_from_protein(protein_id):
    """Derive the genome ID from a protein ID: <genome_id>@<locus_tag>."""
    if isinstance(protein_id, str) and '@' in protein_id:
        return protein_id.split('@')[0]
    return None


def load_resolutions(resolutions_file):
    """
    Load the curated resolutions CSV into a dict:
        (PADLOC, DefenseFinder, BLAST_fwd, BLAST_rev) -> (TYPE, SUBTYPE, OUTCOME)
    Only patterns with at least one non-blank curated field are kept.
    """
    print(f"Loading curated resolutions: {resolutions_file}")
    res_df = pd.read_csv(resolutions_file)

    required = ['PADLOC', 'DefenseFinder', 'BLAST_fwd', 'BLAST_rev', 'TYPE', 'SUBTYPE', 'OUTCOME']
    missing = [c for c in required if c not in res_df.columns]
    if missing:
        print(f"ERROR: resolutions file is missing required columns: {missing}")
        print(f"       Expected columns: {required}")
        sys.exit(1)

    def clean(v):
        if pd.isna(v):
            return ''
        v = str(v).strip()
        return '' if v.lower() == 'nan' else v

    resolutions = {}
    n_blank = 0
    for _, row in res_df.iterrows():
        key = (clean(row['PADLOC']), clean(row['DefenseFinder']),
               clean(row['BLAST_fwd']), clean(row['BLAST_rev']))
        type_v = clean(row['TYPE'])
        subtype_v = clean(row['SUBTYPE'])
        outcome_v = clean(row['OUTCOME'])

        # Skip patterns the curator has not touched at all.
        if not (type_v or subtype_v or outcome_v):
            n_blank += 1
            continue

        resolutions[key] = (
            type_v or TYPE_UNRESOLVED,
            subtype_v or SUBTYPE_UNRESOLVED,
            outcome_v or OUTCOME_UNRESOLVED,
        )

    print(f"   Curated patterns loaded: {len(resolutions)}")
    if n_blank:
        print(f"   WARNING: {n_blank} pattern(s) in the file were left entirely blank "
              f"and will be treated as unresolved.")
    return resolutions


def load_profiles(consensus_dir):
    """
    Load all per-genome profiles. Returns:
        genes_df   - concatenated DataFrame of all defence-gene rows (may be empty)
        all_genomes - sorted list of every genome ID (incl. those with zero defence genes)
    """
    profile_files = sorted(Path(consensus_dir).glob('*_defenceprofile.csv'))
    if not profile_files:
        print(f"ERROR: no '*_defenceprofile.csv' files found in {consensus_dir}")
        sys.exit(1)

    print(f"Found {len(profile_files)} per-genome profile files")

    all_genomes = []
    frames = []
    n_empty = 0

    for fp in profile_files:
        gid = genome_id_from_filename(fp)
        all_genomes.append(gid)

        # Empty profile (genome with no defence genes) -> 0 bytes or header-only.
        if fp.stat().st_size == 0:
            n_empty += 1
            continue
        try:
            df = pd.read_csv(fp)
        except pd.errors.EmptyDataError:
            n_empty += 1
            continue

        if df.empty or 'status' not in df.columns:
            n_empty += 1
            continue

        # Prefer genome ID from filename (robust even if protein_id is odd).
        df['genome_id'] = gid
        frames.append(df)

    all_genomes = sorted(set(all_genomes))

    if frames:
        genes_df = pd.concat(frames, ignore_index=True)
    else:
        genes_df = pd.DataFrame()

    print(f"   Genomes with defence genes: {len(all_genomes) - n_empty}")
    print(f"   Genomes with no defence genes: {n_empty}")
    print(f"   Total defence-gene rows: {len(genes_df)}")
    return genes_df, all_genomes


def apply_resolutions(genes_df, resolutions):
    """
    Apply curated resolutions and assign final TYPE/SUBTYPE/OUTCOME plus a
    resolution_source label to every gene row.
    """
    final_type = []
    final_subtype = []
    final_outcome = []
    source = []

    n_auto = n_curated = n_unresolved = 0

    for _, row in genes_df.iterrows():
        status = row.get('status')

        if status in AUTO_STATUSES:
            final_type.append(row.get('final_system_type'))
            final_subtype.append(row.get('final_system_subtype'))
            final_outcome.append(row.get('final_system_outcome'))
            source.append('auto')
            n_auto += 1

        elif status in CURATION_STATUSES:
            key = build_pattern_key(row)
            if key in resolutions:
                t, s, o = resolutions[key]
                final_type.append(t)
                final_subtype.append(s)
                final_outcome.append(o)
                source.append('curated')
                n_curated += 1
            else:
                final_type.append(TYPE_UNRESOLVED)
                final_subtype.append(SUBTYPE_UNRESOLVED)
                final_outcome.append(OUTCOME_UNRESOLVED)
                source.append('unresolved')
                n_unresolved += 1

        else:
            # Unexpected status; keep whatever is there, flag it.
            final_type.append(row.get('final_system_type'))
            final_subtype.append(row.get('final_system_subtype'))
            final_outcome.append(row.get('final_system_outcome'))
            source.append('auto')
            n_auto += 1

    genes_df = genes_df.copy()
    genes_df['final_type'] = final_type
    genes_df['final_subtype'] = final_subtype
    genes_df['final_outcome'] = final_outcome
    genes_df['resolution_source'] = source

    # Normalise stray whitespace so e.g. "Unknown " and "Unknown" do not become
    # separate matrix columns / summary categories (robust to older profiles).
    for col in ('final_type', 'final_subtype', 'final_outcome'):
        genes_df[col] = genes_df[col].apply(lambda v: str(v).strip() if pd.notna(v) else v)

    print("\nResolution summary:")
    print(f"   Automatic (AGREE/RESOLVED/SINGLE/BLAST): {n_auto}")
    print(f"   Curated (MAPPING/CONFLICT resolved):     {n_curated}")
    print(f"   Unresolved (no curation match):          {n_unresolved}")
    if n_unresolved:
        print(f"   NOTE: {n_unresolved} gene(s) remain *_unresolved. Re-run "
              f"extract_unresolved_patterns.py and complete the template to resolve them.")
    return genes_df


def build_tidy_table(genes_df):
    """Build the tidy/long annotation table (one row per defence gene)."""
    cols = ['genome_id', 'protein_id', 'status', 'resolution_source',
            'final_type', 'final_subtype', 'final_outcome',
            'padloc_original', 'deffind_original', 'fwd_blast', 'rev_blast',
            'final_consensus']
    present = [c for c in cols if c in genes_df.columns]
    tidy = genes_df[present].copy()
    tidy = tidy.sort_values(['genome_id', 'final_type', 'final_subtype', 'protein_id'],
                            kind='stable')
    return tidy


def build_matrix(genes_df, all_genomes, level='full', binary=False):
    """
    Build the wide genome x system matrix of gene counts (or 0/1 if binary).
    Every genome appears as a row, including those with zero defence genes.

    Column header format depends on --level:
      full    -> <type>#<subtype>#<outcome>   (default; preserves all three levels)
      subtype -> <subtype>
      type    -> <type>
    """
    if genes_df.empty:
        return pd.DataFrame(index=pd.Index(all_genomes, name='genome_id'))

    df = genes_df.copy()
    if level == 'subtype':
        df['_col'] = df['final_subtype'].astype(str)
    elif level == 'type':
        df['_col'] = df['final_type'].astype(str)
    else:  # 'full' — composite key, type#subtype#outcome
        df['_col'] = (df['final_type'].astype(str) + '#'
                      + df['final_subtype'].astype(str) + '#'
                      + df['final_outcome'].astype(str))

    counts = (df
              .groupby(['genome_id', '_col'])
              .size()
              .unstack(fill_value=0))

    # Ensure every genome is represented, even zero-defence ones.
    counts = counts.reindex(sorted(all_genomes), fill_value=0)
    counts.index.name = 'genome_id'
    counts = counts.reindex(sorted(counts.columns), axis=1)

    if binary:
        counts = (counts > 0).astype(int)
    return counts


def build_summary(genes_df, all_genomes):
    """Build a per-genome summary table."""
    rows = []
    if not genes_df.empty:
        grouped = dict(tuple(genes_df.groupby('genome_id')))
    else:
        grouped = {}

    for gid in sorted(all_genomes):
        g = grouped.get(gid)
        if g is None or g.empty:
            rows.append({
                'genome_id': gid, 'n_defence_genes': 0,
                'n_unique_subtypes': 0, 'n_unique_types': 0,
                'n_agree': 0, 'n_resolved': 0, 'n_single': 0, 'n_blast': 0,
                'n_curated': 0, 'n_unresolved': 0,
            })
            continue

        status_counts = g['status'].value_counts()
        src_counts = g['resolution_source'].value_counts()
        # Count "real" subtypes/types, excluding unresolved sentinels.
        real_subtypes = g.loc[g['final_subtype'] != SUBTYPE_UNRESOLVED, 'final_subtype']
        real_types = g.loc[~g['final_type'].isin([TYPE_UNRESOLVED, 'UNMAPPED_TYPE']), 'final_type']

        rows.append({
            'genome_id': gid,
            'n_defence_genes': len(g),
            'n_unique_subtypes': real_subtypes.nunique(),
            'n_unique_types': real_types.nunique(),
            'n_agree': int(status_counts.get('AGREE', 0)),
            'n_resolved': int(status_counts.get('RESOLVED', 0)),
            'n_single': int(status_counts.get('SINGLE', 0)),
            'n_blast': int(status_counts.get('BLAST', 0)),
            'n_curated': int(src_counts.get('curated', 0)),
            'n_unresolved': int(src_counts.get('unresolved', 0)),
        })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Apply curation and aggregate per-genome defence profiles into final outputs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Outputs (written next to --output-prefix):
    <prefix>_annotations.csv  tidy table, one row per defence gene
    <prefix>_matrix.csv       wide genome x system matrix (gene counts, or 0/1 with --binary);
                              columns are <type>#<subtype>#<outcome>
    <prefix>_summary.tsv      per-genome summary statistics
        """
    )
    parser.add_argument('--consensus-dir', required=True,
                        help='Directory containing *_defenceprofile.csv files')
    parser.add_argument('--resolutions', default=None,
                        help='Curated resolutions CSV from extract_unresolved_patterns.py '
                             '(optional; MAPPING/CONFLICT genes become *_unresolved if omitted)')
    parser.add_argument('--output-prefix', default='defence',
                        help='Prefix for output files (default: defence)')
    parser.add_argument('--level', choices=['full', 'subtype', 'type'], default='full',
                        help='Matrix column format: full=<type>#<subtype>#<outcome> (default), '
                             'subtype, or type')
    parser.add_argument('--binary', action='store_true',
                        help='Write 0/1 presence-absence instead of gene counts')

    args = parser.parse_args()

    consensus_dir = Path(args.consensus_dir)
    if not consensus_dir.exists():
        print(f"ERROR: directory not found: {consensus_dir}")
        sys.exit(1)

    print("=" * 70)
    print("CREATING FINAL DEFENCE MATRIX")
    print("=" * 70)

    # 1. Load resolutions (optional)
    resolutions = {}
    if args.resolutions:
        res_path = Path(args.resolutions)
        if not res_path.exists():
            print(f"ERROR: resolutions file not found: {res_path}")
            sys.exit(1)
        resolutions = load_resolutions(res_path)
    else:
        print("No --resolutions file supplied: MAPPING/CONFLICT genes will be marked "
              "*_unresolved.")

    # 2. Load profiles
    print()
    genes_df, all_genomes = load_profiles(consensus_dir)

    # 3. Apply resolutions / assign final classifications
    if not genes_df.empty:
        genes_df = apply_resolutions(genes_df, resolutions)
    else:
        print("\nNo defence genes found across any genome.")

    # 4. Build outputs
    out_annotations = Path(f"{args.output_prefix}_annotations.csv")
    out_matrix = Path(f"{args.output_prefix}_matrix.csv")
    out_summary = Path(f"{args.output_prefix}_summary.tsv")

    tidy = build_tidy_table(genes_df) if not genes_df.empty else pd.DataFrame()
    matrix = build_matrix(genes_df, all_genomes, level=args.level, binary=args.binary)
    summary = build_summary(genes_df, all_genomes)

    tidy.to_csv(out_annotations, index=False)
    matrix.to_csv(out_matrix)
    summary.to_csv(out_summary, sep='\t', index=False)

    # 5. Report
    print("\n" + "=" * 70)
    print("DONE")
    print("=" * 70)
    print(f"Genomes:            {len(all_genomes)}")
    print(f"Defence genes:      {len(tidy)}")
    print(f"Systems ({args.level}-level): {matrix.shape[1]} columns")
    value_kind = 'presence/absence (0/1)' if args.binary else 'gene counts'
    print(f"Matrix values:      {value_kind}")
    print(f"\nWrote:")
    print(f"  {out_annotations}  (tidy, one row per gene)")
    print(f"  {out_matrix}  (genome x {args.level} matrix)")
    print(f"  {out_summary}  (per-genome summary)")

    if not summary.empty and 'n_unresolved' in summary.columns:
        total_unresolved = int(summary['n_unresolved'].sum())
        if total_unresolved:
            print(f"\nWARNING: {total_unresolved} gene(s) remain unresolved across the panel.")


if __name__ == "__main__":
    main()
