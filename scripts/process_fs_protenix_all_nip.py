import pandas as pd
from tqdm import tqdm
import os
import glob
import json
import statistics


def process_foldseek_directory(results_path, file_extension, filter_taxname=None):
    """Process foldseek results from a given directory"""
    foldseek_results_list = []
    
    df_cols = ["query","target","alntmscore","lddt","fident","alnlen",
        "mismatch","gapopen","qstart","qend","tstart","tend","evalue",
        "bits","taxid","taxname","taxlineage"]

    for filepath in tqdm(glob.glob(os.path.join(results_path, f"*{file_extension}"))):
        
        transcript_id = os.path.basename(filepath).replace(file_extension, "")
        results_dict = {}
        results_dict['idx'] = transcript_id

        try:
            df = pd.read_csv(filepath, sep='\t', header=None)
        except pd.errors.EmptyDataError:
            for col in df_cols:
                results_dict[col] = 'NA'
            foldseek_results_list.append(results_dict)
            continue

        df.columns = df_cols
        
        # Filter for evalue < 1E-5
        df = df[df['evalue'].astype(float) < 1e-5]
        
        if df.empty:
            for col in df_cols:
                results_dict[col] = 'NA'
            foldseek_results_list.append(results_dict)
            continue
        
        # Filter by taxname if specified (for AFSP results)
        if filter_taxname:
            df = df[df['taxname'] != filter_taxname]
            if df.empty:
                for col in df_cols:
                    results_dict[col] = 'NA'
                foldseek_results_list.append(results_dict)
                continue
        
        # Get top hit with lowest evalue
        best_hit_df = df.nsmallest(1, 'evalue')
        
        for col in df_cols:
            results_dict[col] = best_hit_df[col].values[0]

        foldseek_results_list.append(results_dict)

    return pd.DataFrame(foldseek_results_list)


# Define directories
pdb_dir = "./pdb"
afsp_dir = "./afsp"

# Process PDB results
print('Getting foldseek results from PDB directory')
pdb_results_df = process_foldseek_directory(pdb_dir, ".m8")
pdb_cols_to_rename = {col: f"{col}.pdb" for col in pdb_results_df.columns if col != 'idx'}
pdb_results_df = pdb_results_df.rename(columns=pdb_cols_to_rename)

# Process AFSP results
print('Getting foldseek results from AFSP directory')
afsp_results_df = process_foldseek_directory(afsp_dir, ".m8", filter_taxname="Oryza sativa Japonica Group")
afsp_cols_to_rename = {col: f"{col}.sp" for col in afsp_results_df.columns if col != 'idx'}
afsp_results_df = afsp_results_df.rename(columns=afsp_cols_to_rename)

# Merge dataframes
print('Merging results')
merged_results_df = pd.merge(pdb_results_df, afsp_results_df, on='idx', how='outer')

# Save results
merged_results_df.to_csv('all_protenix_fs_results_df.csv', index=False)
print(f'Saved merged results with {len(merged_results_df)} transcripts')
