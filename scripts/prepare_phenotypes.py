"""Prepare and clean phenotype data for the 1,003 Arabidopsis thaliana accessions."""

import os
import pandas as pd


def main():
    os.makedirs('data/processed', exist_ok=True)

    # 1. Load phenotype data
    pheno_path = 'data/pheno/phenotype.csv'
    df_pheno = pd.read_csv(pheno_path)
    print(f"Loaded raw phenotype data: {len(df_pheno)} rows")

    # 2. Handle duplicates
    # Replicate 7078 and 7079 for accession 8424 both have FT10 = 68.25
    # Average across replicates per accession_id
    df_dedup = df_pheno.groupby('accession_id', as_index=False).agg({
        'FT10': 'mean',
        'FT16': 'mean',
        'replicate_id': 'first',
    })
    print(f"Deduplicated to {len(df_dedup)} unique accessions")

    # 3. Load common accessions (intersection with 1001G VCF)
    common_path = 'data/processed/common_samples.txt'
    with open(common_path) as f:
        common_ids = set(int(line.strip()) for line in f if line.strip())
    print(f"Common samples with VCF: {len(common_ids)}")

    # 4. Filter to the 1,003 accessions
    df_clean = df_dedup[df_dedup['accession_id'].isin(common_ids)].copy()
    df_clean = df_clean.sort_values('accession_id')
    print(f"Cleaned phenotype set: {len(df_clean)} accessions")

    # 5. Merge geographic and accession metadata
    meta_path = 'data/accession_list.csv'
    if os.path.exists(meta_path):
        df_meta = pd.read_csv(meta_path)
        # accession_list id matches accession_id
        df_clean = df_clean.merge(
            df_meta[['id', 'name', 'country', 'sitename', 'latitude', 'longitude', 'elevation', 'csnumber']],
            left_on='accession_id',
            right_on='id',
            how='left'
        ).drop(columns=['id'])

    # Save cleaned CSV
    out_csv = 'data/processed/phenotype_clean.csv'
    df_clean.to_csv(out_csv, index=False)
    print(f"Saved cleaned phenotype table: {out_csv} ({len(df_clean)} rows)")

    # 6. Generate PLINK phenotype files
    # PLINK 2 format: #FID IID FT10 (or #IID FT10)
    plink_pheno = pd.DataFrame({
        '#FID': df_clean['accession_id'],
        'IID': df_clean['accession_id'],
        'FT10': df_clean['FT10']
    })
    plink_out = 'data/processed/ft10_pheno_plink.txt'
    plink_pheno.to_csv(plink_out, sep='\t', index=False)
    print(f"Saved PLINK phenotype file: {plink_out}")


if __name__ == '__main__':
    main()
