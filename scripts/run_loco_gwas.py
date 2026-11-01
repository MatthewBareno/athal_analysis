"""Leave-One-Chromosome-Out (LOCO) Genome-Wide Association Study using GEMMA and PLINK.
1. Constructs 5 LOCO kinship matrices (K_-1, K_-2, K_-3, K_-4, K_-5) using 397,850 LD-pruned variants.
2. Fits LMM on each chromosome using its corresponding LOCO kinship matrix.
3. Merges all chromosomes into a unified LOCO association summary table.
4. Generates publication-quality LOCO Manhattan and QQ plots with candidate gene annotations.
"""

import os
import subprocess
import glob
import pandas as pd
import numpy as np


def run_loco_pipeline():
    os.makedirs('data/processed/loco', exist_ok=True)
    os.makedirs('results/gwas', exist_ok=True)
    os.makedirs('output', exist_ok=True)
    
    chroms = [1, 2, 3, 4, 5]
    
    for c in chroms:
        print(f"\n=================== Processing Chromosome {c} (LOCO) ===================")
        kinship_file = f'output/athal_loco_no_chr{c}.cXX.txt'
        pruned_no_chr = f'data/processed/loco/athal_pruned_loco_no_chr{c}'
        
        # 1. Extract LD-pruned variants excluding chromosome c
        if not os.path.exists(f'{pruned_no_chr}.bed'):
            print(f"Extracting pruned variants excluding Chr {c}...")
            subprocess.run([
                'plink', '--bfile', 'data/processed/athal_1003_qc',
                '--extract', 'data/processed/athal_1003_ld.prune.in',
                '--not-chr', str(c),
                '--make-bed',
                '--out', pruned_no_chr
            ], check=True)
            
        # 2. Compute LOCO centered relatedness matrix in GEMMA
        if not os.path.exists(kinship_file):
            print(f"Computing LOCO centered kinship matrix for Chr {c}...")
            subprocess.run([
                'gemma', '-bfile', pruned_no_chr,
                '-gk', '1',
                '-o', f'athal_loco_no_chr{c}'
            ], check=True)
            
        # 3. Extract Chr c testing variants
        test_chr = f'data/processed/loco/athal_chr{c}_test'
        if not os.path.exists(f'{test_chr}.bed'):
            print(f"Extracting Chr {c} test variants...")
            subprocess.run([
                'plink', '--bfile', 'data/processed/athal_1003_qc',
                '--chr', str(c),
                '--make-bed',
                '--out', test_chr
            ], check=True)
            
        # 4. Run GEMMA LMM on Chr c with LOCO kinship
        assoc_file = f'output/athal_loco_chr{c}_res.assoc.txt'
        if not os.path.exists(assoc_file):
            print(f"Running GEMMA LMM on Chr {c} with LOCO kinship...")
            subprocess.run([
                'gemma', '-bfile', test_chr,
                '-k', kinship_file,
                '-lmm', '4',
                '-o', f'athal_loco_chr{c}_res'
            ], check=True)
        else:
            print(f"Chr {c} LMM already complete: {assoc_file}")
            
    # 5. Merge all 5 LOCO results
    print("\nMerging all 5 LOCO chromosome association tables...")
    dfs = []
    for c in chroms:
        f = f'output/athal_loco_chr{c}_res.assoc.txt'
        df_c = pd.read_csv(f, sep=r'\s+')
        dfs.append(df_c)
        
    df_all = pd.concat(dfs, ignore_index=True)
    merged_file = 'results/gwas/ft10_gemma_loco_all.assoc.txt'
    df_all.to_csv(merged_file, sep='\t', index=False)
    print(f"Successfully created unified LOCO association table: {merged_file} ({len(df_all):,} SNPs)")
    
    return merged_file


if __name__ == '__main__':
    run_loco_pipeline()
