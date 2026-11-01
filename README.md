# athal_analysis
Population and quantitative genomic analyses of _Arabidopsis thaliana_ (1001 Genomes).

## Project Overview
This project performs population genomics and Genome-Wide Association Studies (GWAS) on _Arabidopsis thaliana_ accessions from the 1001 Genomes Project using flowering time phenotypes.

- **Phenotype:** FT10 (days to flowering at 10 °C, AraPheno study 12).
- **Genotypes:** 1001 Genomes v3.1 joint VCF (`1001genomes_snp-short-indel_only_ACGTN.vcf.gz`, 19.2 GB, 1,135 accessions, TAIR10).
- **Target Analysis Set:** **1,003 accessions** with both FT10 phenotype data and 1001 Genomes VCF genotypes (159 accessions from the phenotype study do not exist in any 1001G release and are omitted).

---

## Directory Layout
```
athal_analysis/
├── data/
│   ├── accession_list.csv         # 1001G master list (id, name, country, lat/lon, ...)
│   ├── pheno/
│   │   └── phenotype.csv          # Raw FT10 / FT16 data from AraPheno
│   ├── processed/                 # Prepared datasets (large binaries gitignored)
│   │   ├── common_samples.txt     # List of the 1,003 overlapping accession IDs
│   │   ├── phenotype_clean.csv    # Deduplicated, merged phenotype & metadata
│   │   ├── ft10_pheno_plink.txt   # PLINK-formatted phenotype table (#FID IID FT10)
│   │   ├── athal_1003_qc.pgen     # PLINK 2 binary genotype file (1,694,941 SNPs)
│   │   ├── athal_1003_qc.pvar     # PLINK 2 variant information
│   │   ├── athal_1003_qc.psam     # PLINK 2 sample information
│   │   └── athal_1003_qc.{bed,bim,fam} # PLINK 1.9 binary fileset for GEMMA / ADMIXTURE
│   ├── ref/                       # Reference genome (TAIR10, gitignored)
│   └── vcf/                       # Downloaded joint VCF & MD5 checksum (gitignored)
├── docs/
│   └── log.md                     # Dated lab notes and changelog
├── envs/
│   └── environment.yml            # Conda / mamba environment specification
├── results/
│   ├─- accession_world_map.html   # Professional interactive world map of accessions
│   ├─- accession_world_map.png    # High-resolution (2200x1160) publication map figure
│   └─- popgen/                    # Population genetics figures and evaluation
│       ├─- pca_admixture_composite.png # Nature-style 3-panel figure (PCA, K-means, ADMIXTURE)
│       ├─- pca_pc1_pc2.html       # Interactive WebGL PCA scatter plot
│       ├─- pca_pc1_pc2.png        # Publication PCA plot (300 DPI)
│       ├─- kmeans_evaluation.png  # Silhouette score and Inertia elbow plot
│       ├─- admixture_barplot_k9.png # K=9 ADMIXTURE stacked barplot
│       └─- kmeans_concordance.csv # K-means cluster vs lineage cross-tabulation
└── scripts/
    ├── download_data.sh           # Downloads 1001G v3.1 VCF and verifies MD5
    ├── prepare_phenotypes.py      # Deduplicates & formats phenotypes for the 1,003 accessions
    ├── prep_genotypes.sh          # Subsets VCF, applies QC filters, and creates PLINK files
    ├── plot_accession_map.py      # Generates publication world maps (HTML & PNG)
    └── analyze_popgen.py          # LD pruning, PLINK2 PCA, K-means & ADMIXTURE analysis
```

---

## Environment Setup

### Option 1: Python Virtual Environment (`.venv`)
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy pandas scipy scikit-learn matplotlib plotly
```

### Option 2: Conda / Mamba
```bash
conda env create -f envs/environment.yml
conda activate athal
```

Tools required in PATH:
- `bcftools` (v1.24+)
- `plink2` (v2.0.0+)
- `plink` (v1.9+)

---

## Data Preparation Pipeline

1. **Download Joint VCF:**
   ```bash
   bash scripts/download_data.sh
   ```
2. **Clean & Prepare Phenotypes:**
   ```bash
   .venv/bin/python scripts/prepare_phenotypes.py
   ```
   - Deduplicates accession `8424` replicates (both measured at FT10 = 68.25).
   - Intersects with VCF samples to yield the 1,003 target accessions.
   - Merges geographic coordinates (latitude, longitude, elevation, country).
   - Outputs `data/processed/phenotype_clean.csv` and `data/processed/ft10_pheno_plink.txt`.

3. **Subset & Quality Control Genotypes:**
   ```bash
   bash scripts/prep_genotypes.sh
   ```
   - Subsets to the 1,003 accessions.
   - Filters to strictly biallelic ACGT SNPs (`--snps-only just-acgt --min-alleles 2 --max-alleles 2`).
   - Applies minimum 90% call rate per variant (`--geno 0.10`).
   - Filters by minor allele frequency (`--maf 0.01`).
   - Assigns standardized variant IDs: `chr:pos:ref:alt` (`@:#:$r:$a`).
   - Generates PLINK 2 (`athal_1003_qc.pgen`) and PLINK 1.9 (`athal_1003_qc.bed`) datasets containing **1,694,941 high-quality SNPs**.

---

## Completed Analysis Pipeline & Milestones

1. **[Completed] Data Preparation & QC:**
   - Biallelic SNP filtering, 90% genotyping call rate, MAF ≥ 0.01 → **1,694,941 high-quality SNPs**.
   - Harmonization and deduplication with AraPheno study 12 → **1,003 accessions** with complete FT10 and genotype data.
   - Interactive & publication world accession maps ([`results/accession_world_map.png`](results/accession_world_map.png)).
2. **[Completed] Population Genetics & Ancestry Modeling (`scripts/analyze_popgen.py`):**
   - LD pruning (`50 10 0.2` in PLINK 2) → **397,850 neutral SNPs**.
   - PCA (20 PCs): PC1 (11.67%), PC2 (9.24%), capturing Iberian/Relict and North-South gradients.
   - Unsupervised K-means clustering ($K=9$, Silhouette score = 0.5261).
   - Global ADMIXTURE ancestry modeling ($K=9$ stacked barplot).
   - Pharmacogenetics-style Multi-PC AUC diagnostic matrix ([`results/popgen/pca_auc_matrix.png`](results/popgen/pca_auc_matrix.png)).
3. **[Completed] GWAS & Fine Mapping Methodology (`scripts/analyze_gwas.py`, `scripts/plot_clear_loco_manhattan.py`):**
   - **Model A (GLM + 10 PCs)**: Discovered *FT*, *DOG1*, *VIN3*, but exhibited severe residual inflation ($\lambda_{GC} = 2.44$).
   - **Model B (GEMMA Whole-Genome LMM)**: Controlled inflation ($\lambda_{GC} = 1.09$) via kinship matrix $\mathbf{K}$, but suffered proximal contamination.
   - **Model C (GEMMA LOCO LMM)**: Solved proximal contamination using leave-one-chromosome-out kinship matrices ($\lambda_{GC} = 1.108$).
   - **Model D (Conditioned Stepwise MLMM)**: Unmasked ***FLC* (*FLOWERING LOCUS C*)** at $P = 1.12 \times 10^{-22}$ and a secondary peak on Chr 2 ($P = 6.49 \times 10^{-18}$) after conditioning on lead QTLs.
   - Aligned publication Manhattan & QQ plots ([`results/gwas/ft10_loco_manhattan_clear.png`](results/gwas/ft10_loco_manhattan_clear.png)).
4. **[Completed] Genomic Prediction & Polygenic Scoring (`scripts/compute_polygenic_scores.py`):**
   - Narrow-sense heritability estimation: $h^2 = 88.06\%$.
   - **4-QTL Model (*FT*, *DOG1*, *VIN3*, *FLC*)**: 5-fold CV $R^2 = 42.93\%$, Pearson $r = 0.655$, MAE $= 10.72$ days.
   - **Whole-Genome GBLUP**: 5-fold CV $R^2 = \mathbf{66.64\%}$, Pearson $r = \mathbf{0.8173}$, MAE $= \mathbf{7.72}$ days.
   - Full individual scoring table for all 1,003 accessions ([`data/processed/athal_1003_polygenic_scores.csv`](data/processed/athal_1003_polygenic_scores.csv)).
   - Publication genomic prediction figure ([`results/prs/polygenic_score_evaluation.png`](results/prs/polygenic_score_evaluation.png)).
5. **[Completed] Macro-Regional ANOVA & Contrasts (`scripts/anova_prs_populations.py`):**
   - Fisher's One-Way ANOVA: $F(9, 993) = 136.11$, $P = 1.55 \times 10^{-166}$, $\eta^2 = 55.23\%$.
   - Welch's heteroscedastic ANOVA: $A = 638.03$, $P = 1.43 \times 10^{-131}$.
   - 45 pairwise contrasts: 34 significant under Tukey HSD, 38 significant under Benjamini-Hochberg FDR ($q < 0.05$).
   - Publication contrasts figure ([`results/prs/pgs_population_contrasts.png`](results/prs/pgs_population_contrasts.png)).
6. **[Completed] Evolutionary Tests of Divergent Selection (`scripts/test_divergent_selection.py`):**
   - **$Q_{ST} - F_{ST}$ Test**: $Q_{ST} = 0.5786$ vs weighted neutral $F_{ST} = 0.1210$ ($4.8\times$ excess variance ratio, $P < 10^{-20}$), proving selection over drift.
   - **Berg & Coop (2014) $Q_X$ Polygenic Test**: Overall $Q_X = 296.07$ ($df = 9$), $P_{\text{drift}} < 10^{-50}$ (7 / 9 kinship PCs show significant adaptive divergence).
   - **Causal QTL $F_{ST}$ Outlier Scan**: *FLC* peak locus has $F_{ST} = 0.7578$ (**top 0.04% of genome**, empirical $P = 0.0004$), *VIN3* top 2.8%, *DOG1* top 3.6%.
   - **Latitudinal Adaptive Gradient**: $+1.50$ days / °N latitude ($r = 0.650, P < 10^{-100}$).
   - Publication composite figure ([`results/prs/divergent_selection_tests.png`](results/prs/divergent_selection_tests.png)).

---

## References & Data Sources
- **1001 Genomes Consortium** (2016). *1,135 Genomes Reveal the Global Pattern of Polymorphism in Arabidopsis thaliana*. Cell, 166(2), 481-491.
- **AraPheno:** FT10 flowering time study 12 ([AraPheno citation](https://arapheno.1001genomes.org/faq/cite/)).
