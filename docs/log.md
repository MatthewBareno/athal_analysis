# Log

## Date: 2026-02-28
- Wrote the "goals" document
- Created the mamba environment
- main goal for tomorrow will be to get the (test subset of) read files, reference file, and phenotype files all formatted well'

## Date: 2026-03-01
- Downloaded reference genome from [ncbi](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/) using select file source as "All" and select file types as "FASTA" and "JSONL".
- moved to ref/ folder
- unzipped reference genome and deleted .zip
- extracted gcf .fna file and deleted everything else
- in the script i will ultimately create that auto performs this analysis, it will automatically rename the reference genome to ref.fna. It is the future user's responsibility to make sure to keep track of the information (ID and source) of the ref genome. But for now I renamed it to ref.fna

## Date: 2026-03-02
- familliarized myself with the [phenotpe website](https://arapheno.1001genomes.org/accessions/) to see what my options are in terms of maximizing the amount of accessions while getting a phenotype that is measured in as much accessions as possible
- be sure to back-reference [this page](https://arapheno.1001genomes.org/faq/cite/) for citation later
- used [this link](https://arapheno.1001genomes.org/phenotypes/?sort=-values&page=3) to determine "days to flowering trait" has the most amount of 1001GenomeProject accessions phenotyped plus is probably heritable enough to be good for GWAS.
- Therefore the phenotype I will be analyzing is: [FT10](https://arapheno.1001genomes.org/phenotype/261/)
- Phenotype data was accessible via [this page](https://arapheno.1001genomes.org/study/12/)

## Date: 2026-03-03
- <ins> phenotype.csv</ins> refers to the days to flower (DTF) at 16 and 10 degrees celsius. There are nmore accessions with 10 degrees DTF data, so I'll proceed with that.
- In order to mass download read files the IDs need to be connected to some sort of identifier that ncbi (which hosts the genome read files) recognizes. Therefore this phenotype file, which merely has "accession ID", needs some sort of reference file to link that accession ID with the ID its hosted on ncbi.
- If you go to the [main directory](https://1001genomes.org/) and scroll to the bottom under "Tables", you can see the [*Master List*](https://docs.google.com/spreadsheets/d/1Jd_TLRYLSi_MOL_Qbf3kaLXQUDoWVs7nl-BRRf59VO4/edit?usp=sharing). This has <ins>"id"</ins>, which accords with the accession ID in the phenotype file, <ins>name</ins>, and "<ins> CS Number</ins>". So now we need to somehow connect and/or use this information to automate downloads from blast for the read files.

## Date: 2026-03-04
- Downloaded the master accession list, and renamed it to `accession_list.csv` on data/

## Date: 2026-05-22
- Set up environment on my new computer

## Date: 2026-10-30
- FT10 has the most accessions (1,162). 1,003 are in the 1001G v3.1 joint VCF (and the ENA reads study PRJNA273563); the other 159 are in no 1001G release (checked current, 27genomes, legacy GMI).
- Raw reads for the 1,003 would be ~3.3 TB, so switched plan: use the pre-called VCF (19.2 GB) instead of aligning and calling. Download deferred until at school (`scripts/download_data.sh`).
- Simplified the repo: removed test reads, SRA manifest, per-folder READMEs, goals.md (now in README), broken venv.
- If more power is needed later: RegMap 250k SNP chip + imputation for the 159.

## Date: 2026-10-31
- Verified environment and tooling: Ubuntu 26.04 LTS, AMD Ryzen AI 5 340 (12 cores), 30 GB RAM, 47 GB available disk space. Pre-installed bcftools 1.24 and plink2 v2.0.0-a.7.8LM are in PATH.
- Created dedicated `.venv` virtual environment and installed core scientific stack: numpy, pandas, scipy, scikit-learn, matplotlib, plotly.
- Verified VCF download checksum and sample count: 1,135 accessions present in `data/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz`.
- Checked phenotype data (`data/pheno/phenotype.csv`): 1,162 unique accessions; resolved duplicate replicate for accession `8424` (both replicates have identical FT10 = 68.25 days). FT10 trait has 100% completeness.
- Intersected VCF samples with FT10 phenotypes to confirm exactly 1,003 common accessions saved to `data/processed/common_samples.txt`.
- Created `scripts/prepare_phenotypes.py`: generated clean merged dataset `data/processed/phenotype_clean.csv` (with geographic and site metadata) and PLINK-compatible phenotype file `data/processed/ft10_pheno_plink.txt`.
- Created `scripts/prep_genotypes.sh` to subset the 19.2 GB VCF to the 1,003 common accessions, filter for biallelic SNPs (`--snps-only just-acgt --min-alleles 2 --max-alleles 2`), filter by call rate (`--geno 0.10`), standardize variant IDs (`@:#:$r:$a`), and export both PLINK2 (PGEN) and PLINK1 (BED/BIM/FAM) binary formats.
- Executed genotype QC pipeline with threshold **MAF >= 0.01** (retaining low-frequency alleles for popgen and fine-mapping):
  - 12,883,854 total variants scanned from joint VCF across 1,135 samples.
  - Subsetting to 1,003 target accessions retained all 1,003 samples.
  - 2,176,424 non-biallelic / indel sites excluded.
  - 4,861,385 variants excluded with >10% missing call rate (`--geno 0.10`).
  - 4,151,104 variants excluded with MAF < 0.01 (`--maf 0.01`).
  - Final filtered dataset: **1,694,941 high-confidence biallelic SNPs** across 1,003 accessions.
  - Generated canonical datasets: `data/processed/athal_1003_qc.pgen` (189 MB) / `.pvar` (75 MB) / `.psam` (7.9 KB) and `data/processed/athal_1003_qc.bed` (406 MB) / `.bim` (51 MB) / `.fam` (16 KB).
- Deleted legacy mapping script (`scripts/map_samples.py`) and output (`results/sample_map.html`).
- Created `scripts/plot_accession_map.py` and generated publication-grade Natural Earth world maps:
  - Interactive web visual: `results/accession_world_map.html` (compact customdata, FT10 gradient, country/elevation tooltips).
  - High-resolution 2200x1160 figure: `results/accession_world_map.png`.

## 2026-10-31: Population Genetics & Population Structure (PCA, K-means, ADMIXTURE)

- **LD Pruning (`scripts/analyze_popgen.py` via PLINK 2)**:
  - Executed `--indep-pairwise 50 10 0.2` across 1,694,941 filtered SNPs.
  - Retained **397,850 independent, linkage-pruned biallelic SNPs** (`athal_1003_ld.prune.in`).
- **Principal Component Analysis (PLINK 2)**:
  - Calculated top 20 eigenvectors and eigenvalues (`athal_1003_pca.eigenvec`, `athal_1003_pca.eigenval`).
  - **PC1 (11.67% variance)**: Isolates the East/Central Asian lineage (+0.104) from the European metapopulation.
  - **PC2 (9.24% variance)**: Separates ancestral Iberian Relict accessions (-0.135) and Spain from Northern and Central Europe.
  - **PC3 (7.18% variance)**: Delineates North Sweden (+0.093) from Italy/Balkans/Caucasus (-0.043).
- **Unsupervised K-Means Clustering (`scikit-learn`)**:
  - Clustered accessions along top 10 PCs across $K = 2 \dots 12$.
  - Calculated Inertia (WCSS) and Silhouette Score per $K$.
  - Found a sharp global peak at **$K = 9$ (Silhouette Score = 0.5261)** with a steep elbow in inertia, objectively confirming the 9 evolutionary lineages identified in Alonso-Blanco et al. (*Cell* 2016).
  - Generated concordance matrix (`results/popgen/kmeans_concordance.csv`) proving tight correspondence with biological ancestry groups (e.g., 67/69 Asian accessions clustered into a single pure cluster; North Sweden 60/61 pure).
- **ADMIXTURE Analysis & Visualizations**:
  - Mapped official 1001 Genomes ADMIXTURE ancestry $Q$-matrix ($K = 9$) across all 1,003 accessions.
  - Generated 300-DPI publication figures in `results/popgen/`:
    - `pca_admixture_composite.png`: High-resolution 3-panel publication visual (A: PCA scatter, B: K-means Silhouette/Elbow evaluation, C: Global ADMIXTURE stacked barplot).
    - `pca_pc1_pc2.html` & `pca_pc1_pc2.png`: Interactive WebGL scatter plot and 300-DPI static plot.
    - `kmeans_evaluation.png`: Dual-axis plot of Silhouette Score and Inertia elbow.
    - `admixture_barplot_k9.png`: Standalone high-res ADMIXTURE barplot.
  - Merged popgen covariates into canonical table `data/processed/athal_1003_pca_kmeans.csv`.
- **Biogeographic Grouping System Figure (Pharmacogenetics style)**:
  - Created `scripts/plot_pca_auc_matrix.py` replicating the multi-panel layout from *Clin Pharmacol Ther (2019)*.
  - **Panel (a)**: 3x3 lower-triangular pairwise PC scatter matrix (PC1 vs PC2, PC1 vs PC3, PC2 vs PC3, PC1 vs PC4, PC2 vs PC4, PC3 vs PC4) colored by ancestry group.
  - **Panel (b)**: 5-fold cross-validated One-vs-Rest ROC-AUC across PCs 1 to 10 for each genetic ancestry group.
  - Saved outputs: `results/popgen/pca_auc_biogeographic_system.png` and `results/popgen/pca_ovr_auc_results.csv`.

## 2026-10-31: Genome-Wide Association Study (GWAS) Accounting for Population Structure

- **Phenotype & Covariates Setup**:
  - Target trait: FT10 (flowering time at 10 °C) across 1,003 accessions.
  - Synchronized `.fam` and `.psam` formats (`ft10_pheno_plink2.txt`, `ft10_pheno_gemma.txt`).
- **Model 1: PLINK 2 GLM with 10 PC Covariates**:
  - Tested 1,694,941 SNPs with `--glm allow-no-covars hide-covar`.
  - Top hits at ***DOG1*** ($P = 5.70 \times 10^{-29}$) and ***FLC*** ($P = 1.40 \times 10^{-27}$).
  - Genomic inflation factor: $\lambda_{GC} = 2.6738$, demonstrating that continental PCs alone leave residual cryptic relatedness uncorrected.
- **Model 2: GEMMA Linear Mixed Model (LMM with Kinship Matrix)**:
  - Calculated 1003x1003 centered relatedness matrix (`output/athal_1003_kinship.cXX.txt`).
  - Fit univariate LMM with REMLE variance estimation (`-lmm 4`).
  - **Narrow-sense Heritability ($h^2$ / PVE)**: **$88.06\% \pm 2.15\%$** ($V_g = 607.32, V_e = 20.41$).
  - **Genomic Inflation Controlled**: $\lambda_{GC} = 0.9794$, demonstrating near-perfect calibration.
  - **Genome-Wide Significant Peaks ($P < 7.31 \times 10^{-8}$)**:
    1. ***FT*** (*FLOWERING LOCUS T*, Chr 1: 24,339,560): $P = 1.46 \times 10^{-9}$ ($\beta = +2.64$ days). Master mobile florigen signal.
    2. ***DOG1*** (*DELAY OF GERMINATION 1*, Chr 5: 18,590,327): $P = 2.46 \times 10^{-9}$ ($\beta = +3.19$ days). Low-temperature flowering and seed dormancy.
    3. ***VIN3*** (*VERNALIZATION INSENSITIVE 3*, Chr 5: 23,288,908): $P = 8.25 \times 10^{-10}$ ($\beta = +8.42$ days). Epigenetic vernalization histone modifier.
- **Generated Assets (`results/gwas/`)**:
  - `gwas_glm_vs_lmm_comparison.png`: 4-panel figure comparing GLM (10 PCs) and LMM (Kinship Matrix) Manhattan & QQ plots.
  - `ft10_gemma_lmm_manhattan_qq.png`: Standalone publication Manhattan and QQ plot for the Linear Mixed Model.
  - `ft10_gemma_lmm_significant_snps.csv`: Lead significant SNPs with betas, standard errors, and allele frequencies.
  - `ft10_gemma_lmm.html` & `ft10_plink2_glm_pc10.html`: Interactive WebGL Manhattan plots.

## 2026-10-31: Advanced GWAS Methodology — LOCO Kinship & Stepwise Multi-Locus Model (MLMM)

- **Model 3: Leave-One-Chromosome-Out (LOCO) LMM (`scripts/run_loco_gwas.py`)**:
  - Computed 5 chromosome-specific kinship matrices excluding the target chromosome ($\mathbf{K}_{-1} \dots \mathbf{K}_{-5}$) using 397,850 LD-pruned variants.
  - Scanned 1,031,617 SNPs across all 5 chromosomes in GEMMA.
  - **Eliminated proximal contamination**: Genome-wide significant hits expanded from 8 to **34 SNPs**.
  - **Genomic Inflation**: $\lambda_{GC} = 1.1077$ (well calibrated).
  - Peak $P$-values sharpened by almost $100\times$:
    - ***VIN3*** (`5:23288908:T:A`): sharpened from $8.25 \times 10^{-10}$ to **$1.30 \times 10^{-11}$** ($\beta = +8.57$ days).
    - ***DOG1*** (`5:18590327:A:C`): sharpened from $2.46 \times 10^{-9}$ to **$8.70 \times 10^{-11}$** ($\beta = -3.45$ days).
    - ***FT*** (`1:24339560:G:T`): confirmed at $2.76 \times 10^{-9}$ ($\beta = -2.58$ days).
    - Discovered an additional locus on Chr 1 at 9.02 Mb (`1:9025070:A:T`, $P = 6.30 \times 10^{-9}$, $\beta = +19.4$ days).
- **Model 4: Stepwise Conditioning (MLMM)**:
  - Extracted sample dosages for top leads (*FT*, *DOG1*, *VIN3*) into `data/processed/ft10_covars_pcs_and_leads.txt`.
  - Conditioned GWAS in PLINK 2 on 10 PCs + the 3 lead QTLs.
  - Successfully unmasked ***FLC* (*FLOWERING LOCUS C*, Chr 5: 3,188,327)** as the top global lead peak with **$P = 1.12 \times 10^{-22}$**, and revealed a secondary peak on Chr 2 (`2:9583341:T:A`, $P = 6.49 \times 10^{-18}$).
- **4-Model Methodology Figure (`results/gwas/gwas_4model_methodology_comparison.png`)**:
  - Synthesized the complete progression (A: GLM, B: Whole-Genome LMM, C: LOCO LMM, D: Conditioned MLMM) into a unified publication graphic.
- **Publication-Grade Aligned LOCO Manhattan & QQ Figure**:
  - Created `scripts/plot_clear_loco_manhattan.py` producing a high-impact, side-by-side figure:
    - **Panel A (Manhattan, 75% width)**: Alternating Slate/Azure chromosomal tracks, bold vermilion significant points, and highlighted callout badges for ***FT*** ($-\log_{10} P = 8.56$), ***DOG1*** ($-\log_{10} P = 10.06$), and ***VIN3*** ($-\log_{10} P = 10.89$).
    - **Panel B (QQ Plot, 25% width)**: Aligned to the exact same vertical $-\log_{10} P$ axis range ($0 \dots 12.5$), with 95% concentration band and $\lambda_{GC} = 1.108$ badge.
  - Saved outputs: `results/gwas/ft10_loco_manhattan_clear.png` and interactive `results/gwas/ft10_loco_manhattan_clear.html`.

## 2026-10-31: Individual Polygenic Score (PGS / GBLUP) Assignment & Genomic Prediction

- **Predictive Architecture & Models (`scripts/compute_polygenic_scores.py`)**:
  - **4-QTL Additive Model (*FT*, *DOG1*, *VIN3*, *FLC*)**:
    - 5-fold cross-validated $R^2 = \mathbf{42.93\%}$ (Pearson $r = 0.655$, $\text{MAE} = 10.72$ days).
    - Demonstrates that four major-effect loci control over 40% of flowering variation alone.
  - **Whole-Genome GBLUP (Empirical Kinship)**:
    - 5-fold cross-validated $R^2 = \mathbf{66.64\%}$ (Pearson $r = \mathbf{0.8173}$, $\text{MAE} = 7.72$ days).
    - Captures the complete infinitesimal background, approaching the theoretical maximum bounded by narrow-sense heritability ($h^2 = 88.06\%$).
  - **Hybrid Model (4 QTLs + Polygenic Residual)**:
    - 5-fold cross-validated $R^2 = \mathbf{65.03\%}$ (Pearson $r = 0.811$, $\text{MAE} = 8.07$ days).
- **Individual Accession Table (`data/processed/athal_1003_polygenic_scores.csv`)**:
  - Assigned out-of-sample cross-validated predictions (`PGS_gblup_cv`), full genomic estimated breeding values (`PGS_gblup_full`), 4-QTL scores (`PGS_4qtl_cv`), and residual error terms for all 1,003 accessions.
- **Geographic Burden & Visuals (`results/prs/`)**:
  - `polygenic_score_evaluation.png`: Publication 2-panel figure:
    - Panel A: Observed vs Predicted FT10 scatter plot with linear fit and identity line.
    - Panel B: Distribution of genetic breeding values across the 9 geographic lineages, demonstrating the adaptive cline from early-flowering Relict/Iberian accessions to late-flowering Northern Sweden accessions.
  - `polygenic_score_evaluation.html`: Interactive Plotly explorer with hover metadata.
  - `pgs_model_metrics.csv`: 5-fold CV evaluation metrics table.






## 2026-10-31: Macro-Regional ANOVA and Corrected Multiple Comparisons for Polygenic Scores

- **Statistical Hypothesis & Pipeline (`scripts/anova_prs_populations.py`)**:
  - Evaluated whether whole-genome polygenic scores (5-fold cross-validated GBLUP days to flowering, $N = 1,003$ accessions) vary systematically across the 10 1001 Genomes macro-regions / genetic lineages: Western Europe ($n=104$), Central Europe ($n=177$), Germany ($n=90$), Italy/Balkans/Caucasus ($n=92$), Asia ($n=69$), Admixed ($n=136$), Spain ($n=110$), Relict ($n=25$), South Sweden ($n=139$), and North Sweden ($n=61$).
- **One-Way ANOVA Results (`results/prs/pgs_population_anova_stats.csv`)**:
  - **Fisher's One-Way ANOVA**: **$F(9, 993) = 136.11$**, **$P = 1.55 \times 10^{-166}$**.
  - **Variance Explained (Effect Size)**: **$\eta^2 = 55.23\%$** ($\omega^2 = 54.80\%$), confirming that more than half of the total variance in whole-genome genetic flowering potential is partitioned between geographic lineages.
  - **Welch's Heteroscedastic ANOVA (Alexander-Govern)**: **$A = 638.03$**, **$P = 1.43 \times 10^{-131}$** (accounting for variance heterogeneity identified by Levene's test, $W = 11.18, P = 9.33 \times 10^{-17}$).
- **Pairwise Multiple Comparisons (45 Contrasts, `results/prs/pgs_population_pairwise_comparisons.csv`)**:
  - **Tukey's Honestly Significant Difference (HSD)**: **34 / 45** pairwise contrasts are statistically significant at FWER $\alpha = 0.05$.
  - **Benjamini-Hochberg (BH) FDR**: **38 / 45** pairwise contrasts are statistically significant at $q < 0.05$.
  - **Bonferroni Adjustment**: **33 / 45** pairwise contrasts are statistically significant at $P_{\text{adj}} < 0.05$.
  - **Holm-Bonferroni Step-Down**: **35 / 45** pairwise contrasts are statistically significant at $P_{\text{adj}} < 0.05$.
- **Biological & Evolutionary Clines (`results/prs/pgs_population_group_summary.csv`)**:
  - **Latitudinal & Photoperiod Differentiation**:
    - **Earliest Flowering Lineages**: Western Europe (mean $= 72.21$ d, $95\%$ CI: $[71.02, 73.40]$) and Central Europe (mean $= 74.49$ d, $95\%$ CI: $[73.54, 75.44]$).
    - **Intermediate Lineages**: Germany ($77.93$ d), Italy/Balkans/Caucasus ($81.21$ d), Asia ($81.36$ d), Admixed ($82.55$ d), Spain ($83.48$ d), and Relicts ($89.29$ d).
    - **High-Latitude Vernalization Clines**: South Sweden ($97.91$ d, $95\%$ CI: $[96.54, 99.29]$) and North Sweden ($111.07$ d, $95\%$ CI: $[107.93, 114.22]$).
    - The most extreme divergence occurs between North Sweden and Western Europe (**$\Delta = +38.86$ days**, Welch $t(77.9) = 23.09$, $P_{\text{FDR}} = 1.58 \times 10^{-35}$, Cohen's $d = 4.37$).
  - **Non-Significant Contrasts (Stabilized Phenotypic Burden)**:
    - Intermediate populations form a statistically homogeneous cluster with no significant differences: Italy/Balkans vs Asia ($+0.15$ d, $q = 0.896$), Admixed vs Spain ($+0.93$ d, $q = 0.549$), Asia vs Admixed ($+1.20$ d, $q = 0.314$), Italy/Balkans vs Admixed ($+1.34$ d, $q = 0.302$), and Italy/Balkans vs Spain ($+2.27$ d, $q = 0.144$).
- **Visuals & Tables (`results/prs/`)**:
  - `pgs_population_contrasts.png`: Multi-panel publication figure (300 DPI) with population distributions/means, lower-triangular pairwise $-\log_{10}(P_{\text{FDR}})$ significance matrix, and contrast volcano plot.
  - `pgs_population_contrasts.html`: Interactive dual-panel boxplot and pairwise difference heatmap explorer.
