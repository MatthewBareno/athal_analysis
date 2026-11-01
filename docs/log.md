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





