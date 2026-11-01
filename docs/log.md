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





