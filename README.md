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
│   └── sample_map.html            # Interactive geographic distribution of accessions
└── scripts/
    ├── download_data.sh           # Downloads 1001G v3.1 VCF and verifies MD5
    ├── prepare_phenotypes.py      # Deduplicates & formats phenotypes for the 1,003 accessions
    ├── prep_genotypes.sh          # Subsets VCF, applies QC filters, and creates PLINK files
    └── map_samples.py             # Generates Plotly map of accession locations
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

## Analysis Roadmap

1. **[Completed] Data Prep:** VCF subsetting, biallelic SNP extraction, missingness & MAF filtering, phenotype deduplication.
2. **[Next] Population Structure:**
   - LD pruning (`plink2 --indep-pairwise 50 10 0.2`)
   - Principal Component Analysis (`plink2 --pca`)
   - Population ancestry modeling (ADMIXTURE / k-means)
   - Pairwise $F_{ST}$ and heterozygosity across geographic origins
3. **GWAS:**
   - Linear model with PC covariates (PLINK 2 / statsmodels)
   - Linear Mixed Model / kinship matrix (GEMMA / GCTA)
   - Manhattan and QQ plots, significant peak identification
4. **Quantitative & Evolutionary Genetics:**
   - Polygenic scores (PRS) and additive genetic variance ($V_A$)
   - Narrow-sense heritability ($h^2$)
   - Divergent selection tests ($Q_{ST} - F_{ST}$ / $Q_X$)

---

## References & Data Sources
- **1001 Genomes Consortium** (2016). *1,135 Genomes Reveal the Global Pattern of Polymorphism in Arabidopsis thaliana*. Cell, 166(2), 481-491.
- **AraPheno:** FT10 flowering time study 12 ([AraPheno citation](https://arapheno.1001genomes.org/faq/cite/)).
