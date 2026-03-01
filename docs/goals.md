# Project Goals - Athal Analysis

## Overview
Perform population genetic analysis and a genome-wide association study (GWAS) on *Arabidopsis thaliana* samples.

## Inputs
- Read files (paired-end: R1 and R2)
- Reference genome
- Annotation index file (optional, as config values can be indirectly discovered anyways)
- Phenotype data

## Analysis Goals

### 1. Alignment / Mapping
- I want to take each set of reads and align them to the reference genome to generate a full construction of each accesion genome

### 2. Variant Calling
- With each genome constructed, we will now generate VCF files for each individual, and then combine them. All analyses will run on this combined, population level vcf or package-concordant conversions of it.

### 3. Population Structure & Genetic Analysis
- **PCA**: Principal component analysis for population stratification
- **K-means Clustering**: Identify genetic subpopulations
- **Phylogenetic Tree**: Construct evolutionary relationships between samples
- **Other Population genetics figures**:Fst matrix among populations, heterozygosity by population (and perhaps by longitude/lattitude?, other potential vairables can also be considered)

### 4. GWAS (Genome-Wide Association Study)
- Perform association analysis between SNPs and phenotypes
- **Manhattan Plot**: Visualize significant SNP associations across the genome
- **Summary Table**: Generate results table of SNPs with effect sizes ≥ threshold value
- For each phenotype, generate a polygenic score, and give each individual and population a score
- perform a test for divergent selection among populations in the phenotye
