# athal_analysis
Population and quantitative genomic analyses of _Arabidopsis thaliana_ (1001 Genomes).

## Plan
Use the pre-called 1001 Genomes v3.1 joint VCF instead of aligning reads and calling variants.

- **Phenotype:** FT10 (days to flowering at 10 °C), AraPheno study 12. 1,162 accessions.
- **Genotypes:** `1001genomes_snp-short-indel_only_ACGTN.vcf.gz` (19.2 GB, 1,135 accessions, TAIR10).
- **Analysis set:** the 1,003 accessions with both FT10 and genotypes (159 FT10 accessions are in no 1001G release; ignored for now).

## Analyses
1. **Prep:** subset VCF to the 1,003 accessions, biallelic SNPs, MAF/missingness filter, convert to PLINK.
2. **Population structure:** PCA, k-means / ADMIXTURE, phylogenetic tree, Fst matrix, heterozygosity by population (and latitude/longitude).
3. **GWAS:** compare linear model + PCs vs. mixed model; Manhattan and QQ plots; table of SNPs above an effect-size threshold.
4. **Quantitative genetics:** polygenic scores per individual and population, heritability, test for divergent selection (Qst-Fst / Q_X).

Methods within each step are meant to be compared; nothing is final.

## Layout
| Path | Contents |
|---|---|
| `data/accession_list.csv` | 1001G master list (id, name, country, lat/lon, ...) |
| `data/pheno/phenotype.csv` | FT10 / FT16 per accession |
| `data/vcf/` | downloaded VCF (gitignored) |
| `data/ref/ref.fna` | TAIR10 reference (gitignored) |
| `scripts/` | analysis code; `download_data.sh` fetches the VCF |
| `envs/environment.yml` | software environment |
| `results/` | outputs and figures |
| `docs/log.md` | dated lab notes |

## Getting the data
`bash scripts/download_data.sh` (about 19 GB; run on a fast connection).

## Notes
- Phenotype `accession_id` equals the `id` column in `accession_list.csv` and the VCF sample names.
- `phenotype.csv` has one duplicated accession; dedupe before use.
- Cite AraPheno: https://arapheno.1001genomes.org/faq/cite/
