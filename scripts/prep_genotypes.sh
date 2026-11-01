#!/usr/bin/env bash
# Subset 1001 Genomes VCF to 1,003 accessions with FT10 phenotypes,
# filter to biallelic SNPs (MAF >= 0.01, missingness <= 0.1),
# and produce PLINK2 (PGEN) and PLINK1 (BED) datasets.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${ROOT_DIR}"

VCF="data/vcf/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
SAMPLES="data/processed/common_samples.txt"
OUT_PREFIX="data/processed/athal_1003_qc"

echo "=== Stage 1: Converting VCF and applying QC filters with PLINK2 ==="
echo "Input VCF: ${VCF}"
echo "Samples: ${SAMPLES} ($(wc -l < "${SAMPLES}") accessions)"
echo "Target output prefix: ${OUT_PREFIX}"

plink2 \
  --vcf "${VCF}" \
  --keep "${SAMPLES}" \
  --snps-only just-acgt \
  --min-alleles 2 \
  --max-alleles 2 \
  --geno 0.10 \
  --maf 0.01 \
  --set-all-var-ids '@:#:$r:$a' \
  --allow-extra-chr \
  --threads 8 \
  --memory 12000 \
  --make-pgen \
  --out "${OUT_PREFIX}"

echo "=== Stage 2: Exporting PLINK 1.9 binary fileset (BED/BIM/FAM) ==="
plink2 \
  --pfile "${OUT_PREFIX}" \
  --allow-extra-chr \
  --threads 8 \
  --memory 12000 \
  --make-bed \
  --out "${OUT_PREFIX}"

echo "=== QC and Conversion Complete ==="
ls -lh "${OUT_PREFIX}".*
