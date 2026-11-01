#!/usr/bin/env bash
# Download the 1001 Genomes v3.1 joint VCF and verify its checksum.
set -euo pipefail
cd "$(dirname "$0")/../data/vcf"
B=https://1001genomes.org/data/GMI-MPI/releases/v3.1
F=1001genomes_snp-short-indel_only_ACGTN.vcf.gz
curl -fLC - -O "$B/$F"
curl -fL -O "$B/$F.md5"
md5sum -c "$F.md5"
