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
- Phenotype downloaded and moved to data/pheno/ directory
- confirmed that the first column in this phenotype data refers to the accession ID on the website
- Worked out the logic necessary to perform this project. Will resume later.