# Log

## Date: 2026-02-28
- Downloaded the list of studies on the species
- Generated a list for each accession that is genotyped
- Retrieved the phenotype database for each accession
- Wrote the "goals" document
- Created the mamba environment
- main goal for tomorrow will be to get the (test subset of) read files, reference file, and phenotype files all formatted well'

## Date: 2026-03-01
- Downloaded reference genome from [ncbi](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/) using select file source as "All" and select file types as "FASTA" and "JSONL".
- moved to ref/ folder
- unzipped reference genome and deleted .zip
- extracted gcf .fna file and deleted everything else
- in the script i will ultimately create that auto performs this analysis, it will automatically rename the reference genome to ref.fna. It is the future user's responsibility to make sure to keep track of the information (ID and source) of the ref genome. But for now I renamed it to ref.fna
