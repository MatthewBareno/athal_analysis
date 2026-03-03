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