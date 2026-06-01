# MABC_eDCCs

This repository contains scripts and selected input files used to reproduce analyses of the global population structure, evolution, and transmission dynamics of the *Mycobacterium abscessus* complex.

The workflow is organized around major analytical steps from genome processing to phylogenetic and evolutionary analysis, including pan-genome construction, read mapping, recombination-aware SNP processing, clade assignment, pNS estimation, and downstream visualization.

## Repository structure

```text
MABC_eDCCs/
├── data/
│   ├── BEAST/                         BEAST XML, posterior tree, and MCC tree files
│   ├── Gubbins FASTA/                 Recombination-filtered polymorphic-site FASTA files
│   └── Phylogenetic trees/            Subspecies trees and PastML ancestral-state results
└── scripts/
    ├── 1_Pan_genome_pipline/          Panaroo and Gubbins batch scripts
    ├── 2_MAB_mapping/                 Reference genomes, repeat masks, and BWA-MEM mapping scripts
    ├── 3_pNS/                         SNP annotation formatting and gene-level pNS scripts
    ├── 4_MAB_gene_translate/          Gene annotation and translation helper scripts
    ├── 5_lineage_assign_cns/          Subspecies-specific lineage assignment scripts
    └── 6_R_script/                    R Markdown analysis and plotting workflow
```

## Main workflow

1. Build and filter the pan-genome using the scripts in `scripts/1_Pan_genome_pipline/`.
2. Map sequencing reads to the relevant *M. abscessus* complex reference genome with the scripts in `scripts/2_MAB_mapping/`.
3. Detect and filter recombination with Gubbins. Filtered polymorphic-site FASTA files are stored in `data/Gubbins FASTA/`.
4. Assign lineages using the subspecies-specific SNP marker scripts in `scripts/5_lineage_assign_cns/`.
5. Estimate gene-level pNS values from annotated SNP tables using `scripts/3_pNS/`.
6. Run downstream statistical analysis and visualization with `scripts/6_R_script/abscessus_analyze.Rmd`.

## pNS analysis

The pNS workflow estimates the ratio of observed nonsynonymous to synonymous mutations after correcting for the expected nonsynonymous and synonymous mutation opportunities under a mutation spectrum model.

Relevant scripts:

- `scripts/3_pNS/format-snppar.py`: converts SNP annotation output into the CSV format required by the pNS script.
- `scripts/3_pNS/pNS.py`: original gene-level pNS calculation script.
- `scripts/3_pNS/pNS_setsynto1.py`: modified pNS script that sets the synonymous count to 1 when a gene has no observed synonymous mutation, preventing division by zero for genes with nonsynonymous-only mutations.

Example usage:

```bash
python scripts/3_pNS/format-snppar.py snp_annotation.txt > sample_ann.csv
python scripts/3_pNS/pNS_setsynto1.py sample_ann.csv
```

The pNS script writes a corresponding `*_gene.csv` output file.

## Example pNS output: L3i1.csv

`L3i1.csv` is a gene-level pNS result table generated from annotated SNP data. It contains 3,357 rows, including 3,356 gene entries after excluding the header-like `Gene` row. Across those gene entries, the table summarizes 14,144 coding mutations, including 5,516 synonymous and 8,628 nonsynonymous mutations. A total of 2,913 genes have a finite pNS value.

Columns:

| Column | Description |
| --- | --- |
| `GENE` | Gene identifier, usually an `Rv` locus tag. |
| `pNS` | Observed pNS value, calculated as `(OBSERVED_NSY / EXPECTED_NSY) / (OBSERVED_SYN / EXPECTED_SYN)`. |
| `OBSERVED_SYN` | Number of observed synonymous SNPs in the gene. |
| `OBSERVED_NSY` | Number of observed nonsynonymous SNPs in the gene. |
| `EXPECTED_SYN` | Expected synonymous mutation count under the codon-level mutation model. |
| `EXPECTED_NSY` | Expected nonsynonymous mutation count under the codon-level mutation model. |
| `NEUTRAL_SYN` | Synonymous count sampled from the neutral binomial model. |
| `NEUTRAL_NSY` | Nonsynonymous count sampled from the neutral binomial model. |
| `pNS_NEUTRAL` | Neutral pNS value from the simulated synonymous/nonsynonymous counts. |
| `TOTAL` | Total coding SNP count summarized for the gene. |

Interpretation:

- `pNS > 1` indicates an excess of nonsynonymous changes relative to synonymous changes after accounting for expectation.
- `pNS < 1` indicates fewer nonsynonymous changes than expected relative to synonymous changes.
- Empty or `NaN` pNS values usually arise when the ratio cannot be computed, for example because one side of the observed or expected ratio is zero.

## Dependencies

The scripts use a mixture of Python, Perl, R, and external bioinformatics tools. The pNS scripts require:

- Python 2 syntax for `format-snppar.py`
- Python with `pandas`, `numpy`, and `scipy` for `pNS.py` and `pNS_setsynto1.py`

Other workflow steps may require tools such as BWA-MEM, Panaroo, Gubbins, BEAST, and R packages used by `scripts/6_R_script/abscessus_analyze.Rmd`.

## Notes

- Paths in batch scripts may need to be adjusted before running on a different cluster or workstation.
- Large intermediate files, raw sequencing reads, and some derived outputs may not be included in this repository.
- The pNS implementation was modified from `https://github.com/swisstph/TBRU_serialTB`.
