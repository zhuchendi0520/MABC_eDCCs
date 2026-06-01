# MABC_eDCCs

This repository contains scripts and selected analysis files for the manuscript:

**Extensive Parallel Global Expansion and Ongoing Adaptive Evolution of *Mycobacterium abscessus***

The study analyzes global whole-genome sequencing data from the *Mycobacterium abscessus* complex (MAB) to identify dominant circulating clones (DCCs), emerging DCCs (eDCCs), their geographic spread, temporal expansion, clade assignment markers, and genes under recent adaptive selection.

## Repository overview

```text
MABC_eDCCs/
|-- data/
|   |-- BEAST/
|   |-- Gubbins FASTA/
|   `-- Phylogenetic trees/
`-- scripts/
    |-- 1_Pan_genome_pipline/
    |-- 2_MAB_mapping/
    |-- 3_pNS/
    |-- 4_MAB_gene_translate/
    |-- 5_lineage_assign_cns/
    `-- 6_R_script/
```

## Data directories

- `data/BEAST/`: BEAST XML files, posterior tree files, and summarized tree files used for temporal reconstruction and demographic inference of DCC and eDCC lineages.
- `data/Gubbins FASTA/`: recombination-filtered polymorphic-site FASTA files generated after core-genome alignment and Gubbins filtering. These files are used for clade-level phylogenetic reconstruction.
- `data/Phylogenetic trees/`: phylogenetic trees for the three MAB subspecies and PastML-related phylogeographic reconstruction results.

## Script directories

### `scripts/1_Pan_genome_pipline/`

Scripts in this directory are used for clade-level core-genome construction and recombination filtering, corresponding to the manuscript sections on DCC/eDCC core-genome phylogenies.

- `panaroo.sl`: SLURM workflow for annotating assemblies with Prokka and building core-genome alignments with Panaroo. The script takes a sample list, runs Prokka in parallel, and then runs Panaroo in strict mode to generate a core-gene alignment.
- `gubbins.sl`: SLURM script for running Gubbins on Panaroo core-genome alignments. This step detects and removes recombinant regions before downstream phylogenetic reconstruction.

### `scripts/2_MAB_mapping/`

Scripts and reference files in this directory are used for read mapping and SNP calling against subspecies-specific MAB reference genomes, corresponding to the manuscript Methods sections on SNP calling and repetitive-region filtering.

- `mapping_bwamem.M.abscessus.sh`: generates commands for trimming reads, mapping to the *M. abscessus* reference genome with BWA-MEM, sorting/indexing BAM files, calling SNPs with VarScan, filtering repetitive regions, and extracting fixed SNPs.
- `mapping_bwamem.MAB.subsp.massiliense.sh`: same mapping and SNP-calling workflow for *M. abscessus* subsp. *massiliense*.
- `mapping_bwamem.MAB.subsp.bolletii.sh`: same mapping and SNP-calling workflow for *M. abscessus* subsp. *bolletii*.
- `mapping_bwamem.ATCC19977.sh`: mapping workflow using the ATCC19977 reference.
- `M.abscessus.fna`, `MAB.subsp.massiliense.fna`, `MAB.subsp.bolletii.fna`: reference genome FASTA files used by the mapping scripts.
- `*_repeat_PEPPE_phage_transposase.loci*.txt`: repeat-region masks used to exclude SNPs in repetitive loci, phage regions, PE/PPE-like regions, and transposase-associated regions.

### `scripts/3_pNS/`

Scripts in this directory are used for gene-level pN/pS or pNS analysis, corresponding to the manuscript section on adaptive selection during recent DCC/eDCC expansion.

- `format-snppar.py`: reformats SNPPar annotation output into the input table required by the pNS scripts. It classifies SNPs as synonymous, nonsynonymous, or intergenic and records the wild-type codon.
- `pNS.py`: original pNS calculation script. It estimates expected synonymous and nonsynonymous mutation opportunities under a codon-level mutation model, then calculates gene-level pNS values.
- `pNS_setsynto1.py`: modified pNS script used when a gene has no observed synonymous mutation. In that case, the synonymous count is set to 1 to avoid division-by-zero problems while retaining genes with nonsynonymous-only mutations.
- `pNS_readme.txt`: brief usage note for the pNS scripts.

### `scripts/4_MAB_gene_translate/`

Scripts and annotation tables in this directory are used to annotate SNPs with gene location and coding consequence, corresponding to the manuscript analyses of synonymous/nonsynonymous mutations and candidate adaptive genes.

- `1_M.abscessus_Annotation.py`: Python SNP annotation script for the *M. abscessus* reference. It maps mutations to genes or intergenic regions, determines codon position, and classifies coding mutations as synonymous or nonsynonymous.
- `1_M.abscessus_Annotation.pl`: Perl version of the *M. abscessus* SNP annotation script.
- `1_massiliense_Annotation.pl`: SNP annotation script for subsp. *massiliense*.
- `1_bolletii_Annotation.pl`: SNP annotation script for subsp. *bolletii*.
- `2_M.abscessus_20220725`, `2_MAB.subsp.massiliense_20221024`, `2_MAB.subsp.bolletii_20220725`: gene annotation tables used by the annotation scripts.
- `3_genetic_codes`: codon translation table used to classify mutations as synonymous or nonsynonymous.

### `scripts/5_lineage_assign_cns/`

Scripts and marker files in this directory implement the evolutionary path-based typing framework described in the manuscript. They assign isolates to DCC/eDCC lineages using clade-defining SNP markers.

- `subtype_assign_abs_pathway.py`: assigns subsp. *abscessus* isolates to DCC/eDCC lineages by comparing an isolate SNP file against clade-defining SNP markers.
- `subtype_assign_mas_pathway.py`: same lineage assignment workflow for subsp. *massiliense*.
- `subtype_assign_bol_pathway.py`: same lineage assignment workflow for subsp. *bolletii*.
- `typing_SNP_ABS_pathway.txt`: clade-defining SNP marker set for subsp. *abscessus*.
- `typing_SNP_MAS_pathway.txt`: clade-defining SNP marker set for subsp. *massiliense*.
- `typing_SNP_BOL_pathway.txt`: clade-defining SNP marker set for subsp. *bolletii*.
- `trans.py`: helper script for cleaning marker tables by replacing missing reference-base fields with the available base information and writing a tab-delimited output file.

### `scripts/6_R_script/`

This directory contains the R Markdown workflow used to generate summary plots and statistical analyses for the manuscript.

- `abscessus_analyze.Rmd`: R analysis and plotting notebook for manuscript figures and supplementary figures, including global DCC/eDCC distribution heatmaps, expansion-time plots, BEAST skyline summaries, population-size regressions, lineage marker heatmaps, positive-selection plots, sample metadata summaries, rarefaction analyses, recombination pattern plots, unfixed mutation summaries, and nonsynonymous mutation burden comparisons.

## Relationship to the manuscript

The scripts correspond to the main analytical components of the study:

- Genome processing and SNP calling: `scripts/2_MAB_mapping/`
- Core-genome alignment and recombination filtering: `scripts/1_Pan_genome_pipline/`
- Subspecies and clade-level phylogenetic analyses: `data/Phylogenetic trees/`, `data/Gubbins FASTA/`
- Temporal and demographic reconstruction: `data/BEAST/`
- Evolutionary path-based clade assignment: `scripts/5_lineage_assign_cns/`
- Mutation annotation and selection analysis: `scripts/3_pNS/`, `scripts/4_MAB_gene_translate/`
- Manuscript statistical analysis and plotting: `scripts/6_R_script/`

## Notes

- Some scripts contain absolute paths from the original analysis environment and should be edited before running on a new system.
- Several workflows were designed for an HPC/SLURM environment.
- Raw sequencing reads and large intermediate files are not included in this repository.
- The pNS implementation was modified from the script in `https://github.com/swisstph/TBRU_serialTB`.
