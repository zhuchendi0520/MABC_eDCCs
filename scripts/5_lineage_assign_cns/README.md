# Evolutionary Path-Based Lineage Assignment from VarScan CNS Files

This directory contains scripts and clade-defining SNP marker sets for assigning *Mycobacterium abscessus* complex isolates to dominant circulating clones (DCCs) and emerging dominant circulating clones (eDCCs).

The workflow consists of four steps:

1. Map sequencing reads to the appropriate subspecies-specific reference genome.
2. Generate a consensus (`.cns`) file for each isolate with VarScan.
3. Convert each CNS file to a three-column typing SNP file.
4. Run the lineage-assignment script and marker set corresponding to the isolate's subspecies.

## Reference genomes

Each isolate must be mapped to the reference genome corresponding to its subspecies.

| Subspecies | Short name | Reference accession | Reference genome |
|---|---:|---|---|
| *M. abscessus* subsp. *abscessus* | ABS | NZ_CP034181.1 | *Mycobacteroides abscessus* strain GZ002 chromosome, complete genome |
| *M. abscessus* subsp. *massiliense* | MAS | NZ_AP014547.1 | *Mycobacteroides abscessus* subsp. *massiliense* CCUG 48898 = JCM 15300 chromosome, complete genome |
| *M. abscessus* subsp. *bolletii* | BOL | NZ_CP065265.1 | *Mycobacteroides abscessus* subsp. *bolletii* strain GD91 chromosome, complete genome |

The corresponding FASTA files are available in `../2_MAB_mapping/`:

```text
../2_MAB_mapping/M.abscessus.fna
../2_MAB_mapping/MAB.subsp.massiliense.fna
../2_MAB_mapping/MAB.subsp.bolletii.fna
```

Do not mix reference coordinate systems. The CNS file and typing marker file used for an isolate must be based on the same subspecies-specific reference genome.

## Requirements

- Python 3
- BWA
- SAMtools
- VarScan 2
- Paired-end Illumina FASTQ files

The commands below assume that the reads have already undergone adapter and quality trimming.

## Step 1: Map reads to the subspecies-specific reference

Set the sample name, FASTQ paths, reference genome, and VarScan JAR path:

```bash
SAMPLE="sample01"
R1="/path/to/sample01_R1.fastq.gz"
R2="/path/to/sample01_R2.fastq.gz"
VARSCAN_JAR="/path/to/VarScan.v2.3.9.jar"
THREADS=8
```

Select the correct reference:

```bash
# ABS
REF="../2_MAB_mapping/M.abscessus.fna"

# MAS
# REF="../2_MAB_mapping/MAB.subsp.massiliense.fna"

# BOL
# REF="../2_MAB_mapping/MAB.subsp.bolletii.fna"
```

Index each reference once before processing samples:

```bash
bwa index "${REF}"
samtools faidx "${REF}"
```

Map the paired-end reads, sort the alignments, and index the BAM file:

```bash
bwa mem \
  -t "${THREADS}" \
  -c 100 \
  -M \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "${REF}" "${R1}" "${R2}" \
  | samtools sort -@ "${THREADS}" -o "${SAMPLE}.sorted.bam" -

samtools index "${SAMPLE}.sorted.bam"
```

## Step 2: Generate a VarScan CNS file

Generate a pileup using minimum mapping quality 30 and minimum base quality 20:

```bash
samtools mpileup \
  -q 30 \
  -Q 20 \
  -B \
  -O \
  -f "${REF}" \
  "${SAMPLE}.sorted.bam" \
  > "${SAMPLE}.pileup"
```

Generate the consensus file with VarScan:

```bash
java -jar "${VARSCAN_JAR}" mpileup2cns "${SAMPLE}.pileup" \
  --min-coverage 3 \
  --min-avg-qual 20 \
  --min-var-freq 0.75 \
  --min-reads2 2 \
  --strand-filter 0 \
  > "${SAMPLE}.cns"
```

The resulting `${SAMPLE}.cns` file is the input for the next step.

## Step 3: Convert the CNS file to the typing format

Run:

```bash
python3 cns_to_typesnp.py "${SAMPLE}.cns"
```

The script writes:

```text
sample01.cns.typesnp
```

The output is a tab-delimited file with three columns:

```text
Position    Ref    Var
```

For positions at which the VarScan `Var` field is `.`, the script uses the reference allele as the sample allele.

Example:

```text
Position    Ref    Var
30579       A      G
32836       G      A
44839       C      T
```

## Step 4: Assign the DCC/eDCC lineage

Run the script and marker set corresponding to the isolate's subspecies.

### ABS isolates

```bash
python3 subtype_assign_abs_pathway.py \
  typing_SNP_ABS_pathway.txt \
  "${SAMPLE}.cns.typesnp" \
  > "${SAMPLE}.ABS.lineage.tsv"
```

The ABS script assigns the best-matching lineage when at least 95% of its defining SNPs are detected. Otherwise, the final category is `Non_DCC`.

### MAS isolates

```bash
python3 subtype_assign_mas_pathway.py \
  typing_SNP_MAS_pathway.txt \
  "${SAMPLE}.cns.typesnp" \
  > "${SAMPLE}.MAS.lineage.tsv"
```

The MAS script assigns the best-matching lineage when at least 95% of its defining SNPs are detected. Otherwise, the final category is `Non_DCC`.

### BOL isolates

```bash
python3 subtype_assign_bol_pathway.py \
  typing_SNP_BOL_pathway.txt \
  "${SAMPLE}.cns.typesnp" \
  > "${SAMPLE}.BOL.lineage.tsv"
```

The BOL script assigns the best-matching lineage when at least 95% of its defining SNPs are detected. Otherwise, the final category is `Non_DCC`.

## Assignment output

Each output line contains five tab-delimited fields:

```text
Sample    Best_match    Match_percentage    Matched/Total_markers    Assignment
```

Example:

```text
sample01.cns    DCC1    98.52%    (133/135)    DCC1
```

Field descriptions:

- `Sample`: input filename without its final extension.
- `Best_match`: lineage with the highest marker match percentage.
- `Match_percentage`: percentage of lineage-defining SNPs detected in the isolate.
- `Matched/Total_markers`: number of detected markers divided by the total number of markers for the best-matching lineage.
- `Assignment`: assigned DCC/eDCC name if the subspecies-specific threshold is met; otherwise `Non_DCC`.

If multiple lineages have exactly the same highest match percentage, the script reports one output line for each tied lineage.

## Batch processing

The following examples assume that all CNS files in a directory belong to the same subspecies.

### Convert all CNS files

```bash
for cns in /path/to/cns_files/*.cns; do
  python3 cns_to_typesnp.py "${cns}"
done
```

### Batch assignment for ABS

```bash
printf "Sample\tBest_match\tMatch_percentage\tMatched/Total_markers\tAssignment\n" \
  > ABS_lineage_assignments.tsv

for snp in /path/to/cns_files/*.cns.typesnp; do
  python3 subtype_assign_abs_pathway.py \
    typing_SNP_ABS_pathway.txt \
    "${snp}" \
    >> ABS_lineage_assignments.tsv
done
```

### Batch assignment for MAS

```bash
printf "Sample\tBest_match\tMatch_percentage\tMatched/Total_markers\tAssignment\n" \
  > MAS_lineage_assignments.tsv

for snp in /path/to/cns_files/*.cns.typesnp; do
  python3 subtype_assign_mas_pathway.py \
    typing_SNP_MAS_pathway.txt \
    "${snp}" \
    >> MAS_lineage_assignments.tsv
done
```

### Batch assignment for BOL

```bash
printf "Sample\tBest_match\tMatch_percentage\tMatched/Total_markers\tAssignment\n" \
  > BOL_lineage_assignments.tsv

for snp in /path/to/cns_files/*.cns.typesnp; do
  python3 subtype_assign_bol_pathway.py \
    typing_SNP_BOL_pathway.txt \
    "${snp}" \
    >> BOL_lineage_assignments.tsv
done
```

## Complete single-sample example

Example for an ABS isolate:

```bash
SAMPLE="sample01"
R1="/path/to/sample01_R1.fastq.gz"
R2="/path/to/sample01_R2.fastq.gz"
REF="../2_MAB_mapping/M.abscessus.fna"
VARSCAN_JAR="/path/to/VarScan.v2.3.9.jar"
THREADS=8

bwa mem \
  -t "${THREADS}" \
  -c 100 \
  -M \
  -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  "${REF}" "${R1}" "${R2}" \
  | samtools sort -@ "${THREADS}" -o "${SAMPLE}.sorted.bam" -

samtools index "${SAMPLE}.sorted.bam"

samtools mpileup \
  -q 30 -Q 20 -B -O \
  -f "${REF}" \
  "${SAMPLE}.sorted.bam" \
  > "${SAMPLE}.pileup"

java -jar "${VARSCAN_JAR}" mpileup2cns "${SAMPLE}.pileup" \
  --min-coverage 3 \
  --min-avg-qual 20 \
  --min-var-freq 0.75 \
  --min-reads2 2 \
  --strand-filter 0 \
  > "${SAMPLE}.cns"

python3 cns_to_typesnp.py "${SAMPLE}.cns"

python3 subtype_assign_abs_pathway.py \
  typing_SNP_ABS_pathway.txt \
  "${SAMPLE}.cns.typesnp" \
  > "${SAMPLE}.ABS.lineage.tsv"
```

## Important notes

- Determine the subspecies before lineage assignment.
- Use the ABS reference, marker set, and assignment script only for ABS isolates; use the corresponding MAS or BOL files for the other subspecies.
- Marker positions are reference-specific and cannot be transferred between reference genomes.
- Low sequencing depth, contamination, mixed infection, or poor mapping can reduce marker recovery and lead to a `Non_DCC` result.
- The assignment threshold is 95% for ABS, MAS, and BOL, as implemented in the scripts.
- The assignment scripts report results to standard output. Use `>` or `>>` to save them to a file.
