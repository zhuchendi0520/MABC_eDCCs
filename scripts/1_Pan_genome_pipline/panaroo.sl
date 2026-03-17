#!/bin/bash

# ==============================================================================
# SLURM DIRECTIVES (Adjust these based on your cluster)
# ==============================================================================
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 24            # Request 24 total cores (e.g., for parallel jobs)
#SBATCH -t 1-00:00:00
#SBATCH --mem=50G
#SBATCH --job-name=Panaroo
#SBATCH --error=slurm_%j.err
#SBATCH --output=slurm_%j.out

# ==============================================================================
# SCRIPT SETUP
# ==============================================================================
set -euo pipefail # Exit immediately if a command exits with a non-zero status
# Uncomment the following lines if needed for R or grep preparation steps:
# grep -v NonDCC ../NonDCC_retype/Mab_random30.txt > Mab_random30_DCC.txt
# Rscript sample_flt.R Mab_random30_DCC.txt ../sample_info/bam_list_depth.txt

# Usage: ./panaroo_workflow.sh Mab_150.txt
export SAMPLE_LIST="$1"     # e.g., Mab_150.txt (The file containing sample names)
export OUTDIR="results_panaroo"
# THREADS is the total threads for Panaroo (Step 2)
export THREADS=24 
export CORE_T=0.95

# Prokka Parallel Settings
# Note: Ensure PROKKA_THREADS * PARALLEL_JOBS <= SBATCH -n (24 in this example)
export PROKKA_THREADS=4     # Threads for a single Prokka instance
export PARALLEL_JOBS=6     # Number of simultaneous Prokka jobs

# Derive a suffix based on the sample list filename
export LIST_SUFFIX="${SAMPLE_LIST%.txt}"

echo "[INFO] Input Sample List: $SAMPLE_LIST"
echo "[INFO] Output base directory: $OUTDIR"
echo "[INFO] Panaroo total threads: $THREADS"
echo "[INFO] Prokka parallel jobs: $PARALLEL_JOBS (using $PROKKA_THREADS cpus each)"

# Create necessary directories
# We use LIST_SUFFIX to group outputs from this specific run
export PROKKA_BASE_DIR="${OUTDIR}/prokka_out"
export PROKKA_SPECIFIC_DIR="${OUTDIR}/prokka_out_${LIST_SUFFIX}"
export PANAROO_OUT_DIR="${OUTDIR}/panaroo_out_${LIST_SUFFIX}"


mkdir -p "$PROKKA_BASE_DIR" "$PANAROO_OUT_DIR"

# ==============================================================================
# 1) Prokka Annotation (Parallel Execution)
# ==============================================================================
echo ""
echo "=================================================="
echo "[STEP 1] Prokka Annotation (Parallel)"
echo "=================================================="

# Get the list of sample IDs (assuming $2 contains the sample name)
SAMPLES=$(awk '{print $2}' "$SAMPLE_LIST")

# Define the function to be run in parallel
process_sample() {
    # Localize variables to the function's subshell
    local sample=$1
    local fa="../assemble/${sample}/scaffolds.fasta"
    local o="${PROKKA_BASE_DIR}/${sample}" # Use the consolidated directory
    local rename_fa="${sample}_rename.fa"
    local flt_fa="${sample}_flt.fa"
    
    echo "--- Starting job for $sample ---"

    # 1. Check if the assembly file exists
    if [[ ! -s "$fa" ]]; then
        echo "[WARN] Assembly not found: $fa"
        return 1
    fi
    
    # 2. Check if annotation is already complete (skip)
    if [[ -s "$o/${sample}.gff" ]]; then
        echo "[SKIP] $sample already annotated."
        cp -r "$o" "$PROKKA_SPECIFIC_DIR/"
        return 0
    fi

    # 3. Create sample-specific output directory
    mkdir -p "$o" 

    # 4. Filter FASTA by length (requires filter_fasta_by_length.py)
    python filter_fasta_by_length.py "$fa" "$flt_fa" 

    # 5. Sanitize sequence IDs (requires sanitize_fa_ids.py)
    python sanitize_fa_ids.py "$flt_fa" "$rename_fa" "${sample}_idmap.tsv"

    # 6. Clean up intermediate file
    rm "$flt_fa"

    # 7. Run Prokka
    prokka --outdir "$o" \
           --prefix "$sample" \
           --cpus "$PROKKA_THREADS" \
           --kingdom Bacteria \
           --force "$rename_fa" \
           --quiet
          
           
    # 8. Clean up temporary files generated in the working directory
    rm "$rename_fa" 
    cp -r "$o" "$PROKKA_SPECIFIC_DIR/"
    echo "--- Finished job for $sample ---"
}

# Export the function definition to make it available to parallel subshells
export -f process_sample

# Execute the function in parallel
echo "$SAMPLES" | parallel -j "$PARALLEL_JOBS" process_sample {} || true

echo "All Prokka jobs have been completed or skipped."

# ==============================================================================
# 2) Panaroo Pangenome Analysis
# ==============================================================================
echo ""
echo "=================================================="
echo "[STEP 2] Panaroo Pangenome Analysis"
echo "=================================================="

# Panaroo input GFF files are located in the consolidated directory
# The structure is $PROKKA_BASE_DIR/sample_id/sample_id.gff
shopt -s nullglob
GFF_INPUTS=( "$PROKKA_SPECIFIC_DIR"/*/*.gff )

if (( ${#GFF_INPUTS[@]} == 0 )); then
  echo "[ERROR] No GFF files found in $PROKKA_SPECIFIC_DIR"
  exit 1
fi

panaroo -i "${GFF_INPUTS[@]}" \
        -o "$PANAROO_OUT_DIR" \
        -t "$THREADS" \
        --clean-mode strict \
        --remove-invalid-genes \
        --merge_paralogs \
        -a core \
        --core_threshold "$CORE_T" \
        

CORE_ALN="${PANAROO_OUT_DIR}/core_gene_alignment.aln"
if [[ ! -s "$CORE_ALN" ]]; then
    echo "[ERROR] core_gene_alignment.aln not produced by Panaroo."
    exit 1
fi

echo ""
echo "=================================================="
echo "[DONE] Workflow Finished"
echo "=================================================="
echo " - Prokka results: $PROKKA_BASE_DIR"
echo " - Panaroo output: $PANAROO_OUT_DIR"
echo " - Final core alignment: $CORE_ALN"