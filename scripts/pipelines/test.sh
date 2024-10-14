#!/bin/bash

# ---------------------------------------------------------
# K-CHOPORE (Keen Comprehensive High-throughput Omics Pipeline Organizer)
# An automated tool for analyzing Nanopore sequencing data.
#
# Integrates various essential steps from basecalling to epitranscriptomic modification detection.
#
# Developed by:
# Pelayo González de Lena Rodríguez, MSc
# FPI Severo Ochoa Fellow
# Cancer Epigenetics and Nanomedicine Lab | FINBA
# Systems Biology Lab | University of Oviedo
#
# Contact:
# LinkedIn: https://www.linkedin.com/in/biopelayo/
# GitLab: https://gitlab.com/bio.pelayo/
# Phone: +34 660 74 11 39
#
# License:
# This tool is distributed under the MIT License. See LICENSE file for details.
#
# ---------------------------------------------------------

# Config file path (adjust if located elsewhere)
CONFIG_FILE="/workspace/config/config.yml"

# Function to parse the YAML configuration file
parse_yaml() {
    local prefix=$2
    local s='[[:space:]]*'
    local w='[a-zA-Z0-9_]*'
    local fs=$(echo @|tr @ '\034')
    sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|$fs\2$fs\3|p"  $1 |
    awk -F$fs '{
        indent = length($1)/2;
        vname[indent] = $2;
        for (i in vname) {if (i > indent) {delete vname[i]}}
        if (length($3) > 0) {
            vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
            printf("%s%s%s=\"%s\"\n", "'$prefix'", vn, $2, $3);
        }
    }'
}

# Pretty print function for headers
pretty_print() {
    local message=$1
    local border_length=${#message}
    local border=$(printf '%*s' "$border_length" '' | tr ' ' '-')
    echo
    echo " $border"
    echo " $message"
    echo " $border"
    echo
}

# Parse the configuration file
if [ -f "$CONFIG_FILE" ]; then
    eval $(parse_yaml $CONFIG_FILE "config_")
else
    echo "Error: Configuration file not found at $CONFIG_FILE"
    exit 1
fi


# ------------------------------
# Print all parsed configuration variables
# ------------------------------
echo "DEBUG: Configuration Variables Parsed from config.yml:"
env | grep "^config_"  # Prints all variables starting with "config_"
echo "---------------------------------------------------------"


# Depuración: Verificar que la ruta esté siendo leída correctamente
echo "DEBUG: Path to sequencing summary file: $config_input_files_sequencing_summary"

# Debug Print: Configuration Variables
echo "DEBUG: Configuration Variables Loaded"
echo "DEBUG: config_input_files_reference_genome = $config_input_files_reference_genome"
echo "DEBUG: config_options_run_flair = $config_options_run_flair"
echo "---------------------------------------------------------"

# ------------------------------
# Welcome Message
# ------------------------------
pretty_print "K-CHOPORE Pipeline - Build 2024"
echo "Keen Comprehensive High-throughput Omics Pipeline Organizer"
echo "An automated tool for analyzing Nanopore sequencing data."
echo
echo "Developed by: Pelayo González de Lena Rodríguez, MSc"
echo "FPI Severo Ochoa Fellow | Cancer Epigenetics and Nanomedicine Lab | FINBA"
echo "Systems Biology Lab | University of Oviedo"
echo "Contact: LinkedIn - https://www.linkedin.com/in/biopelayo/"
echo "         GitLab - https://gitlab.com/bio.pelayo/"
echo "---------------------------------------------------------"

# ------------------------------
# Check R Version
# ------------------------------
pretty_print "Checking R Version"
R_VERSION=$(R --version | head -n 1 | awk '{print $3}')
echo "DEBUG: Installed R version: $R_VERSION"

REQUIRED_R_VERSION=$config_options_r_version_required

if [[ $(echo -e "$R_VERSION\n$REQUIRED_R_VERSION" | sort -V | head -n 1) == "$REQUIRED_R_VERSION" ]]; then
    echo "R version is compatible. Proceeding with package installation..."
    if ! Rscript -e 'install.packages("samplesizeCMH", repos="https://cloud.r-project.org")'; then
        echo "Error installing R package 'samplesizeCMH'. Please check your R installation and internet connection."
        exit 1
    fi
else
    echo "Incompatible R version ($R_VERSION). R version $REQUIRED_R_VERSION or higher is required."
    exit 1
fi

# ------------------------------
# Install and Configure Necessary Tools
# ------------------------------
pretty_print "Installing and Configuring Tools"

echo "Configuring PATH..."
export PATH=$PWD:$PWD/Scripts:$PATH || { echo "Error configuring PATH. Please check the directory path."; exit 1; }
echo "DEBUG: PATH configured as: $PATH"

# ------------------------------
# Normalize Reference Genome
# ------------------------------
pretty_print "Normalizing Reference Genome"

echo "Normalizing the reference genome..."
echo "DEBUG: Reference genome input: $config_input_files_reference_genome"
echo "DEBUG: Normalized genome output: $config_input_files_reference_genome_norm"
java -jar /usr/local/bin/picard.jar NormalizeFasta I=$config_input_files_reference_genome O=$config_input_files_reference_genome_norm

# ------------------------------
# Mapping and Modifications using minimap2 and ELIGOS
# ------------------------------
pretty_print "Mapping and Epitranscriptomic Modifications"

echo "Performing mapping with minimap2..."
echo "DEBUG: Mapping with reference genome: $config_input_files_reference_genome_norm and FASTQ: $config_input_files_fastq_file"
/usr/local/bin/minimap2 -ax map-ont $config_input_files_reference_genome_norm $config_input_files_fastq_file > $config_input_files_map_file
echo "DEBUG: Mapping output file: $config_input_files_map_file"

echo "Sorting and indexing BAM file..."
samtools sort $config_input_files_map_file > ${config_input_files_map_file%.sam}_sorted.bam
samtools index ${config_input_files_map_file%.sam}_sorted.bam
echo "DEBUG: Sorted BAM file: ${config_input_files_map_file%.sam}_sorted.bam"

echo "Executing epitranscriptomic modifications with ELIGOS..."
./eligos2 rna_mod -i ${config_input_files_map_file%.sam}_sorted.bam -reg $config_input_files_bed_file -ref $config_input_files_reference_genome -o results --pval 0.05 --oddR 5 --esb 0.2
echo "DEBUG: ELIGOS2 results directory: $(pwd)/results"

echo "Performing multi-sample analysis with ELIGOS..."
./eligos2 multi_samples_test --test_mods ./results/*_baseExt0.txt --ctrl_mods *_baseExt0.txt --prefix output_prefix_
echo "DEBUG: Multi-sample test output prefix: output_prefix_"

# ------------------------------
# Check if FLAIR Should Be Executed
# ------------------------------
if [ "${config_options_run_flair}" == "true" ]; then
    pretty_print "Running FLAIR for Transcriptome Construction"
    
    cd $config_input_files_flair_repo
    echo "DEBUG: Changed to FLAIR repository directory: $(pwd)"
    
    echo "Activating FLAIR virtual environment..."
    source k-chopo/bin/activate

    echo "Aligning samples with FLAIR..."
    python flair.py align -g $config_input_files_reference_genome_norm -r $config_input_files_fastq_file
    echo "DEBUG: FLAIR aligned samples"

    echo "Generating BED file..."
    samtools sort flair.aligned.unsorted.bam -o flair.aligned.sorted.bam
    samtools index flair.aligned.sorted.bam
    python ./bin/bam2Bed12.py -i flair.aligned.sorted.bam > flair.aligned.sorted.bed
    echo "DEBUG: Generated BED file: flair.aligned.sorted.bed"
    
    echo "Collapsing alignments with FLAIR..."
    python flair.py collapse -g $config_input_files_reference_genome_norm -r $config_input_files_fastq_file -q flair.aligned.sorted.bed -f $config_input_files_gtf_file
    echo "DEBUG: Collapsed alignments"

    echo "Quantifying isoforms..."
    python flair.py quantify -g $config_input_files_reference_genome_norm -r $config_input_files_fastq_file -q flair.aligned.sorted.bed -f $config_input_files_gtf_file --output $config_input_files_flair_output_dir
    echo "DEBUG: Isoform quantification output directory: $config_input_files_flair_output_dir"
    
    echo "Summarizing counts and performing Fisher test..."
    python flair.py merge --inputs $config_input_files_counts_matrix --output $config_input_files_summarized_counts_matrix
    python flair.py fisher --inputs $config_input_files_summarized_counts_matrix --output $config_input_files_diff_isos_file
    echo "DEBUG: Fisher test output file: $config_input_files_diff_isos_file"
fi

# ------------------------------
# Run pycoQC
# ------------------------------
pretty_print "Running pycoQC for Quality Analysis"

if [ -f "$config_input_files_sequencing_summary" ]; then
    pycoQC -i $config_input_files_sequencing_summary -o /path/to/pycoQC_output
    echo "DEBUG: pycoQC output directory: /path/to/pycoQC_output"
else
    echo "Sequencing summary file not found. Check the file and its path."
    exit 1
fi

# ------------------------------
# Additional Isoform and Multivariable Analysis
# ------------------------------
pretty_print "Additional Isoform and Multivariable Analysis"

echo "Executing additional isoform analysis..."
python /path/to/isoformswitch_analyzer/isoformswitch_analyzer.py -i $config_input_files_iso_usage_file -o $config_input_files_isoformswitch_analyzer_output
echo "DEBUG: Isoform switch analysis output directory: $config_input_files_isoformswitch_analyzer_output"

echo "Performing PCA..."
python /path/to/pca_analysis.py -i $config_input_files_counts_matrix -o $config_input_files_pca_output
echo "DEBUG: PCA analysis output directory: $config_input_files_pca_output"

echo "Performing SPLS analysis..."
python /path/to/spls_analysis.py -i $config_input_files_counts_matrix -o $config_input_files_spls_output
echo "DEBUG: SPLS analysis output directory: $config_input_files_spls_output"

# ------------------------------
# Important Note on FASTA and GTF Name Matching
# ------------------------------
pretty_print "Important Note"

echo "Ensure that the names in the reference genome FASTA file match those in the GTF file to avoid issues during analysis."

# ------------------------------
# Completion Message
# ------------------------------
pretty_print "K-CHOPORE Pipeline - Completed"
echo "All steps have been successfully executed."

