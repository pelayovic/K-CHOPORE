# K-CHOPORE 🌟
**Keen Comprehensive High-throughput Omics Pipeline Organizer**

Just like a classic Asturian cachopo, K-CHOPORE is a hearty and satisfying tool for Nanopore sequencing data analysis! 🧬🍴 Dive in and enjoy the data feast with K-CHOPORE! 📊💻

---

## 📜 Overview

**K-CHOPORE** is an automated pipeline for comprehensive analysis of Nanopore sequencing data, designed to handle the entire process from basecalling to epitranscriptomic modification detection. By integrating various specialized tools, K-CHOPORE allows you to perform in-depth RNA-Seq data analysis, covering all critical steps and ensuring robust, reproducible results.

---

## 🔧 Pipeline Components

### 1. Basecalling with Guppy
- **Description**: Converts raw Nanopore sequencing signals into readable DNA sequences.
- **Command**:
  ```bash
  guppy_basecaller -i input_directory -s output_directory -c dna_r9.4.1_450bps_hac.cfg
### 2. Sequence Quality Check with pycoQC
- **Description**: Assesses the quality of generated sequences.
- **Installation and Execution**:
  ```bash
  pip3 install pycoQC

  pycoQC -f ./pycoQC/anac017-1_C_R1_sequencing_summary.txt -a ./pycoQC/anac017-1_C_R1_sorted.bam -o anac017-1_C_R1.html --report_title "anac017-1_C_R1"

### 3. Transcriptome Construction with FLAIR
- **Description**: Constructs and aligns transcriptomes using FLAIR. This step is optional and user-configurable.
- **Installation and Execution**:

  ```bash
  conda env create -n FLAIR_original -f FLAIR.yml
  cd Desktop/flair_master
  python flair.py align -g TAIR10_chr_all.fas -r WT_C_R1.fq.gz
  samtools sort flair.aligned.unsorted.bam -o flair.aligned.sorted.bam
  samtools index flair.aligned.sorted.bam
  python ./bin/bam2Bed12.py -i flair.aligned.sorted.bam > flair.aligned.sorted.bed
- **Reference Genome Normalization**:

  ```bash
  java -jar picard.jar NormalizeFasta I=TAIR_chr_all.fas O=TAIR_chr_all_norm.fas
- **Concatenation and Analysis**:

  ```bash
  cat *.bed > WT_flair_all_corrected_concatenated.bed
  python flair.py collapse -g TAIR10_chr_all.fas -r WT_concatenated.fastq -q WT_flair_all_corrected_concatenated.bed -f AtRTD2_19April2016.gtf
  python flair.py quantify -r reads_manifest.tsv -i flair.collapse.isoforms.fa --salmon --tpm
  python flair.py diffExp -q counts_matrix.tsv -o ./output/diffExp_salmon
  python flair.py diffSplice -i ./DRS/flair.collapse.isoforms.bed -q ./DRS/counts_matrix.tsv --test
  python ./bin/diff_iso_usage.py ./DRS/WT_counts_matrix_sumadas.tsv WT_C WT_AA diff_isos.txt
### 4. Additional Analysis
- **Description**: In-depth analysis of the results, including multivariable analysis.
- **Objective**: Identify significant isoforms and perform comprehensive analysis.

- **Proposed Commands**:

- **Result Analysis**:

  ```python

  # Evaluate differences with DEXSeq and DRIMSeq
  # References:
  # https://towardsdatascience.com/a-large-sample-crisis-or-not-640224020757
  # https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4
- ** Multivariable Analysis (PCA, sPLS)**:

  ```python
    # Perform PCA, sPLS to determine isoforms with specific characteristics.
    # Count the number of new isoforms, DEGs, classify splicing events, run IsoformSwitchAnalyzeR

## ⚙ Configuration

Configuration files are located in the config directory. Key files include:

- **config.yaml**: The main configuration file where you specify paths, parameters, and options for the analysis.
- **environment.yml**: The Conda environment configuration file with required dependencies.
- **Snakefile**: Defines the workflow rules and dependencies.

## 🏗️ Pipeline Structure

The pipeline is implemented using Snakemake, which manages the workflow and dependencies. The structure includes:

- **Basecalling**: Converts raw sequencing signals to DNA/RNA sequences.
- **Normalization**: Prepares the reference genome for mapping.
- **Mapping**: Aligns reads to the reference genome.
- **Modification Detection**: Identifies epitranscriptomic modifications.
- **Isoform Analysis (optional)**: Constructs and quantifies isoforms.
- **Quality Control**: Evaluates data quality and generates metrics.
- **Multivariable Analysis**: Performs PCA and SPLS for data interpretation.

## 📊 Examples and Use Cases

Explore the following use cases to see how K-CHOPORE can be applied:

- **Case Study 1**: Analysis of m6A modifications in RNA samples.
- **Case Study 2**: Isoform quantification and differential expression analysis.

Detailed examples and results are available in the results directory.

## 🛠️ Troubleshooting

If you encounter issues, check the following:

- Ensure all dependencies are correctly installed.
- Verify paths and file names in configuration files.
- Consult the logs in the logs directory for detailed error messages.
- For further assistance, refer to the FAQ or open an issue on GitHub.

## 🧬 Citations and References

If you use K-CHOPORE in your research, please cite the tools and references listed below to acknowledge the developers' contributions:

**Tools Integrated in K-CHOPORE**:

- **Guppy**: The official basecaller from Oxford Nanopore Technologies.
  - Citation: Oxford Nanopore Technologies. "Guppy basecaller software."
- **pycoQC**: A tool for quality control of Oxford Nanopore sequencing data.
  - Citation: Leger, A., & Leonardi, T. (2019). pycoQC, interactive quality control for Oxford Nanopore sequencing. Bioinformatics, 35(23), 5243-5245. doi:10.1093/bioinformatics/btz535
- **FLAIR**: A tool for the correction and quantification of isoforms from direct RNA sequencing.
  - Citation: Tang, A. D., Soulette, C. M., van Baren, M. J., Hart, K., Hrabeta-Robinson, E., Wu, C. J., & Brooks, A. N. (2020). Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. Nature Communications, 11(1), 1438. doi:10.1038/s41467-020-15213-2
- **Picard**: A set of command-line tools for manipulating high-throughput sequencing (HTS) data.
  - Citation: Broad Institute. "Picard toolkit." Broad Institute, GitHub repository.
- **IsoformSwitchAnalyzeR**: A tool for detecting and visualizing isoform switches.
  - Citation: Vitting-Seerup, K., & Sandelin, A. (2017). The landscape of isoform switches in human cancers. Molecular Cancer Research, 15(9), 1206-1220. doi:10.1158/1541-7786.MCR-16-0459


## 🤝 Contributing

We welcome contributions from the community! To contribute:

1. Fork the repository.
2. Create a feature branch.
3. Make your changes and test them.
4. Submit a pull request with a detailed description of your changes.
5. Pray for a second
Please follow the guidelines outlined in our CONTRIBUTING.md.

## 📞 Contact and Support

For questions or support:

- **Email**: pelayo.gonzalez@ispasturias.es
- **GitHub Issues**: Submit an issue
- **Community Forum**: Link to forum

## 🙏 Acknowledgments

We thank the Cancer Epigenetics and Nanomedicine Lab and the Systems Biology Lab at the University of Oviedo for their support. 

## 📂 Directory Structure

The directory structure of the K-CHOPORE project is designed to keep all resources organized and easily accessible:

  ```bash
├── Dockerfile
├── envir.yml
├── setup_directories.sh
├── src
│   ├── README.md                  # Overview of the source code and pipeline components.
│   ├── config                     # Configuration files and settings.
│   │   └── README.md              # Detailed explanation of each configuration file.
│   ├── data                       # Data storage directory.
│   │   ├── processed              # Processed data outputs.
│   │   │   ├── README.md          # Description of processed data.
│   │   │   ├── alignments         # BAM and SAM files post-alignment.
│   │   │   ├── counts             # Count matrices and related files.
│   │   │   └── quality_reports    # Quality control reports.
│   │   ├── raw                    # Raw input data.
│   │   │   ├── README.md          # Description of raw data and sources.
│   │   │   ├── bam                # Raw BAM files.
│   │   │   ├── fastq              # Raw FASTQ files.
│   │   │   └── metadata           # Metadata files for the raw data.
│   │   └── reference              # Reference genome and annotations.
│   │       ├── README.md          # Description of reference genome data.
│   │       ├── annotations        # GTF/GFF and related files.
│   │       └── genome             # Reference genome sequences.
│   ├── docs                       # Documentation for the project.
│   │   └── README.md              # Main documentation overview.
│   ├── envs                       # Environment management.
│   │   └── README.md              # Details on environment setup and configurations.
│   ├── logs                       # Log files generated during execution.
│   │   └── README.md              # Explanation of log files and troubleshooting.
│   ├── notebooks                  # Jupyter notebooks for data analysis.
│   │   ├── README.md              # Overview of the analysis notebooks.
│   │   └── main.txt               # Primary notebook or notebook-related notes.
│   ├── pipelines                  # Pipeline scripts and workflow definitions.
│   │   └── pipeline.sh            # Main pipeline execution script.
│   ├── publication                # Files related to publication.
│   │   ├── README.md              # Instructions and details for publication.
│   │   ├── figures                # Figures for the publication.
│   │   │   └── README.md          # Description of figures and sources.
│   │   ├── manuscript             # Manuscript drafts and revisions.
│   │   │   ├── drafts             # Initial drafts.
│   │   │   │   └── README.md      # Details about the drafts.
│   │   │   ├── revisions          # Revisions based on feedback.
│   │   │   │   └── README.md      # Details about the revisions.
│   │   │   └── submission         # Final submission files.
│   │   │       └── README.md      # Submission instructions and final files.
│   │   ├── supplementary          # Supplementary data and materials.
│   │   │   └── README.md          # Description of supplementary materials.
│   │   └── tables                 # Data tables for the publication.
│   │       └── README.md          # Explanation of tables and data sources.
│   ├── results                    # Results from the analysis.
│   │   ├── README.md              # Overview of the results directory.
│   │   ├── figures                # Generated figures from analysis.
│   │   ├── reports                # Analysis reports and summaries.
│   │   └── tables                 # Data tables generated from analysis.
│   ├── reviews                    # Peer review responses and revisions.
│   │   ├── README.md              # Overview of the review process.
│   │   ├── response               # Responses to reviewer comments.
│   │   ├── review1                # Review 1 specific documents.
│   │   │   └── README.md          # Details about Review 1 feedback.
│   │   └── review2                # Review 2 specific documents.
│   │       └── README.md          # Details about Review 2 feedback.
│   └── scripts                    # Custom scripts for various tasks.
│       ├── README.md              # Overview of available scripts.
│       ├── alignment              # Scripts related to alignment tasks.
│       ├── analysis               # Scripts for data analysis.
│       ├── postprocessing         # Post-processing scripts.
│       ├── preprocessing          # Pre-processing scripts.
│       ├── quality_control        # Quality control scripts.
│       └── utilities              # Utility scripts.
└── start_kchopo.sh                # Script to start the K-CHOPORE pipeline.
