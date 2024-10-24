# K-CHOPORE 🌟
### Keen Comprehensive High-throughput Omics Pipeline Organizer

Just like the iconic Asturian **cachopo**, K-CHOPORE is a layered and satisfying bioinformatics pipeline designed for analyzing **Nanopore RNA-seq** data with a focus on **epitranscriptomics**! 🧬🍴 Dive into the data feast, where each component works together to deliver high-quality, reproducible results. 📊💻

---

## 📜 Overview
**K-CHOPORE** is an open-source pipeline for the comprehensive analysis of Nanopore sequencing data, tailored to handle every step from basecalling to epitranscriptomic modification detection. It integrates multiple cutting-edge bioinformatics tools, including **Snakemake**, **Docker**, **Python**, and well-established tools such as **Guppy**, **Minimap2**, **FLAIR**, and **ELIGOS2**. The pipeline emphasizes **FAIR-compliance** (Findable, Accessible, Interoperable, and Reusable) to ensure reproducibility and scalability across diverse research settings.

### Key Features
- **FAIR-Compliant**: Adheres to the FAIR principles, ensuring Findable, Accessible, Interoperable, and Reusable workflows.
- **Automated Workflow**: Uses Snakemake to streamline each step of the analysis, reducing manual intervention and errors.
- **Containerized Environment**: Employs Docker to guarantee reproducibility across different systems, solving dependency issues.
- **Comprehensive Tools**: Integrates well-known tools such as **Guppy** (basecalling), **Minimap2** (alignment), **FLAIR** (isoform quantification), and **ELIGOS2** (RNA modification detection).

---
## Table of Contents
- [Key Features](#key-features)
- [Getting Started](#getting-started)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Pipeline Structure](#pipeline-structure)
- [Running the Pipeline](#running-the-pipeline)
- [Customization](#customization)
- [Contributing](#contributing)
- [License](#license)
- [Hands-On Guide](#hands-on-guide)
- [Troubleshooting](#troubleshooting)

## 📖 Getting Started
To get started with **K-CHOPORE**, you'll need to meet a few prerequisites and follow the installation steps below. The pipeline is designed to be **modular** and **flexible**, allowing you to analyze **RNA-seq data** from plants or any other organisms, performing everything from **basecalling** to **epitranscriptomic modification detection**.

## 🔧 Prerequisites
Ensure you have the following installed:

### 1. Docker
Docker is required to run **K-CHOPORE** in a fully isolated environment. Install Docker following the instructions for your system:
- [Docker for Windows](https://docs.docker.com/docker-for-windows/install/)
- [Docker for Mac](https://docs.docker.com/docker-for-mac/install/)
- [Docker for Linux](https://docs.docker.com/engine/install/)

Verify that Docker is installed:

```bash
docker --version
```
### 2. Snakemake
#### Install Snakemake for workflow management:

```bash
pip install snakemake
```
#### Check that Snakemake is installed correctly:

```bash
snakemake --version
```

### 3. Git
#### Install Git to clone the K-CHOPORE repository:


```bash
git --version
```

### 4. Python 3.8+
#### K-CHOPORE requires Python 3.8+ to run custom scripts. Ensure that the correct version of Python is installed:

```bash
python --version
```

## 🔥 Installation
Follow these steps to install **K-CHOPORE** and set up the necessary environment:

### Step 1: Clone the K-CHOPORE Repository
Open a terminal and clone the **K-CHOPORE** GitHub repository:

```bash
git clone https://github.com/pelayovic/K-CHOPORE.git
cd K-CHOPORE
```
### Step 2: Build the Docker Image
Inside the cloned K-CHOPORE directory, build the Docker image:

```bash
sudo docker build -t k-chopore .
```
This process will pull in all necessary dependencies, ensuring a reproducible environment.

### Step 3: Configure the Pipeline
Edit the config.yml file located in the config/ folder to point to your input files (e.g., FAST5 files for Nanopore sequencing) and the reference genome:

```yaml
input_files:
  fast5_dir: "/workspace/data/raw/fast5"
  reference_genome: "/workspace/data/reference/genome.fasta"
output_dir: "/workspace/results"
```
## 🏗️ Pipeline Structure
K-CHOPORE’s structure follows a modular approach, ensuring clarity and scalability. The core components are as follows:

```bash

K-CHOPORE/
│
├── data/                      # Input data (FAST5, FASTQ, etc.)
├── results/                   # Pipeline outputs
├── config/                    # Configuration files for the pipeline
├── scripts/                   # Python scripts used in the pipeline
├── Snakefile                  # Snakemake workflow definition
├── Dockerfile                 # Dockerfile for building the containerized environment
└── README.md                  # Project documentation
```

## 🚀 Running the Pipeline
#### Basic Command
To execute the full K-CHOPORE pipeline, run the following command. This will start with basecalling and proceed through RNA modification detection:

```bash

sudo docker run -it --rm -v /path/to/your/local/data:/workspace k-chopore \
    snakemake --snakefile /workspace/Snakefile --configfile /workspace/config/config.yml \
    --cores 10 --latency-wait 30 --printshellcmds
```
#### Execute a Specific Rule
If you want to run only a specific rule (e.g., basecalling or alignment), use:

```bash
sudo docker run -it --rm -v /path/to/your/local/data:/workspace k-chopore \
    snakemake --snakefile /workspace/Snakefile --cores 10 basecalling
```

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

## Tools You Can Add to K-CHOPORE for RNA-Seq in Plants
K-CHOPORE can be further enhanced with additional tools designed for RNA-seq in plants:

**Rattle**: For isoform clustering and annotation in Nanopore RNA-seq data.
**IsoSeq3**: Alternative analysis for detecting long isoforms, especially useful in plant transcriptomics.
**NanoComp**: For comparative analysis of sequencing quality between samples or conditions.
**NanoPlot**: Provides visual summaries of read length, quality, and distribution.
**SQANTI3**: Validates and classifies isoforms generated by long-read technologies.
**TALON**: Identifies and quantifies novel transcripts in plant transcriptomics.
**DESeq2**: Key for differential expression analysis, useful in transcriptomics projects.
**NanoMod**: Complements ELIGOS2 by detecting base modifications via current signal differences in Nanopore data.


## Limitations of K-CHOPORE
### Basecalling Bottlenecks
Basecalling using ONT’s Guppy can be a bottleneck, especially when data access from FAST5 files outpaces the processing power of the CPU/GPU. This is a potential friction point in any Nanopore workflow (Amarasinghe et al. 2020).


## 🙏 Acknowledgments

We thank the Cancer Epigenetics and Nanomedicine Lab and the Systems Biology Lab at the University of Oviedo for their support. 

