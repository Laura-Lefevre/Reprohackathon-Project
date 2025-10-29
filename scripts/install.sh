sudo apt update
sudo apt install -y snakemake

# Install SRA Toolkit for downloading FASTQ files
sudo apt install -y sra-toolkit

# Install Cutadapt v1.11 for trimming FASTQ files
pip install --user cutadapt==1.11
export PATH=$HOME/.local/bin:$PATH