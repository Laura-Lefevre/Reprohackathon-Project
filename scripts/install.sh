sudo apt update
sudo apt install -y snakemake

# Install SRA Toolkit for downloading FASTQ files
sudo apt install -y sra-toolkit

# Install Cutadapt v1.11 for trimming FASTQ files
pip install --user cutadapt==1.11
export PATH=$HOME/.local/bin:$PATH

# sudo apt update
# $ sudo apt install -y wget
# $ cd /tmp
# $ wget https://github.com/apptainer/apptainer/releases/download/v1.4.3/apptainer_1.4.3_
# ˓→amd64.deb
# $ sudo apt install -y ./apptainer_1.4.3_amd64.deb
# $ rm ./apptainer_1.4.3_amd64.deb

