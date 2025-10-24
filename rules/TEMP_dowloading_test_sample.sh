#!/bin/bash

sudo apt install -y sra-toolkit

SRR_FILE=SRR10510434

# Output folder
OUTDIR=~/data/mydatalocal
#OUTDIR=~/data/mydatalocal


# Number of threads to use
THREADS=4

# Loop on each SRR
echo "Downloading $SRR_FILE"
    
# Delete the existing file if there is one
if [ -f "$OUTDIR/$SRR_FILE.fastq" ]; then
    echo "Existing file detected: deleting $OUTDIR/$SRR_FILE.fastq"
    rm "$OUTDIR/$SRR_FILE.fastq"
fi

# Start download
fasterq-dump --threads $THREADS --progress --outdir "$OUTDIR" "$SRR_FILE"

echo "Fastq file download complete"