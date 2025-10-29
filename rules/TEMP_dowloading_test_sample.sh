#sudo apt install -y sra-toolkit

# Get parameters from config
THREADS = config["THREADS"]

# Output folder
OUTDIR=~/data/mydatalocal/output/fastq

# SRR file to download
SRR_FILE=SRR10510434

rule download_fastq:

    output:
        "$OUT_DIR/$SRR_FILE.fastq",

    #container:

    threads:THREADS

    shell:

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