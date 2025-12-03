#!/usr/bin/env bash
set -e

# Absolute path to the project
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$PROJECT_DIR"

# Create logs directory if missing
mkdir -p logs

# Run the workflow
snakemake \
    --use-singularity \
    --cores all \
    --printshellcmds \
    --directory "$PROJECT_DIR"
