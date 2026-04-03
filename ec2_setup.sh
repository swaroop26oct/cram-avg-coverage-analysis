#!/usr/bin/env bash
# =============================================================================
# ec2_setup.sh — Run once after SSHing into a fresh Ubuntu 22.04 EC2 instance.
# Installs Jupyter + samtools, downloads data, starts Jupyter.
# =============================================================================
set -euo pipefail

echo "==== Installing dependencies ===="
sudo apt-get update -qq
sudo apt-get install -y samtools python3-pip wget
pip3 install --quiet jupyter matplotlib

echo ""
echo "==== Downloading data (~21 GB) ===="
BASE="https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW"
for FILE in COLO829T_TEST.cram COLO829T_TEST.cram.crai GCA_000001405.15_GRCh38_no_alt_analysis_set.fa; do
    if [ -f "$FILE" ]; then
        echo "  Already exists: $FILE"
    else
        echo "  Downloading $FILE ..."
        wget -q --show-progress -c "${BASE}/${FILE}"
    fi
done

echo ""
echo "==== Indexing reference FASTA ===="
if [ ! -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai ]; then
    samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
fi
echo "  Done."

echo ""
echo "==== Starting Jupyter ===="
echo "  Access it from your Mac browser at: http://localhost:8888"
echo "  (requires SSH tunnel — see instructions)"
echo ""
jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser
