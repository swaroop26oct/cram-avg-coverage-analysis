# CRAM Coverage Analysis — COLO829T

Calculate average per-base sequence coverage for each chromosome in a CRAM file.

## How to Run

SSH into an Ubuntu 22.04 EC2 instance (us-east-2, 30 GB disk), upload `ec2_setup.sh` and `coverage_analysis.ipynb`, then:

```bash
bash ec2_setup.sh
```

Access Jupyter via SSH tunnel: `ssh -L 8888:localhost:8888 -i key.pem ubuntu@<IP>`

Open the printed URL in your browser and run all cells in `coverage_analysis.ipynb`.

## Data

Files must be in the working directory before running the notebook:

```bash
wget https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/COLO829T_TEST.cram
wget https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/COLO829T_TEST.cram.crai
wget https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
```

## Method

`samtools coverage` runs a pileup over every position in the CRAM and reports the true mean depth per chromosome (positions with zero depth included in the denominator). This is exact — unlike the approximation `(reads x read_length) / chr_length`, which overcounts clipped bases.

## Output

`coverage_report.tsv` - This file can be used by Bioinformaticians for downstream analysis for easy uploading to pandas, R and other downstream scripts
`coverage_report.json`- This file  consists of JSON with metadata and can be used by Software engineers in Data Portal and REST API calls.
`coverage_report.png`- This is a Log-scale bar chart, can be used in publications and presentations.

## Unit Tests

21 tests covering sorting, filtering, weighted mean calculation, samtools output parsing, output file structure, and biological sanity checks.

```bash
pip3 install pytest pytest-csv
python3 -m pytest test_coverage.py -v --csv=test_results.tsv --csv-columns=name,status
```

See `TESTING.md` for a full explanation of what each test checks and why.

## Requirements

- samtools >= 1.12
- Python >= 3.8
- matplotlib

---

> **Note:** A `Dockerfile` is included for reproducibility. To run in a containerised environment instead: `docker build -t cram-coverage:1.0 . && docker run -p 8888:8888 -v $(pwd):/data cram-coverage:1.0`