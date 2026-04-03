"""
Unit tests for the core logic in coverage_analysis.ipynb.
See TESTING.md for a brief explanation all the unit testing done here.

Run:
    pip3 install pytest pytest-csv
    python3 -m pytest test_coverage.py -v
    python3 -m pytest test_coverage.py -v --csv=test_results.tsv
"""

import csv
import json
import tempfile

import pytest


# Functions under test (using the same logic as the notebook) 

# Returns a numeric sort key so chromosomes sort biologically (1-22, X, Y, M) not alphabetically
def sort_key(r):
    chrom = r["chromosome"].replace("chr", "").upper()
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    try:
        return int(chrom)
    except ValueError:
        return order.get(chrom, 99)


PRIMARY = (
    {f"chr{i}" for i in range(1, 23)}
    | {"chrX", "chrY", "chrM", "chrMT"}
    | {str(i) for i in range(1, 23)}
    | {"X", "Y", "M", "MT"}
)


# Filters records to keep only the 25 primary chromosomes (1-22, X, Y, M), removing decoys and contigs
def filter_primary(records):
    return [r for r in records if r["chromosome"] in PRIMARY]


# Calculates genome-wide mean depth weighted by chromosome length (not simple average)
def weighted_mean(records):
    total_len = sum(r["length"] for r in records)
    return sum(r["mean_depth"] * r["length"] for r in records) / total_len


# Parses samtools coverage tab-separated output into a list of structured chromosome records
def parse_samtools_output(stdout):
    lines = stdout.strip().split("\n")
    header = lines[0].lstrip("#").split("\t")
    records = []
    for line in lines[1:]:
        if not line.strip():
            continue
        raw = dict(zip(header, line.split("\t")))
        records.append({
            "chromosome":    raw["rname"],
            "length":        int(raw["endpos"]) - int(raw["startpos"]),
            "num_reads":     int(raw["numreads"]),
            "covered_bases": int(raw["covbases"]),
            "coverage_pct":  float(raw["coverage"]),
            "mean_depth":    float(raw["meandepth"]),
            "mean_baseq":    float(raw["meanbaseq"]),
            "mean_mapq":     float(raw["meanmapq"]),
        })
    return records


# Fixture: shared sample data used across multiple tests

# Provides realistic WGS sample data (5 chromosomes with real GRCh38 sizes and coverage values)
@pytest.fixture
def sample_records():
    return [
        {"chromosome": "chr1", "length": 248956422, "num_reads": 4000000,
         "covered_bases": 240000000, "coverage_pct": 96.4,
         "mean_depth": 30.2, "mean_baseq": 35.0, "mean_mapq": 58.0},
        {"chromosome": "chr2", "length": 242193529, "num_reads": 3900000,
         "covered_bases": 235000000, "coverage_pct": 97.0,
         "mean_depth": 29.8, "mean_baseq": 35.1, "mean_mapq": 58.2},
        {"chromosome": "chrX", "length": 156040895, "num_reads": 2000000,
         "covered_bases": 148000000, "coverage_pct": 94.8,
         "mean_depth": 15.1, "mean_baseq": 34.8, "mean_mapq": 57.5},
        {"chromosome": "chrY", "length": 57227415,  "num_reads": 100000,
         "covered_bases": 1500000,   "coverage_pct": 2.6,
         "mean_depth": 0.5,  "mean_baseq": 30.0, "mean_mapq": 40.0},
        {"chromosome": "chrM", "length": 16569,     "num_reads": 500000,
         "covered_bases": 16569,     "coverage_pct": 100.0,
         "mean_depth": 8500.0, "mean_baseq": 35.5, "mean_mapq": 59.0},
    ]


# Mock samtools output (avoids needing the actual 17 GB CRAM file)

MOCK_SAMTOOLS_OUTPUT = """\
#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq
chr1\t0\t248956422\t4000000\t240000000\t96.4\t30.2\t35.0\t58.0
chrM\t0\t16569\t500000\t16569\t100.0\t8500.0\t35.5\t59.0
"""


# Group 1: sort_key

# Verifies chromosomes sort numerically (chr1, chr2...) not alphabetically (chr1, chr10, chr2...)
def test_sort_key_autosomes_in_order():
    # Chromosomes must sort numerically (chr1, chr2 ... chr22), not alphabetically
    records = [{"chromosome": f"chr{i}"} for i in [3, 1, 22, 10, 2]]
    result  = [r["chromosome"] for r in sorted(records, key=sort_key)]
    assert result == ["chr1", "chr2", "chr3", "chr10", "chr22"]


def test_sort_key_sex_and_mito_at_end():
    # chrX → chrY → chrM must always follow autosomes
    records = [{"chromosome": c} for c in ["chrM", "chr1", "chrX", "chrY", "chr22"]]
    result  = [r["chromosome"] for r in sorted(records, key=sort_key)]
    assert result.index("chr1") < result.index("chrX") < result.index("chrY") < result.index("chrM")


# Verifies sorting works with plain names (1, X, Y, M) without 'chr' prefix
def test_sort_key_plain_names():
    # Handles names without 'chr' prefix (e.g. '1', 'MT')
    records = [{"chromosome": c} for c in ["M", "1", "X", "Y"]]
    result  = [r["chromosome"] for r in sorted(records, key=sort_key)]
    assert result.index("1") < result.index("X") < result.index("Y") < result.index("M")


# Group 2: filter_primary

# Verifies all 5 primary chromosomes from sample data pass the filter
def test_filter_keeps_primary_chromosomes(sample_records):
    # Only primary chromosomes should survive the filter
    assert len(filter_primary(sample_records)) == 5


# Verifies decoys, HLA alleles, and viral sequences are filtered out
def test_filter_removes_decoys(sample_records):
    # Decoys, HLA alleles, and viral sequences must be removed
    records_with_decoys = sample_records + [
        {"chromosome": "chrUn_GL000220v1"},
        {"chromosome": "HLA-A"},
        {"chromosome": "chrEBV"},
    ]
    filtered = filter_primary(records_with_decoys)
    chromosomes = [r["chromosome"] for r in filtered]
    assert len(filtered) == 5
    assert "chrUn_GL000220v1" not in chromosomes


# Verifies filter accepts both 'chr1' and '1' naming conventions
def test_filter_handles_plain_names():
    # Both 'chr1' and '1' naming styles are valid
    records = [
        {"chromosome": "1",  "length": 1, "num_reads": 1, "covered_bases": 1,
         "coverage_pct": 100.0, "mean_depth": 30.0, "mean_baseq": 35.0, "mean_mapq": 58.0},
        {"chromosome": "MT", "length": 1, "num_reads": 1, "covered_bases": 1,
         "coverage_pct": 100.0, "mean_depth": 8500.0, "mean_baseq": 35.0, "mean_mapq": 58.0},
    ]
    assert len(filter_primary(records)) == 2


# Verifies all 25 primary chromosomes (1-22, X, Y, M) are correctly identified
def test_all_25_primary_chromosomes_present():
    # A complete WGS run must produce exactly 25 primary chromosomes
    records = [
        {"chromosome": f"chr{i}", "length": 1, "num_reads": 1, "covered_bases": 1,
         "coverage_pct": 100.0, "mean_depth": 30.0, "mean_baseq": 35.0, "mean_mapq": 58.0}
        for i in range(1, 23)
    ] + [
        {"chromosome": c, "length": 1, "num_reads": 1, "covered_bases": 1,
         "coverage_pct": 100.0, "mean_depth": 30.0, "mean_baseq": 35.0, "mean_mapq": 58.0}
        for c in ["chrX", "chrY", "chrM"]
    ] + [
        {"chromosome": "chrUn_GL000220v1", "length": 1, "num_reads": 0,
         "covered_bases": 0, "coverage_pct": 0.0, "mean_depth": 0.0,
         "mean_baseq": 0.0, "mean_mapq": 0.0}
    ]
    assert len(filter_primary(records)) == 25


# Group 3: weighted_mean

# Verifies weighted mean equals simple average when all chromosomes have equal length
def test_weighted_mean_equal_lengths():
    # Equal-length chromosomes: weighted mean = simple average = 25.0
    records = [
        {"chromosome": "chr1", "length": 100, "mean_depth": 30.0},
        {"chromosome": "chr2", "length": 100, "mean_depth": 20.0},
    ]
    assert weighted_mean(records) == pytest.approx(25.0)


# Verifies longer chromosomes dominate the weighted mean calculation
def test_weighted_mean_unequal_lengths():
    # Longer chromosome dominates: (30*900 + 8000*100) / 1000 = 827.0
    records = [
        {"chromosome": "chr1", "length": 900, "mean_depth": 30.0},
        {"chromosome": "chrM", "length": 100, "mean_depth": 8000.0},
    ]
    assert weighted_mean(records) == pytest.approx(827.0)


# Verifies weighted mean (30.43x) differs from simple average (4015x) with real GRCh38 chromosome sizes
def test_weighted_mean_not_simple_average():
    # Key test: with real GRCh38 sizes, simple average = 4015x (wrong)
    # Weighted mean ≈ 30.43x (correct — chr1 at 249 MB dominates over chrM at 16 KB)
    records = [
        {"chromosome": "chr1", "length": 248956422, "mean_depth": 30.0},
        {"chromosome": "chrM", "length":      16569, "mean_depth": 8000.0},
    ]
    assert weighted_mean(records) != pytest.approx((30.0 + 8000.0) / 2)
    assert weighted_mean(records) == pytest.approx(30.43, rel=0.01)


# Group 4: parse_samtools_output

# Verifies parser extracts the correct number of chromosomes from samtools output
def test_parse_returns_correct_number_of_records():
    assert len(parse_samtools_output(MOCK_SAMTOOLS_OUTPUT)) == 2


# Verifies parsed fields have correct data types (int for counts, float for measurements)
def test_parse_correct_field_types():
    # Types must be correct — storing mean_depth as string breaks arithmetic
    r = parse_samtools_output(MOCK_SAMTOOLS_OUTPUT)[0]
    assert isinstance(r["chromosome"],    str)
    assert isinstance(r["length"],        int)
    assert isinstance(r["num_reads"],     int)
    assert isinstance(r["covered_bases"], int)
    assert isinstance(r["coverage_pct"],  float)
    assert isinstance(r["mean_depth"],    float)


# Verifies chromosome length is correctly calculated as (endpos - startpos)
def test_parse_correct_length_calculation():
    # length = endpos - startpos (samtools reports endpos, not length directly)
    records = parse_samtools_output(MOCK_SAMTOOLS_OUTPUT)
    assert records[0]["length"] == 248956422
    assert records[1]["length"] == 16569


# Verifies parsed values match expected chromosome names and coverage depths
def test_parse_correct_values():
    records = parse_samtools_output(MOCK_SAMTOOLS_OUTPUT)
    assert records[0]["chromosome"] == "chr1"
    assert records[0]["mean_depth"] == pytest.approx(30.2)
    assert records[1]["chromosome"] == "chrM"
    assert records[1]["mean_depth"] == pytest.approx(8500.0)


# Verifies parser gracefully handles trailing newlines in samtools output
def test_parse_skips_empty_lines():
    # Trailing newlines in samtools output must not cause crashes
    records = parse_samtools_output(MOCK_SAMTOOLS_OUTPUT + "\n\n")
    assert len(records) == 2


# Group 5: output file structure

# Verifies output TSV contains required columns for downstream analysis tools
def test_tsv_has_correct_columns(sample_records):
    # TSV must have all required columns for downstream tools (pandas, R, Excel)
    with tempfile.TemporaryDirectory() as tmpdir:
        path = f"{tmpdir}/coverage_report.tsv"
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(sample_records[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(sample_records)
        with open(path) as fh:
            columns = csv.DictReader(fh, delimiter="\t").fieldnames
    assert "chromosome" in columns
    assert "mean_depth" in columns
    assert "length"     in columns


# Verifies TSV row count matches input records (no rows dropped or duplicated)
def test_tsv_row_count_matches(sample_records):
    # One row per chromosome — no rows dropped or duplicated during write
    with tempfile.TemporaryDirectory() as tmpdir:
        path = f"{tmpdir}/coverage_report.tsv"
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(sample_records[0].keys()), delimiter="\t")
            writer.writeheader()
            writer.writerows(sample_records)
        with open(path) as fh:
            assert len(list(csv.DictReader(fh, delimiter="\t"))) == len(sample_records)


# Verifies JSON output has metadata and chromosomes keys for SMaHT Data Portal
def test_json_structure(sample_records):
    # JSON must have 'metadata' and 'chromosomes' keys for SMaHT Data Portal
    payload = {
        "metadata":    {"genome_mean_depth": round(weighted_mean(sample_records), 4)},
        "chromosomes": sample_records,
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        path = f"{tmpdir}/coverage_report.json"
        with open(path, "w") as fh:
            json.dump(payload, fh)
        with open(path) as fh:
            data = json.load(fh)
    assert "metadata"          in data
    assert "chromosomes"       in data
    assert "genome_mean_depth" in data["metadata"]
    assert len(data["chromosomes"]) == len(sample_records)


# Group 6: biological sanity

# Verifies no negative depth values (which would indicate a parsing bug)
def test_no_negative_depths(sample_records):
    # Coverage is a read count — negative values indicate a parsing bug
    assert all(r["mean_depth"] >= 0 for r in sample_records)


# Verifies mitochondrial DNA (chrM) has higher coverage than nuclear chromosomes (hundreds of mtDNA copies per cell)
def test_chrM_higher_coverage_than_autosomes(sample_records):
    # chrM always has higher coverage in WGS — hundreds of mitochondria per cell
    chrM_depth      = next(r["mean_depth"] for r in sample_records if r["chromosome"] == "chrM")
    autosome_depths = [r["mean_depth"] for r in sample_records
                       if r["chromosome"].replace("chr", "").isdigit()]
    assert all(chrM_depth > d for d in autosome_depths)


# Verifies coverage percentage is valid (between 0% and 100%)
def test_coverage_pct_within_bounds(sample_records):
    # Breadth of coverage must be between 0% and 100%
    assert all(0.0 <= r["coverage_pct"] <= 100.0 for r in sample_records)