# Unit Testing Guide — coverage_analysis

## To run

We are using **pytest** — It automatically finds all functions starting with `test_` and runs them. A green dot means passed, `F` means failed.


```bash
# Install dependencies (one time only)
pip3 install pytest pytest-csv

# Run tests with verbose output and save results as TSV
python3 -m pytest test_coverage.py -v --csv=test_results.tsv --csv-columns=name,status
```

## What the TSV output looks like

The `--csv` flag (from the `pytest-csv` plugin) produces a `test_results.tsv` file with one row per test:

| name | status |
|------|--------|
| test_sort_key_autosomes_in_order | passed |
| test_chrM_higher_coverage_than_autosomes | passed |

This is useful for sharing test results with the team or including in CI/CD pipelines.

## What each test group checks

### Group 1: sort_key (3 tests)
Chromosomes must sort in biological order (chr1→chr2→...→chrX→chrY→chrM), not alphabetically (which would give chr1, chr10, chr11, chr2...). If this fails, the plot and TSV will have chromosomes in the wrong order.

### Group 2: filter_primary (4 tests)
GRCh38 has 195 contigs — we only want 25 primary ones. Tests verify decoys (`chrUn_GL000220v1`), HLA alleles, and viral sequences are removed, and that both `chr`-prefixed and plain chromosome names work.

### Group 3: weighted_mean (3 tests)
The most important calculation. The key test (`test_weighted_mean_not_simple_average`) proves that a simple average would give **4015x** while the correct weighted answer is **~30.43x** — because chrM (16 KB) is so tiny compared to chr1 (249 MB) that it barely affects the genome-wide mean.

### Group 4: parse_samtools_output (5 tests)
samtools returns plain text. Tests verify our parser produces correct Python types (int/float/str). If `mean_depth` is stored as a string instead of float, arithmetic like `mean_depth * length` silently crashes.

### Group 5: output file structure (3 tests)
Uses `tempfile.TemporaryDirectory` so test files are auto-deleted after each test — no leftover files on disk. Verifies TSV has the right columns and JSON has both `metadata` and `chromosomes` sections.

### Group 6: biological sanity (3 tests)
These don't test code logic — they test that the biology makes sense. If `chrM_higher_coverage_than_autosomes` fails on real data, something is wrong with the input file (sample swap, missing chrM in reference, very low-quality data), not the code.

## Definitions

### `@pytest.fixture`
A **fixture** is used as a reusable piece of test setup. Instead of copy-pasting the same sample data into every test, we define it once as `sample_records` and pytest automatically passes it to any test that lists it as a parameter. 
