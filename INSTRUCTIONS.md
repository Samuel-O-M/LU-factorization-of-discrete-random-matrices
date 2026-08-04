# Instructions for Running LU Factorization of Discrete Random Matrices

This document provides detailed instructions for running all components of this project from the repository root.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Project Overview](#project-overview)
3. [Running Locally](#running-locally)
4. [Running on a Supercomputer (SLURM)](#running-on-a-supercomputer-slurm)
5. [Utility Scripts](#utility-scripts)
6. [Monte Carlo Experiments](#monte-carlo-experiments)
7. [Output Files](#output-files)

---

## Prerequisites

### Julia Installation

1. Install Julia 1.8 or later from [julialang.org](https://julialang.org/downloads/)
2. Ensure Julia is available in your PATH

### Package Dependencies

From the repository root, run:

```bash
julia -e 'using Pkg; Pkg.instantiate()'
```

This will install the required packages:
- `StaticArrays` - for efficient static arrays
- `Printf` - for formatted output
- `SortingAlgorithms` - for efficient radix sort

---

## Project Overview

The project computes LU factorization properties of discrete random matrices. It consists of two main computational pipelines:

| Pipeline | Matrix Size | Location | Description |
|----------|-------------|----------|-------------|
| **Local** | 8×8 | `src/local/` | Runs on local machine |
| **Supercomputer** | 9×9 | `src/super/` | Runs on SLURM cluster |
| **Monte Carlo** | 30+ | `experiment/` | Exact sampling with progressive expansion |

Each pipeline has two stages:
1. **Families** - Computes equivalence classes of matrices
2. **Aggregation** - Aggregates results by signature (row/column sums)

---

## Running Locally

The local pipeline computes results for 8×8 matrices. You can run these directly with Julia:

### Step 1: Compute Families

```bash
julia src/local/families.jl
```

Or with a specific number of threads:

```bash
JULIA_NUM_THREADS=8 julia src/local/families.jl
```

**What it does:**
- Computes equivalence classes (families) of 8×8 binary matrices
- Uses multithreading for parallel computation
- Outputs binary files to `data/local/binary/families_NxN.bin`

**Runtime:** Several minutes to hours depending on hardware.

### Step 2: Run Aggregation

```bash
julia src/local/aggregation.jl
```

**What it does:**
- Aggregates families by their signatures (row and column sums)
- Outputs counts to `data/local/binary/counts_8x8.bin`
- Generates summary to `results/summary_local.txt`

---

## Running on a Supercomputer (SLURM)

The supercomputer pipeline computes results for 9×9 matrices, which requires significantly more computational resources.

### Configuration

The SLURM scripts use environment variables for configuration. Default values are provided, but you can override them:

| Variable | Default | Description |
|----------|---------|-------------|
| `SBATCH_PARTITION` | `compute` | SLURM partition name |
| `SBATCH_MEM` | `200G` | Memory per node |
| `SBATCH_CPUS` | `100` | CPUs per task |
| `SBATCH_TIME` | `90-00:00:00` | Time limit (DD-HH:MM:SS) |
| `JULIA_MODULE` | *(none)* | Julia module to load (e.g., `julia/1.8.5`) |

### Step 1: Submit Families Job

Using defaults:

```bash
sbatch scripts/run_families.sh
```

With custom configuration:

```bash
SBATCH_PARTITION=gpu SBATCH_MEM=400G sbatch scripts/run_families.sh
```

If your cluster requires loading a Julia module:

```bash
JULIA_MODULE=julia/1.8.5 sbatch scripts/run_families.sh
```

**Output:**
- Job output: `families_<jobid>.out`
- Job errors: `families_<jobid>.err`
- Data files: `data/super/binary/families_NxN.bin`

### Step 2: Submit Aggregation Job

```bash
sbatch scripts/run_aggregation.sh
```

Same configuration options apply as for families.

**Output:**
- Job output: `aggregation_<jobid>.out`
- Job errors: `aggregation_<jobid>.err`
- Counts file: `data/super/binary/counts_9x9.bin`
- Summary: `results/summary_super.txt`

---

## Utility Scripts

All utility scripts are located in `src/utils/`.

### Analyze Ones Distribution

Computes the distribution of ones (non-zero entries) for n=1,...,9 from all available data.
Uses families files for n=1..7 and counts files for n=8,9:

```bash
julia src/utils/analyze_ones.jl
```

**Output:** `results/ones.json`

### Compute Strong Upper Bound

Computes upper bounds on probabilities for larger matrix dimensions.
Automatically processes both 8×8 and 9×9 data if available:

```bash
julia src/utils/compute_strong_upper_bound.jl       # Target N=30 (default)
julia src/utils/compute_strong_upper_bound.jl 50    # Target N=50
julia src/utils/compute_strong_upper_bound.jl 100   # Target N=100
```

**Parameters:**
- First argument: `N` - target matrix dimension (default: `30`)

**Output:**
- `results/strong_upper_bound_8x8.json`
- `results/strong_upper_bound_9x9.json`

### Convert Binary to JSON

Converts binary data files to human-readable JSON format:

```bash
julia src/utils/convert_to_json.jl local
julia src/utils/convert_to_json.jl super
```

**Output:** JSON files in `data/<mode>/readable/`

---

## Monte Carlo Experiments

Estimates P(strongly non‑singular at N) for Bernoulli(p) matrices using exact `BigInt` arithmetic. Samples from the exact 4×4 SNS database, then progressively extends to the target dimension. Each 4×4 matrix is weighted by its probability under Bernoulli(p).

### Step 1: Build the 4×4 SNS Database

Enumerates every 4×4 binary matrix with top‑left corner = 1 and collects those that are strongly non‑singular (all leading principal minors ≠ 0):

```bash
julia experiment/base_4x4.jl
```

**What it does:**
- Enumerates 32768 matrices
- Checks all leading principal minors (1×1, 2×2, 3×3, 4×4)
- Writes SNS matrices to `experiment/matrices_4x4.json`

### Step 2: Bernoulli(p) Probability Estimation

Loads the 4×4 SNS database and estimates P(SNS at N) for p = 0.01, 0.02, ..., 1.00 using progressive expansion through milestone sizes (6, 12, 20, 30, 45, 70):

```bash
julia experiment/ber(p).jl        # N = 30 (default)
julia experiment/ber(p).jl 50     # N = 50
julia experiment/ber(p).jl 100    # N = 100
```

**Parameters:**
- First argument: `N` — target matrix dimension (default: `30`)

**Output:** `experiment/ber_n<N>.json`

---

## Output Files

### Data Files

| File | Format | Description |
|------|--------|-------------|
| `families_NxN.bin` | Binary | Equivalence classes with counts |
| `counts_NxN.bin` | Binary | Aggregated counts by signature |

### Binary Format

**Families file:** Each record is 16 bytes:
- 8 bytes: Matrix representative (UInt64, packed bits)
- 8 bytes: Count (UInt64)

**Counts file (8×8):** Each record is 24 bytes:
- 8 bytes: Row sums (8 × UInt8)
- 8 bytes: Column sums (8 × UInt8)
- 8 bytes: Count (UInt64)

**Counts file (9×9):** Each record is 34 bytes:
- 9 bytes: Row sums (9 × UInt8)
- 9 bytes: Column sums (9 × UInt8)
- 16 bytes: Count (UInt128)

### Results

| File | Description |
|------|-------------|
| `summary_local.txt` | Summary statistics (local 8×8 pipeline) |
| `summary_super.txt` | Summary statistics (super 9×9 pipeline) |
| `ones.json` | Combined ones distribution for n=1..9 |
| `strong_upper_bound_8x8.json` | Strong upper bound calculations (8×8 base) |
| `strong_upper_bound_9x9.json` | Strong upper bound calculations (9×9 base) |
| `matrices_4x4.json` | Exact 4×4 SNS matrix database |
| `ber_n<N>.json` | Bernoulli(p) probability estimates for target N |

---

## Troubleshooting

### "Package not found" errors

Run package installation:
```bash
julia -e 'using Pkg; Pkg.instantiate()'
```

### Out of memory errors

Reduce the matrix size or increase available memory. For SLURM:
```bash
SBATCH_MEM=400G sbatch scripts/run_families.sh
```

### Thread-related issues

Explicitly set thread count:
```bash
JULIA_NUM_THREADS=4 julia src/local/families.jl
```

### Julia module not found on cluster

Specify the Julia module name:
```bash
JULIA_MODULE=julia/1.9 sbatch scripts/run_families.sh
```

---

## Quick Reference

```bash
# Local (8×8) - run directly with Julia
julia src/local/families.jl
julia src/local/aggregation.jl

# Supercomputer (9×9)
sbatch scripts/run_families.sh
sbatch scripts/run_aggregation.sh

# Utilities
julia src/utils/analyze_ones.jl
julia src/utils/compute_strong_upper_bound.jl [N]
julia src/utils/convert_to_json.jl local

# Monte Carlo
julia experiment/base_4x4.jl
julia experiment/ber(p).jl [N]
```
