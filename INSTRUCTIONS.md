# Instructions for Running LU Factorization of Discrete Random Matrices

This document provides detailed instructions for running all components of this project from the repository root.

## Table of Contents

1. [Prerequisites](#prerequisites)
2. [Project Overview](#project-overview)
3. [Running Locally](#running-locally)
4. [Running on a Supercomputer (SLURM)](#running-on-a-supercomputer-slurm)
5. [Utility Scripts](#utility-scripts)
6. [Output Files](#output-files)

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
- Generates summary to `results/local/summary.txt`

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
- Summary: `results/super/summary.txt`

---

## Utility Scripts

All utility scripts are located in `src/utils/` and accept a mode argument (`local` or `super`).

### Analyze Ones Distribution

Computes the distribution of ones (non-zero entries) in the matrices:

```bash
julia src/utils/analyze_ones.jl local
julia src/utils/analyze_ones.jl super
```

**Output:** `results/<mode>/plots/ones_NxN.json`

### Compute Upper Bound

Computes upper bounds on probabilities for larger matrix dimensions:

```bash
julia src/utils/compute_upper_bound.jl local           # Based on 8×8 data, N=30 (default)
julia src/utils/compute_upper_bound.jl super           # Based on 9×9 data, N=30 (default)
julia src/utils/compute_upper_bound.jl local 50        # Based on 8×8 data, N=50
julia src/utils/compute_upper_bound.jl super 100       # Based on 9×9 data, N=100
```

**Parameters:**
- First argument: `local` or `super` (default: `local`)
- Second argument: `N` - target matrix dimension (default: `30`)

**Output:** `results/<mode>/plots/upper_bound_NxN.json`

### Convert Binary to JSON

Converts binary data files to human-readable JSON format:

```bash
julia src/utils/convert_to_json.jl local
julia src/utils/convert_to_json.jl super
```

**Output:** JSON files in `data/<mode>/readable/`

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

**Counts file (9×9):** Each record is 27 bytes:
- 9 bytes: Row sums (9 × UInt8)
- 9 bytes: Column sums (9 × UInt8)
- 9 bytes: Count (UInt128, padded)

### Results

| File | Description |
|------|-------------|
| `summary.txt` | Summary statistics |
| `plots/ones_NxN.json` | Ones distribution data |
| `plots/upper_bound_NxN.json` | Upper bound calculations |

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
julia src/utils/analyze_ones.jl local
julia src/utils/compute_upper_bound.jl local [N]
julia src/utils/convert_to_json.jl local
```
