#!/bin/bash
#SBATCH --job-name=aggregation
#SBATCH --output=aggregation_%j.out
#SBATCH --error=aggregation_%j.err
#
# Default SLURM configuration for aggregation computation
# Override these by setting environment variables before submission:
#   SBATCH_PARTITION - partition name (default: compute)
#   SBATCH_MEM - memory per node (default: 200G)
#   SBATCH_CPUS - CPUs per task (default: 100)
#   SBATCH_TIME - time limit (default: 90-00:00:00)
#   SBATCH_NODES - node list (optional, e.g., "node[01-10]")
#   JULIA_MODULE - Julia module to load (optional, e.g., "julia/1.8.5")
#
# Usage:
#   sbatch scripts/run_aggregation.sh                    # Use defaults
#   SBATCH_PARTITION=gpu sbatch scripts/run_aggregation.sh  # Override partition
#

#SBATCH -p ${SBATCH_PARTITION:-compute}
#SBATCH --nodes=1
#SBATCH --mem=${SBATCH_MEM:-200G}
#SBATCH --cpus-per-task=${SBATCH_CPUS:-100}
#SBATCH --time=${SBATCH_TIME:-90-00:00:00}

# Load Julia module if specified
if [ -n "$JULIA_MODULE" ]; then
    module load "$JULIA_MODULE"
fi

# Set node list if specified
if [ -n "$SBATCH_NODES" ]; then
    # This needs to be set before sbatch parses the script
    # Use: SBATCH_NODES="node[01-10]" sbatch scripts/run_aggregation.sh
    : # Node list already set via SBATCH_NODES
fi

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK:-${SBATCH_CPUS:-100}}

julia src/super/aggregation.jl
