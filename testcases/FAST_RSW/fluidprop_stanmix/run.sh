#!/bin/sh
#
#SBATCH --job-name=FP
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --account=research-ae-fpt

# Save current directory
CURRENT_DIR=$(pwd)

# Install FluidProp
cd /home/fneri/FluidProp
./Install_FluidProp
. /home/fneri/FluidProp/SetEnv_FluidProp.sh

# Return to original directory
cd "$CURRENT_DIR"

# Run Python case
python main.py
