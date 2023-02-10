#!/usr/bin/zsh
sbatch -t 24:00:00 --cpus-per-task=1 --mem=64G run_batch_simple.sh
