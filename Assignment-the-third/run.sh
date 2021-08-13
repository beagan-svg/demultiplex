#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task 10

/usr/bin/time -v python3 twoplex.py
