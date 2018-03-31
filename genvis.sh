#!/bin/bash
#SBATCH --time=0:0:40
#SBATCH --job-name=fractalStitch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=10
#SBATCH --partition=YOUR-PARTITION

ffmpeg -loglevel fatal -framerate 10 -i %04d.png -vf "scale=1000:1000,format=yuv420p" -c:v libx264  out.mp4