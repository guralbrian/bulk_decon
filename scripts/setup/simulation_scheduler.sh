#!/bin/bash
# Loop over models
for model in $(seq $1 $2); do
  # Loop over chunks
  for chunk in $(seq 1 20); do 
    sbatch scripts/setup/manual_simulations.sh $model $chunk
  done
done