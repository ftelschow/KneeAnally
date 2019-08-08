#!/bin/bash

for alpha in 0.15 0.1 0.05
do
for noise in 1 2 3
do
Rscript SO3CovRateSimulations.R $alpha $noise & 
done
done
