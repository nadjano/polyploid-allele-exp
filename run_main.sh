#!/bin/bash

# Run the main script
module load nextflow
nextflow run main.nf  -bg -resume -with-report report.html 