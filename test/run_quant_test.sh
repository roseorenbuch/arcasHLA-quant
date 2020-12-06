#!/bin/bash
set -Eeuxo pipefail

## This script downloads one example BAM and runs the entire processing including quantification.
## It is meant to be used during development.
## Instructions: 
## - Checkout the repository and the required branch
## - Build the Docker image by running this from the repo directory: docker build -t arcashla-gen -f Docker/Dockerfile . 
## - Change to some empty directory outside of the repo and execute this script by its full path from the repo source

this_dir=$(cd $(dirname $0) && pwd)
HOST_SRC=$(cd $this_dir/.. && pwd)

## Download one sample BAM from 1000 genomes Geuvadis dataset used in the arcasHLA publication https://doi.org/10.1093/bioinformatics/btz474
curl -o HG00096.1.M_111124_6.bam https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/HG00096.1.M_111124_6.bam

## Expected HLA calls can be found, for example, here:
## https://github.com/genevol-usp/phasetools/blob/master/input/hla_genotypes.tsv

## where host pwd will be bind-mounted inside docker container
DOCKER_DIR=/home/work

## outputs of different stages inside container
DOCKER_EXTRACT_DIR=$DOCKER_DIR/extract
DOCKER_GENOTYPE_DIR=$DOCKER_DIR/genotype
DOCKER_QUANTREF_DIR=$DOCKER_DIR/quantref

## final quantification outputs inside docker - corresponds to `pwd`/quant
DOCKER_QUANT_DIR=$DOCKER_DIR/quant

## this is development script, so we bind mount `scripts` directory from host to avoid rebuilding the image on every code edit
DOCKER_BASE_CMD="docker run --rm -v `pwd`:$DOCKER_DIR -v $HOST_SRC/scripts:/home/arcasHLA-quant/scripts arcashla-gen:latest arcasHLA"

## extract chrom 6 from BAM
$DOCKER_BASE_CMD extract --paired $DOCKER_DIR/HG00096.1.M_111124_6.bam --outdir $DOCKER_EXTRACT_DIR

## genotype
$DOCKER_BASE_CMD genotype --outdir $DOCKER_GENOTYPE_DIR \
 $DOCKER_EXTRACT_DIR/HG00096.1.M_111124_6.extracted.1.fq.gz \
 $DOCKER_EXTRACT_DIR/HG00096.1.M_111124_6.extracted.2.fq.gz

## build custom reference for quantification
$DOCKER_BASE_CMD customize --genes A,B,C --transcriptome chr6 --genotype $DOCKER_GENOTYPE_DIR/HG00096.genotype.json --outdir $DOCKER_QUANTREF_DIR

## quantify
$DOCKER_BASE_CMD quant --ref $DOCKER_QUANTREF_DIR/HG00096 --outdir $DOCKER_QUANT_DIR \
 $DOCKER_EXTRACT_DIR/HG00096.1.M_111124_6.extracted.1.fq.gz \
 $DOCKER_EXTRACT_DIR/HG00096.1.M_111124_6.extracted.2.fq.gz

