#!/bin/bash

input=$1

bcftools norm -cx -f hg19.fa -Ov $input | bgzip > $input.fixed.gz

mv $input.fixed.gz $input

tabix -f -p vcf $input