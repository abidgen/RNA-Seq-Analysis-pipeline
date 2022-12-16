#!/bin/bash

#make directories
mkdir fastqc_files kallisto_output quality_control reference_genome 


# move fastqc.gz (big files) files to fastqc_files folder

mv *fastq.gz fastqc_files/

# move QC files (*fastqc.zip and *fastqc.html) into quality_control folder

mv *fastqc.zip *fastqc.html quality_control/

# move reference genome  files into reference_genome folder

mv Homo_sapiens.GRCh38.cdna.all.index Homo_sapiens.GRCh38.cdna.all.fa.gz Homo_sapiens.GRCh38.cdna.all.fa reference_genome\


# move kallisto outputs into kallisto_output folder

mv *.log kallisto_output\
mv HS01 kallisto_output\
mv HS02 kallisto_output\
mv HS03 kallisto_output\
mv HS04 kallisto_output\
mv HS05 kallisto_output\
mv CL08 kallisto_output\
mv CL10 kallisto_output\
mv CL11 kallisto_output\
mv CL12 kallisto_output\
mv CL13 kallisto_output\


## when done
echo "Done moving"