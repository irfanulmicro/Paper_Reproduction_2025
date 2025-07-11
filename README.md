# Paper_Reproduction_2025
# Generating paper reproduction
Name of the paper: Akinyemi, K.O., Fakorede, C.O., Linde, J. et al. Whole genome sequencing of Salmonella enterica serovars isolated from humans, animals, and the environment in Lagos, Nigeria. BMC Microbiol 23, 164 (2023). https://doi.org/10.1186/s12866-023-02901-1

# Making directory for raw file
mkdir raw_seq

# Downloading raw sequence file:
1.Raw sequence was donwloaded from EMBL-EBL. download script was found from the respective site.
2.bash script: ena-file-download-read_run-PRJEB56537-fastq_ftp-20250710-0909.sh (found from EMBL-ENA website)
3.chmod +x ena-file-download-read_run-PRJEB56537-fastq_ftp-20250710-0909.sh
4. ./ena-file-download-read_run-PRJEB56537-fastq_ftp-20250710-0909.sh

# Deleting unnecessary files
rm wget-log* (wget log file removal)
find . -name "*:Zone.Identifier" -delete (Zone Identifier file removal)

# Installation of the tools
conda install bioconda::bbmap
conda install bioconda::fastqc
conda install bioconda::fastp
conda install bioconda::spades
conda install bioconda::bandage

# reformating / resampling raw sequence data (sample from my raw fastq samples, to elevate the next operation smoothly)
1.(genomics) irfan@User:~/raw_seq$ nano batch_sample.sh
[#!/bin/bash
for f in *_1.fastq.gz; do
  base=$(basename "$f" _1.fastq.gz)
  reformat.sh in1="${base}_1.fastq.gz" in2="${base}_2.fastq.gz" out1="${base}_sampled_R1.fq.gz" out2="${base}_sampled_R2.fq.gz" samplerate=0.1
done]

2.(genomics) irfan@User:~/raw_seq$ ./batch_sample.sh
[Error in ERR10359958_2.fastq.gz, line 2492185, with these 4 lines:
correction: 
a. repair.sh in1=ERR10359958_1.fastq.gz in2=ERR10359958_2.fastq.gz out1=ERR10359958_fixed_1.fq.gz out2=ERR10359958_fixed_2.fq.gz outs=unpaired_ERR10359958.fq.gz
b. reformat.sh in1=ERR10359958_fixed_1.fq.gz in2=ERR10359958_fixed_2.fq.gz out1=ERR10359958_sampled_R1.fq.gz out2=ERR10359958_sampled_R2.fq.gz samplerate=0.1]
3.(genomics) irfan@User:~/raw_seq$ mkdir sampled_files
4.(genomics) irfan@User:~/raw_seq$ mv *_sampled_*.fq.gz sampled_files/
5.(genomics) irfan@User:~/raw_seq$ ls sampled_files/

# Contamination check
1.(genomics) irfan@User:~/raw_seq/sampled_files$ nano run_sendsketch.sh
[#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p sketches

# Loop through all sampled files
for infile in *_sampled_*.fq.gz; do
    # Strip path and extension to create a dataset name
    dataset=$(basename "$infile" .fq.gz)

    echo "Processing $infile -> sketches/${dataset}.sketch"

    # Run sendsketch
    sendsketch.sh -da in="$infile" out="sketches/${dataset}.sketch"
done]
2.(genomics) irfan@User:~/raw_seq/sampled_files$ chmod +x run_sendsketch.sh
3.(genomics) irfan@User:~/raw_seq/sampled_files$ ./run_sendsketch.sh
4.Checking each sketch files: (genomics) irfan@User:~/raw_seq/sampled_files/sketches$ nano ERR10359916_sampled_R1.sketch (as for example)
5.Summary: 
a.(genomics) irfan@User:~/raw_seq/sampled_files/sketches$ nano summary.sh
[#!/bin/bash

# Create the header for summary file
echo 'strainID  ANI     gSize   taxName' > sketch_summary.tsv

# Loop through sketch report and extract necessary column from the first hit
for sketch in `ls -d *.sketch`;
do
        strainid=`basename "$sketch" | sed 's/.sketch//g'`
        ani=`tail -n +4 $sketch | head -n 1 | cut -f 3`
        gSize=`tail -n +4 $sketch | head -n 1 | cut -f 10`
        taxa=`tail -n +4 $sketch | head -n 1 | cut -f 12`
        echo "${strainid}       ${ani}  ${gSize}        ${taxa}" >> sketch_summary.tsv
done]
b.(genomics) irfan@User:~/raw_seq/sampled_files/sketches$ chmod +x summary.sh
c.(genomics) irfan@User:~/raw_seq/sampled_files/sketches$ ./ summary.sh
d.(genomics) irfan@User:~/raw_seq/sampled_files/sketches$ nano sketch_summary.tsv]

# Quality check 
1.(genomics) irfan@User:~$ cd quality
2.(genomics) irfan@User:~/quality$ cp ~/raw_seq/sampled_files/*_sampled_*.fq.gz ~/quality/
3.(genomics) irfan@User:~/quality$ ls
4.(genomics) irfan@User:~/quality$ fastqc *.fq.gz
5.(genomics) irfan@User:~/quality$ mkdir -p fastqc_zips fastqc_htmls
6.(genomics) irfan@User:~/quality$ mv *_fastqc.zip fastqc_zips/
7.(genomics) irfan@User:~/quality$ mv *_fastqc.html fastqc_htmls/
8.Visualization:
a. (genomics) irfan@User:~/quality/fastqc_htmls$ python3 -m http.server 8000
b. On browser address bar- http://localhost:8000
9.All fastqc reports in a single file: multiqc .

# Processing before assembly
1.(genomics) irfan@User:~/quality$ nano run_fastp.sh
[#!/bin/bash

for sample in *_sampled_R1.fq.gz; do
    base=${sample%_sampled_R1.fq.gz}
    echo "Processing sample $base"
    fastp \
      -i "${base}_sampled_R1.fq.gz" \
      -I "${base}_sampled_R2.fq.gz" \
      -o "${base}_sampled_R1.fastp.fastq.gz" \
      -O "${base}_sampled_R2.fastp.fastq.gz" \
      -e 28
done]
2.(genomics) irfan@User:~/quality$ chmod +x run_fastp.sh
3.(genomics) irfan@User:~/quality$ ./run_fastp.sh
4.Recheck the quality- fastqc *fastp.fastq.gz
5.Move the zip file and html file as follows:
a. (genomics) irfan@User:~/quality$ mkdir -p fastp_fastqc_htmls fastp_fastqc_zips
b.(genomics) irfan@User:~/quality$ mv *.fastp_fastqc.html fastp_fastqc_htmls/
c.(genomics) irfan@User:~/quality$ mv *.fastp_fastqc.zip fastp_fastqc_zips/




# Creating AutoMLSA Phylogeny
