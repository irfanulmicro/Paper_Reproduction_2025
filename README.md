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
conda install bioconda::shovill
sudo apt install rename

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
# Assembly
[Note: 
1.Assembly can be done by shovill and/or spades. Spades are integrated into shovill , in spades graph files are generated. 
2.Shovill Supports other assemblers like SKESA, Velvet, and Megahit.]
1.Shovill:
A. For single file: shovill [options] --outdir DIR (Directory name) --R1 R1.fq.gz --R2 R2.fq.gz 
B. For multiple file: 
(i) nano shovill.sh
{#!/bin/bash

# Loop through all *_R1.fastp.fastq.gz files in the current directory
for r1 in *_R1.fastp.fastq.gz
do
    # Derive corresponding R2 file and sample name
    sample=$(basename "$r1" _R1.fastp.fastq.gz)
    r2="${sample}_R2.fastp.fastq.gz"
    outdir="shovill_${sample}"

    echo "Running Shovill for sample: $sample"

    # Run Shovill with fixed genome size to avoid KMC memory error
    shovill --R1 "$r1" --R2 "$r2" --outdir "$outdir" --ram 3 --cpus 4 --gsize 5m [--ram 3 --cpus 4 --gsize 5m: needs when your pc is ram/processor sensitive]
done}
(ii) chmod +x shovill.sh
(iii) ./shovill.sh

2.Spades
A. For sigle file:
(i) Default: spades.py -o output_directory --phred-offset 33 --careful -1 strain1_R1.fastp.fastq.gz -2 strain1_R2.fastp.fastq.gz -k 21,33,55,77,99 
(ii) In case of limited resources: spades.py -o output_directory --phred-offset 33 -1 strain1_R1.fastp.fastq.gz -2 strain1_R2.fastp.fastq.gz -k 21,33,55 -t 4 -m 10 --only-assembler
B. For multiple files:
(i) nano run_spades_all.sh
{#!/bin/bash

# Parameters for SPAdes
KMER="21,33,55" [add "99" if you have high resource]
THREADS=4
MEMORY=10
PHRED=33

# Output base directory
OUTPUT_BASE="spades_outputs"
mkdir -p "$OUTPUT_BASE"

# Loop through each R1 file
for R1 in *_R1.fastp.fastq.gz; do
    # Extract sample name
    SAMPLE=$(basename "$R1" _R1.fastp.fastq.gz)
    
    # Match the corresponding R2 file
    R2="${SAMPLE}_R2.fastp.fastq.gz"

    # Check if R2 file exists
    if [[ ! -f "$R2" ]]; then
        echo "Missing R2 file for sample $SAMPLE. Skipping..."
        continue
    fi

    # Output directory for this sample
    OUTDIR="${OUTPUT_BASE}/${SAMPLE}"
    mkdir -p "$OUTDIR"

    echo "Running SPAdes for sample: $SAMPLE"
    
    # Run SPAdes
    spades.py -o "$OUTDIR" \
        --phred-offset $PHRED \
        -1 "$R1" -2 "$R2" \
        -k $KMER -t $THREADS -m $MEMORY \
        --only-assembler [no need when you have high resources]
done

echo "All SPAdes runs are done."}
(ii) chmod +x run_spades_all.sh
(iii) ./run_spades_all.sh








# Creating AutoMLSA Phylogeny
