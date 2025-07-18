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
conda install -c bioconda quast -y
sudo apt install rename
wget https://raw.githubusercontent.com/KorfLab/Assemblathon/refs/heads/master/assemblathon_stats.pl
wget https://raw.githubusercontent.com/ucdavis-bioinformatics/assemblathon2-analysis/refs/heads/master/FAlite.pm

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
1.Assembly execution:
a.Shovill:
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

b.Spades
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

2.Assebly stat check 
I.Aseembly statistics check execution:
a.For single file: 
(i) perl -I /home/irfan ~/assemblathon_stats.pl ~/assembly/selected/shovill_ERR10359916/contigs.fa (for shovill)
(iI) perl -I /home/irfan ~/assemblathon_stats.pl ~/assembly/selected/spades_outputs/ERR10359916/contigs.fasta (for spades)

b.For multiple file:
A. for shovill: 
(i) nano stat.sh
{for dir in ~/assembly/selected/shovill_ERR*/; do
    contig="$dir/contigs.fa"
    if [ -f "$contig" ]; then
        echo "Processing $contig"
        perl -I /home/irfan ~/assemblathon_stats.pl "$contig" > "${dir%/}/assemblathon_stats.txt"
    else
        echo "No contigs.fa found in $dir"
    fi
done}
(ii) chmod +x stat.sh
(iii) ./stat.sh

B. for spades:
(i) nano run_assemblathon_stats.sh
{#!/bin/bash

# Paths
PERL_SCRIPT=~/assemblathon_stats.pl
MODULE_PATH=/home/irfan
INPUT_BASE=~/assembly/selected/spades_outputs

# Loop over each ERR* directory
for dir in "$INPUT_BASE"/ERR*/; do
    contig_file="${dir}/contigs.fasta"
    output_file="${dir}/assemblathon_stats.txt"
    
    if [ -f "$contig_file" ]; then
        echo "Processing $contig_file"
        perl -I "$MODULE_PATH" "$PERL_SCRIPT" "$contig_file" > "$output_file"
    else
        echo "contigs.fasta not found in $dir"
    fi
done}
(ii) chmod +x run_assemblathon_stats.sh
(iii) ./run_assemblathon_stats.sh

II.If combine all output in a single tsv or excel file (optional):
A. for shovill:
(i) nano combine.sh
{output_file=~/assembly/selected/all_assemblathon_stats.tsv
> "$output_file"  # empty the file first

for dir in ~/assembly/selected/shovill_ERR*/; do
    contig="$dir/contigs.fa"
    sample=$(basename "$dir")
    if [ -f "$contig" ]; then
        echo ">>> $sample" >> "$output_file"
        perl -I /home/irfan ~/assemblathon_stats.pl "$contig" >> "$output_file"
        echo -e "\n" >> "$output_file"
    fi
done}
(ii) chmod +x combine.sh
(iii) ./combine.sh
B. for spades:
(i) nano combine_assemblathon_stats.sh
{#!/bin/bash

# Output file path (summary of all)
output_file=~/assembly/selected/spades_outputs/all_assemblathon_stats.tsv

# Empty the output file first
> "$output_file"

# Loop over each ERR* directory
for dir in ~/assembly/selected/spades_outputs/ERR*/; do
    contig="${dir}/contigs.fasta"
    sample=$(basename "$dir")

    if [ -f "$contig" ]; then
        echo "Processing $sample..."
        echo ">>> $sample" >> "$output_file"
        perl -I /home/irfan ~/assemblathon_stats.pl "$contig" >> "$output_file"
        echo -e "\n" >> "$output_file"
    else
        echo "⚠️  No contigs.fasta in $sample"
    fi
done

echo "✅ All stats written to: $output_file"}
(ii) chmod +x combine_assemblathon_stats.sh
(iii) ./combine_assemblathon_stats.sh

III. Visualization of statistics:
(i) cat assemblathon_stats.txt (for single view)
(ii) less assemblathon_stats.txt (for scrollable view)
(iii) nl assemblathon_stats.txt | less (with line numbers)

IV.Visualization of graph files originated from spades:
(i) (spades) irfan@User:~/assembly/selected/spades_outputs/ERR10359916$ Bandage image assembly_graph.fastg graph_output.png --height 1000 (graph to png conversion)
(ii) (spades) irfan@User:~/assembly/selected/spades_outputs/ERR10359916$ explorer.exe graph_output.png (Visualization for Windows)

# Quality assessment of assembly
1.Copying all directory (stored in spades output directory into new directory-quast)
(i) (base) irfan@User:~/assembly/selected/spades_outputs$ mkdir quast
(ii) (base) irfan@User:~/assembly/selected/spades_outputs$ ls
ERR10359916  ERR10359918  ERR10359921  ERR10359938  ERR10359946  ERR10359954  quast
(iii) (base) irfan@User:~/assembly/selected/spades_outputs$ cp -r ERR*/ quast/
[note: if i want to copy specific directory from spades output to quast then the command will be:
(base) irfan@User:~/assembly/selected/spades_outputs$ cp -r ERR10359916 ERR10359918 quast/]
(iv) (base) irfan@User:~/assembly/selected/spades_outputs$ cd quast
(v) (base) irfan@User:~/assembly/selected/spades_outputs/quast$ ls
ERR10359916  ERR10359918  ERR10359921  ERR10359938  ERR10359946  ERR10359954

2.Download a reference FASTA sequence and after download delete the metadata file (zone identifier)
(i) rm 'reference_fasta.fasta:Zone.Identifier'
3.Quast Exectuion:
[Note:
(a) Basic command: ./quast.py test_data/contigs_1.fasta \
           test_data/contigs_2.fasta \
        -r test_data/reference.fasta.gz \
        -g test_data/genes.txt \
        -1 test_data/reads1.fastq.gz -2 test_data/reads2.fastq.gz \
        -o quast_test_output
        
 (b) Additional parameter: 
 (I) -g <genes.txt>: If you have a gene annotation file in GFF or BED format. This helps QUAST assess gene completeness.
(II) -1 <reads_1.fastq> and -2 <reads_2.fastq>: If you have paired-end read files and want QUAST to align them to the assemblies for more in-depth evaluation (e.g., coverage analysis, misassemblies).
so, using -g and -1 can be used as follows:
quast.py ERR10359916/contigs.fasta \
         ERR10359918/contigs.fasta \
         -r reference_fasta.fasta \
         -g genes.gff \
         -1 reads_1.fastq.gz -2 reads_2.fastq.gz \
         -o quast_output
  (c) Confusion: 
  (i) quast.py ...	When quast.py is installed and in your PATH (like in your conda env) — most common case.
(ii) ./quast.py ...	When you are in the folder that contains the quast.py script file itself.]

(i) quast.py ERR*/contigs.fasta \
         -r reference_fasta.fasta \
         -o quast_output = Here ERR is the folder/directory and each directory contain contigs.fasta


(ii) (quast) irfan@User:~/assembly/selected/spades_outputs/quast/quast_output$ ls
aligned_stats    genome_stats    quast.log    report.tex  transposed_report.tex
basic_stats      icarus.html     report.html  report.tsv  transposed_report.tsv
contigs_reports  icarus_viewers  report.pdf   report.txt  transposed_report.txt
[see report.pdf file]













# Creating AutoMLSA Phylogeny
