# Paper_Reproduction_2025
# Generating paper reproduction
Name of the paper: Akinyemi, K.O., Fakorede, C.O., Linde, J. et al. Whole genome sequencing of Salmonella enterica serovars isolated from humans, animals, and the environment in Lagos, Nigeria. BMC Microbiol 23, 164 (2023). https://doi.org/10.1186/s12866-023-02901-1

# 1.Making directory for raw file
mkdir raw_seq

# 2.Downloading raw sequence file:
1.Raw sequence was donwloaded from EMBL-EBL. download script was found from the respective site.
2.bash script: ena-file-download-read_run-PRJEB56537-fastq_ftp-20250710-0909.sh (found from EMBL-ENA website)
3.chmod +x ena-file-download-read_run-PRJEB56537-fastq_ftp-20250710-0909.sh
4. ./ena-file-download-read_run-PRJEB56537-fastq_ftp-20250710-0909.sh

# 2.1 Deleting unnecessary files
rm wget-log* (wget log file removal)
find . -name "*:Zone.Identifier" -delete (Zone Identifier file removal)

# 3. Installation of the tools
[1.mamba is faster than conda
2.If you want to use mamba first install anaconda or miniconda and then mamba into it as follows-
conda install mamba -n base -c conda-forge]
conda or mamba install bioconda::bbmap
conda or mamba install bioconda::fastqc
conda or mamba install bioconda::fastp
conda or mamba install bioconda::spades
conda or mamba install bioconda::bandage
conda or mamba install bioconda::shovill
conda or mamba install -c bioconda quast -y
conda or mamba install -c conda-forge -c bioconda bakta
mamba create -n beav(environment name) beav(tool)
sudo apt install rename
wget https://raw.githubusercontent.com/KorfLab/Assemblathon/refs/heads/master/assemblathon_stats.pl
wget https://raw.githubusercontent.com/ucdavis-bioinformatics/assemblathon2-analysis/refs/heads/master/FAlite.pm
python3 -m pip install -U automlsa2
sudo apt install curl -y [curl for edirect installation, entrez-direct for query file preparation]
[entrez-direct installation code:
cd ~
curl -O https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/edirect.tar.gz
tar -xzf edirect.tar.gz
echo 'export PATH=$HOME/edirect:$PATH' >> ~/.bashrc
source ~/.bashrc]
pip install aniclustermap
conda install orthofinder -c bioconda
conda install pirate 
[dependencies: conda install r==3.5.1 r-ggplot2==3.1.0 r-dplyr==0.7.6 bioconductor-ggtree==1.14.4 r-phangorn==2.4.0 r-gridextra]

# optional dependencies for plotting figures in R
conda install r==3.5.1 r-ggplot2==3.1.0 r-dplyr==0.7.6 bioconductor-ggtree==1.14.4 r-phangorn==2.4.0 r-gridextra
mamba install -c conda-forge -c bioconda -c defaults 'panaroo>=1.3'
pip install scoary-2
pip install scoary

# 4.reformating / resampling raw sequence data (sample from my raw fastq samples, to elevate the next operation smoothly) [optional]
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

# 5.Contamination check
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

# 6.Quality check 
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

# 7.Processing before assembly and re-quality check 
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

# 8.Assembly
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

# 9. Quality assessment of assembly
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
(i) NCBI-Nucleotide-organism name-FASTA or NCBI-Nucleotide-organism name-WGS-link click-FASTA
(ii) rm 'reference_fasta.fasta:Zone.Identifier'
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

# 10.Annotation:
1. Database download: 
(A) By using bakta:
(i) script for database download to be saved in particular directory:
(a) (beav) irfan@User:~$ nano ~/bin/bakta-db
[#!/bin/bash
bakta_db download --output /home/irfan --type full/light]
(b) (beav) irfan@User:~$ chmod +x ~/bin/bakta-db
(c) (beav) irfan@User:~$ echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc
(d) (beav) irfan@User:~$  bakta-db
[Note: 
bakta-db is a script where path is defined by command ("echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc") . That's why it is executable from any folder/directory and the command is "beav-db" in stead of "./beav-db"]
                                      or
(ii) By tar file download and then extraction:
(a) wget https://zenodo.org/record/14916843/files/db-light.tar.xz
(b) bakta_db install -i db-light.tar.xz
(c) mkdir -p ~/bakta_db_light
(d) tar -xf db-light.tar.xz -C ~/bakta_db_light]
                                    or
(iii) Direct download (More applicable for full database):  bakta_db download --output /home/irfan --type full
[Note:
(a) sometimes it not works then follow above (i) or (ii)}
(b) after --output the path (/home/irfan) is the path where database is saved. 

(B).By using beav:
(i) beav_db

2.Execution of Annotation:
(i) For Bakta: 
(a) If i want to single organism's annotation separately: 
bakta --db /home/irfan/bakta_db_light/db-light \
      --output ~/assembly/selected/spades_outputs/ERR10359916/bakta_out \
      --prefix ERR10359916 \
      ~/assembly/selected/spades_outputs/ERR10359916/contigs.fasta
(b) If i want to multiple organism's annotation simultaneously:
(I) nano run_bakta_all.sh
[#!/bin/bash

BAKTA_DB="/home/irfan/bakta_db_light/db-light"
BASE_DIR=~/assembly/selected/spades_outputs

for sample_dir in $BASE_DIR/ERR*/; do
    sample=$(basename "$sample_dir")
    contig_file="$sample_dir/contigs.fasta"
    out_dir="$sample_dir/bakta_out"

    if [[ -f "$contig_file" ]]; then
        echo "Running Bakta on $sample..."
        bakta --db $BAKTA_DB --output $out_dir --prefix $sample $contig_file
    else
        echo "No contigs.fasta found in $sample_dir, skipping."
    fi
done]
(II) chmod +x run_bakta_all.sh
(III) ./run_bakta_all.sh

(ii) For Beav:
beav --input /home/irfan/assembly/selected/spades_outputs --threads 8 --tiger_blast_database home/irfan/blast/refseq_genomic.fna

3.Visualization:
(i) Bakta provided the follwoing file- bakta_plot contigs.json
(ii) Circos plot should be visualized in the follwoing server-  https://bakta.computational.bio/viewer

# 11. Alignment and phylogenetic tree construction (preferring Automlsa2 web server because it automatically select outgroup)
1.Rename of contigs.fasta of each directory and copy to particular directory for alignment: 
(i) (phylogeny) irfan@User:~/assembly/selected/spades_outputs$ ls
ERR10359916  ERR10359918  ERR10359921  ERR10359938  ERR10359946  ERR10359954 [Here every ERRXXXXX is a directory and each directory contain contigs.fasta file]
(ii) (phylogeny) irfan@User:~/assembly/selected/spades_outputs$ mkdir -p all_contigs
(iii) (phylogeny) irfan@User:~/assembly/selected/spades_outputs$ for dir in ERR*/; do
    cp "$dir/contigs.fasta" "all_contigs/${dir%/}_contigs.fasta"
done
2.Preparation of query file: The reference genome/protein sequence with which we will align our genome
[Note: in query file if two fasta have same sequence one should be excluded and replaced with new one]
(i) Manually:
(a) NCBI-Genome/Protein database- "gene name" and "organism name" in search box
(b) Go to FASTA-Copy amino acid sequence-paste into single file in Notepad (wsl supported ubuntu) 
                                          or 
                      I.nano ref_faa.sh or nano ref_faa.txt
                      II.paste all the sequencce with header (e.g- >bp xxx) and save
(ii) Command line:
(a) nano download_mlsa_salmonella.sh
[#!/bin/bash

# Define genes and organism
genes=("recA" "dnaA" "gyrB"  "16S rRNA")
organism="Salmonella enterica"
limit=20  # number of sequences per gene

# Output directory
mkdir -p mlsa_fasta
cd mlsa_fasta

# Download each gene
for gene in "${genes[@]}"; do
  echo "Downloading $gene for $organism..."
  outfile="${gene// /_}_Salmonella_complete.fasta"
  
  esearch -db protein -query "\"$gene\"[Gene] AND \"$organism\"[Organism] AND complete genome[Title]" | \
    efetch -format fasta | \
    awk '/^>/ {n++} n<='$limit' {print}' > "$outfile"

  echo "Saved: $outfile"
done

# Concatenate all into one file
cat *_Salmonella_complete.fasta > all_MLSA_Salmonella.fasta

echo "✅ All done. Combined FASTA: all_MLSA_Salmonella.fasta"]

3.Alignment and tree file generation :
(i) without Outgroup:(phylogeny) irfan@User:~/assembly/selected/spades_outputs$  automlsa2 --query ref_faa --dir /home/irfan/assembly/selected/spades_outputs/all_contigs -t 8 --allow_missing 4 --missing_check salmonella_paper 
Output:  salmonella_paper.nex.iqtree/  salmonella_paper.nex.treefile
(ii) With Outgroup:
(a) Outgroup bacteria selection: BVBRC-similar genome browser-In search box write your organism name-pick up a bacterial genus that is close to my organism and have genes of query file
(b) Environment preparation:
(I) keep outgroup sequence in that folder where assembled contigs.fasta are located. In this case 6 testing assembled seqeunce and 1 outgroup sequence were kept in all_contigs.folder
(II) keep ref_faa.txt(query files) another folder (in case of mine it was spades_output)
(c) Command execution:   (phylogeny) irfan@User:~/assembly/selected/spades_outputs$ automlsa2   --query ref_faa.txt   --dir all_contigs   --files all_contigs/ERRiiiiiiiiiiiii_contigs.fasta   --outgroup ERRiiiiiiiiiiiii_contigs   --allow_missing 7   --missing_check   -t 8   test_salmonella_with_ERRiiiiiiiiiiiii_outgroup2 
[I.If you work with single file then --files (flag) will be used
II. If you work with multiple file then --dir will be used (directory where assembled files are saved)

                                                                  Or 
(phylogeny) irfan@User:~/assembly/selected/spades_outputs$ automlsa2   --query ref_faa.txt   --dir all_contigs  --outgroup ERRiiiiiiiiiiiii_contigs   --allow_missing 7   --missing_check   -t 8   test_salmonella_with_ERRiiiiiiiiiiiii_outgroup3
Output: tree file

# 12.phylogenetic tree visualization (in R)
# 1.Necessary packages
#library(ggtree)
#library(phangorn)
#library(ape)
#library(BiocManager)
#library(Biocversion)
#library(ggplot2)
#library(ggtree)
#library(mgcv)
#library(phangorn)
#library(treeio)
#library(viridis)
#library(ggnewscale)

# 2.Creating a basic phylogenetic tree

sal.tree = read.tree('test_salmonella_with_ERRiiiiiiiiiiiii_outgroup3.nex.treefile')
sal.meta = read.table('metadata.csv', sep=',', header=T)

sal.tree = midpoint(sal.tree)

sal.tree = read.iqtree ('test_salmonella_with_ERRiiiiiiiiiiiii_outgroup3.nex.treefile')
sal.tree@phylo = midpoint(sal.tree@phylo)
plot(sal.tree@phylo)

# 3.Ornamenting phylogenetic tree with metadata
t1 = ggtree(sal.tree) %<+% sal.meta +
  geom_tippoint(aes(color = Source), size = 2)
t1

# 4.Circular phylogenetic tree
p <- ggtree(sal.tree, layout = "circular", size = 0.8) %<+% sal.meta +
  geom_tippoint(aes(color = Source), size = 2.5) +
  theme(legend.position = "right")
p

# 5.Add bootstrap support
p$data$bootstrap <- '0'
p$data[which(p$data$SH_aLRT >= 70 & p$data$UFboot  >= 70),]$bootstrap <- '1'

# 6.Add bootstrap overlay
p <- p + new_scale_color() +
  geom_tree(aes(color = bootstrap == "1"), size = 0.9) +
  scale_color_manual(name = "Bootstrap", values = c("TRUE" = "black", "FALSE" = "gray70"), guide = "none")
p

# 7.Showing metadata in heatmap
# 7.1 dataframe generation
(i) meta.sus <- data.frame(Susceptibility = sal.meta$Susceptiblity)
    rownames(meta.sus) <- sal.meta$Strain
                   or 
meta.sus <- as.data.frame(sal.meta[,'Susceptiblity'])
colnames(meta.Sus) <- 'Susceptiblity'
rownames(meta.Sus) <- sal.meta$Strain

(ii) meta.year <- data.frame(Year = sal.meta$Year)
    rownames(meta.year) <- sal.meta$Strain
              or 
meta.year <- as.data.frame(sal.meta[,'Year'])
colnames(meta.year) <- 'Year'
rownames(meta.year) <- sal.meta$Strain

# 7.2 Heat-map generation
(i) Three-stage strategy:
(a) For Categorical Variable:
 after defining [meta.sus <- data.frame(Susceptibility = sal.meta$Susceptiblity)
    rownames(meta.sus) <- sal.meta$Strain] write-
    p1 <- gheatmap(p, meta.sus, offset = 1.9, width = 0.32, colnames = FALSE) +
  scale_fill_viridis_d(option = "A", name = "Susceptibility") +
  new_scale_fill()

(b) For Continous Variable:
after defining [meta.year <- data.frame(Year = sal.meta$Year)
    rownames(meta.year) <- sal.meta$Strain] write-
    p2 <- gheatmap(p1, meta.year, offset = 2.6, width = 0.24, colnames = FALSE) +
  scale_fill_viridis(option = "D", name = "Year")

  (c) final_plot <- p2 +
  geom_tiplab2(aes(label = Country), offset = 0.60, size = 2.5, color = "gray40")
  (d) final_plot

  (ii) Two-stage strategy:
  (a) Categorical variable
p1 <- gheatmap(p, meta.sus, width = 9, offset = 0.51) + 
  scale_fill_viridis_d(option="A", name="Susceptiblity") +
  new_scale_fill()

(b) Continuous variable
p1 <- gheatmap(p1, meta.year, width= 7, offset=0.81) + 
  scale_fill_viridis(option="D", name="Year")
(c) p1 + geom_tiplab2(aes(label=Country), color="gray40", offset=0.06, size=3)
  
# 13.ANI clustermap
(phylogeny) irfan@User:~$ ANIclustermap -i  all_contigs -o ANI2

# 14.Orthologue clustering
1.Environment preparation:
(i)make a directory named "proteome" (any name)
(ii)From bakta annotation, all related faa files copy and paste into proteome directory
2.Execution: 
(i)(orthofinder) irfan@User:~/assembly/selected/spades_outputs$ orthofinder -f proteome/
output: 
(i) (orthofinder) irfan@User:~/assembly/selected/spades_outputs/proteome/OrthoFinder/Results_Jul26/Phylogenetic_Hierarchical_Orthogroups$ls
(ii) under the directory you will get- No.tsv, N1.tsv, N2.tsv, N3.tsv, N4.tsv and traits.txt
(iii) No.tsv file is more important

# 15.Pan-genome analysis
1.via Pirate (preferred for proper visualization, statitstics and presence-absence determination)
(A) Pan-genome analysis execution:
(i)Collection of gff3 files (should be qualityfull) and fna files from annotation:
agat_input_gff= gff files containing directory
agat_input_fna= fna files containing directory
(ii)stadardization of gff3 files into gff files via (agat) 
(note:bash script: save it before the gff3 containing folder)
(a)nano pangenome.sh
[#!/bin/bash
# Paths to input folders
gff_dir="agat_input_gff"
fna_dir="agat_input_fna"

# Output folder (optional)
mkdir -p agat_output

# Loop over all .gff files in the gff input directory
for infile in ${gff_dir}/*.gff; do
    # Extract strain name by removing path and extension
    strain=$(basename "$infile" .gff)
    
    echo "Processing: $strain"

    # Convert to GFF3 using AGAT
    agat_convert_sp_gxf2gxf.pl -g "$infile" -o "agat_output/${strain}.gff" --gff_version_output 3

    # Append the FASTA content
    echo "##FASTA" >> "agat_output/${strain}.gff"
    cat "${fna_dir}/${strain}.fna" >> "agat_output/${strain}.gff"
done]
[note: if your multiple sequence have multiple format like .fna, .fa, .fasta then you should modify the last line as follows:-
cat "${fna_dir}/${strain}.f*" >> "agat_output/${strain}.gff"]

(b)chmod +x pangenome.sh
(c)./pangenome.sh
(iii) mkdir log
(iv) mv *agat.log log
 (v) (orthologue) irfan@User:~/exp$ PIRATE -i  agat_output/ -o output/ -t 16
 here, agat_output = is the directory where standardized gff files are located
 output: Speficially prefer on- 
 (a) PIRATE.log    
 (b) PIRATE.gene_families.ordered.tsv 
(B) Pan-genome visualization :via R
(i)Environment preparation: copy "PIRATE.gene_families.ordered.tsv" and save it to the respective "R-studio" working directory
(ii)Open a new "R-script" and follow the belowing code:
library(data.table)
library(dplyr)
library(ggplot2)
library(doParallel)
library(ggpubr)
PIRATE.gene_families = as.data.frame(fread("PIRATE.gene_families.ordered.tsv"))
PIRATE.gene_families$number_genomes <- 0 
PIRATE.gene_families$number_genomes <- rowSums( !(PIRATE.gene_families[ , 23:length(colnames(PIRATE.gene_families))] == "") )
PIRATE.gene_families <- PIRATE.gene_families[ which(PIRATE.gene_families$number_genomes > 0), ]
head(PIRATE.gene_families)
gene_total <- colSums( !(PIRATE.gene_families[ , 23:length(colnames(PIRATE.gene_families))] == ""))
summary(as.vector(gene_total)) 
sd(gene_total)
hist(gene_total)

n_gene_fams<-  nrow(PIRATE.gene_families)
n_gene_fams

cols_to_exclude <- c("allele_name", "gene_family", "consensus_gene_name", "consensus_product", "threshold", "alleles_at_maximum_threshold", "number_genomes", "average_dose", "min_dose", "max_dose", "genomes_containing_fissions", "genomes_containing_duplications", "number_fission_loci",                    "number_duplicated_loci", "no_loci", "products",                    "gene_names", "min_length(bp)", "max_length(bp)",                    "average_length(bp)", "cluster", "cluster_order")

strains_only <- PIRATE.gene_families[, !names(PIRATE.gene_families) %in% cols_to_exclude]
strain_names <- names(strains_only)

ngenomes <- length(unique(strain_names))
ngenomes

n_gene_fams_core_all <- sum(PIRATE.gene_families$number_genomes == ngenomes)
n_gene_fams_core_all

n_gene_fams_core_all <- sum(PIRATE.gene_families$number_genomes == ngenomes)
n_gene_fams_core_all

core_all_percentage <- (n_gene_fams_core_all * 100) / n_gene_fams
core_all_percentage

cutoff <- round(.95* ngenomes)

n_gene_fams_core_w95per <- sum(PIRATE.gene_families$number_genomes >= cutoff)
n_gene_fams_core_w95per

core_gene_percentage <- (n_gene_fams_core_w95per * 100) / n_gene_fams
core_gene_percentage

n_gene_fams_singletons <- sum(PIRATE.gene_families$number_genomes == 1)
n_gene_fams_singletons

accesory_gene_percentage <- (n_gene_fams_singletons*100)/n_gene_fams
accesory_gene_percentage

n_accessory <- n_gene_fams - (n_gene_fams_singletons + n_gene_fams_core_w95per)
n_accessory

(n_accessory*100)/n_gene_fams

n_accessory / ngenomes


# Visualize the distribution of gene presence in a gene family (distribution of core to accessory genes)

# Prepare dataframe for plotting
gene_fam_totals <-as.data.frame(PIRATE.gene_families$number_genomes)

colnames(gene_fam_totals) <- 'count'

# Categorize gene groups
gene_fam_totals$group = 0

for (i in 1:nrow(gene_fam_totals)){  
  if (gene_fam_totals$count[i] == 1) {
    gene_fam_totals$group[i] = "Singleton"  
  } else if (gene_fam_totals$count[i] >= .95*ngenomes){
    gene_fam_totals$group[i] = "Core"  
  } else {
    gene_fam_totals$group[i] = "Accessory"  
  }
}

p1 <- gene_fam_totals %>%  ggplot(aes(x = count)) +  
  geom_bar(aes(fill = group), position = "identity") +  
  ggtitle("Pangenome Distribution by gene ortholog group") +  
  ylab("n gene families") + 
  xlab("n genomes in family") +  
  theme_classic() +  
  theme(plot.title = element_text(size=15), legend.title = element_blank(), legend.position.inside = c(0.85, 0.85)) +  
  scale_fill_manual(name = "", values = c("blueviolet", "coral",  "darkcyan")) +          
  theme_minimal()
p1

#make df of totals
fam_dist_df <- data.frame(
  category=c("Singleton", "Accessory", "Core"),
  count=c(n_gene_fams_singletons, n_gene_fams -(n_gene_fams_singletons + n_gene_fams_core_w95per), n_gene_fams_core_w95per)
)

# Compute percentages
fam_dist_df$fraction <- fam_dist_df$count / sum(fam_dist_df$count)
fam_dist_df$fraction

#Compute the cumulative percentages (top of each rectangle)
fam_dist_df$ymax <- cumsum(fam_dist_df$fraction)
fam_dist_df$ymax

# Compute the bottom of each rectangle
fam_dist_df$ymin <- c(0, head(fam_dist_df$ymax, n=-1))
fam_dist_df$ymin


# Compute label position
fam_dist_df$labelPosition <- (fam_dist_df$ymax + fam_dist_df$ymin) / 2
fam_dist_df$labelPosition

# Compute a good label
fam_dist_df$label <- paste0(fam_dist_df$category, "\n", fam_dist_df$count)
fam_dist_df$label



# Plot
p2 <- ggplot(fam_dist_df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +  
  geom_rect() +  
  geom_text( x=1.5, aes(y=labelPosition, label=label, color=category), size=4) + 
  # x here controls label position (inner / outer)  
  scale_fill_manual(values=c("blueviolet", "coral",  "darkcyan")) +  
  scale_color_manual(values=c("blueviolet", "coral",  "darkcyan")) +  
  coord_polar(theta="y") +  xlim(c(-1, 4)) +  theme_void() +
  theme(legend.position = "none")

p2

# Gene Accumulation Curve

binary_df <- strains_only %>%   mutate_all(~ as.numeric(nzchar(.)))
binary_df

count_all_ones <- function(df) {  
  num_cols <- ncol(df)  
  count_cases <- df %>%    
    filter(rowSums(.) == num_cols) %>%                                    
    nrow()  
  return(count_cases)}

count_all_ones

find_common_genes <- function(df, cols) {  
  common_genes <- count_all_ones(df[cols])
  return(common_genes)}

# Initialize dataframe to store results
results <- data.frame(  NumGenomes = integer(),  CommonGenes = integer(),  TotalGenes = integer())

# Initialize an empty list to store results
result_list <- list()

library(parallel)

# Setup parallel backend to use many processorscores=detectCores()

cores <- detectCores()
cl <- makeCluster(cores - 2) #not to overload your computer
registerDoParallel(cl) 

# Generate combinations and count common genes
# Initialize an empty list to store results

for(num_genomes in 2:ncol(binary_df)) {  
  rand_combs <- replicate(1000, sample(names(binary_df), num_genomes), simplify = FALSE)
  print(num_genomes)  
  
  local_results <- foreach(i = 1:length(rand_combs), .combine = rbind, .packages = 'dplyr') %dopar% {    
    common_genes_count <- find_common_genes(binary_df, rand_combs[[i]])    
    tota_genes_count <- sum(rowSums(binary_df[rand_combs[[i]]]) > 0)    
    data.frame(NumGenomes = num_genomes, CommonGenes = common_genes_count, TotalGenes = tota_genes_count)  }  
  result_list[[num_genomes - 1]] <- local_results  # Store results in list
}

# Stop the parallel backend
stopCluster(cl)

# Combine all results into one data frame
results <- do.call(rbind, result_list)

# Fit power law

# Log-transform the data

results <- results %>%   
  mutate(log_NumGenomes = log(NumGenomes), log_CommonGenes = log(CommonGenes), log_TotalGenes = log(TotalGenes))

head(results)

#Fit a linear model to the log-transformed data
lm_model_common <- lm(log_CommonGenes ~ log_NumGenomes, data = results)
lm_model_all <- lm(log_TotalGenes ~ log_NumGenomes, data = results)

# Get the coefficients
coefficients_common <- coef(lm_model_common)
intercept_common <- coefficients_common[1]
slope_common <- coefficients_common[2]

coefficients_all <- coef(lm_model_all)
intercept_all <- coefficients_all[1]
slope_all <- coefficients_all[2]

# Create a data frame with the fitted power-law curve
fitted_curve <- results %>%
  mutate(FittedCommonGenes = exp(intercept_common) * NumGenomes^slope_common)

# Estimate core genome size based on the accumulation curve
core_genome_size <- exp(intercept_common) * max(results$NumGenomes)^slope_common

fitted_curve <- fitted_curve %>%  
  mutate(FittedAllGenes = exp(intercept_all) * NumGenomes^slope_all)

p3 <- ggplot(results, aes(x = NumGenomes, y = CommonGenes)) +  
  geom_point(aes(color = "Common Genes"), size = 3) +  
  geom_line(data = fitted_curve, aes(x = NumGenomes, y = FittedCommonGenes, color = "Fitted Curve"), size = 1) +  
  scale_color_manual(name = "",                     values = c("Common Genes" = "grey", "Fitted Curve" = "red")) +  
  labs(title = "Core Gene Accumulation Curve",       x = "Number of Genomes Compared",       y = "Number of Common Genes") +  
  theme_minimal() +  
  theme(legend.position = c(0.95, 0.95),  legend.justification = c("right", "top"))

p3



p4 <- ggplot(results, aes(x = NumGenomes, y = TotalGenes)) +  
  geom_point(aes(color = "Total Genes"), size = 3) +  
  geom_line(data = fitted_curve, aes(x = NumGenomes, y = FittedAllGenes, color = "Fitted Curve"), size = 1) +  
  scale_color_manual(name = "",                     values = c("Total Genes" = "grey", "Fitted Curve" = "red")) +  
  labs(title = "Pan Genome Accumulation Curve",       x = "Number of Genomes Compared",       y = "Number of Total Genes") +  
  theme_minimal() +  
  theme(legend.position = c(0.95, 0.05), legend.justification = c("right", "bottom"))

p4

2.via Panaroo (preferred for only pangenome statistics)
(A) collect gff3 files of annotation from prokka (better) or bakta
(B) run Panaroo: panaroo -i *.gff -o results --clean-mode strict --remove-invalid-genes
[output: final_graph.gml]
(C) Visualization: https://cytoscape.org/

3.Web server: https://nmdc.cn/ipga/

# 16.gene presence-absence determination
1.Option1/no trait file:
I.Preparation: collect the No.tsv file from orthofinder2 output
II.Execution:  scoary2 \
    --genes N0.tsv \
    --gene-data-type 'gene-list:\t' \
    --n-permut 300 \
    --n-cpus 7 \
    --outdir out
*Here, -n-permut and --n-cpus are customizable on basis of your pc's ram+memory consideration

2.Option2/trait file (gaussian/continuous data):
I. Preparation:
(i)  collect the No.tsv file from orthofinder2 output
(ii) manually create a trait_numeric.tsv file (strain name should match with column name of N0.tsv file)
*No.tsv file (column name) and trait_numeric.tsv file (strain name) should match otherwise can match with python script
II.Execution:
(i) (scoary-2) irfan@User:~/assembly/selected/spades_outputs/proteome/OrthoFinder/Results_Jul26/Phylogenetic_Hierarchical_Orthogroups$ scoary2 \
    --genes N0.tsv \
    --gene-data-type 'gene-list:\t' \
    --traits traits_numeric.tsv \
    --trait-data-type 'gaussian:\t' \
    --n-permut 100 \
    --n-cpus 4 \
    --outdir out
III. output:
(i) (scoary-2) irfan@User:~/assembly/selected/spades_outputs/proteome/OrthoFinder/Results_Jul26/Phylogenetic_Hierarchical_Orthogroups/out$ ls
     app  binarized_traits.tsv  logs  traits  tree.nwk
*Here-
(A) binarized_traits.tsv (file)-
(a) Shows how Scoary2 converted your numeric traits into binary 0/1 values
(b) Based on its Gaussian model and threshold (cutoff=0.85)
(c) You can inspect this to see which samples were considered “high” or “low” for each trait
(B) tree.nwk (file)-
(a) Newick-formatted phylogenetic tree built from your gene presence/absence matrix
(b) You can open this in a tree viewer like:
-iTOL
-FigTree
-Dendroscope
(C) traits/ (directory/folder)
(a) Contains per-trait analysis files, including intermediate results
_(b) If any traits had passed filtering, Scoary would generate:
-*.sig.tsv = significant genes
=*.pairs.tsv = gene-gene pairs for co-association
(c) 
-Under triats directory there will be subdirectories (subdirectory number = number of traits used in traits_numeric.tsv file. In my case i used two traits: Resistant Score and Growth rate) e.g-
irfan@User:~/assembly/selected/spades_outputs/proteome/OrthoFinder/Results_Jul26/Phylogenetic_Hierarchical_Orthogroups/out/traits/ResistanceScore$ ls
result.tsv
-in each subdirectory there contains results.tsv file and the file contains as follows-
contingency_table	The 2x2 contingency table in format (g+t+, g+t-, g-t+, g-t-), counts of gene presence/trait presence combinations
fisher_p	Fisher’s exact test p-value assessing association between gene presence and trait
odds_ratio	Odds ratio measuring strength of association; inf means a perfect association
Gene	The gene or orthogroup ID being tested
g+t+	Number of strains with gene present AND trait present
g+t-	Number of strains with gene present AND trait absent
g-t+	Number of strains with gene absent AND trait present
g-t-	Number of strains with gene absent AND trait absent
* Unlike binary traits where Scoary gives multiple test corrections like Bonferroni and FDR by default, for continuous traits Scoary does not always output corrected p-values automatically.

# 17.snp calling
1.Select a good reference (in FASTA format) which have good quality (hybrid assembly, N50, L50). I choose the following organism-
>CP196174.1 Salmonella enterica strain W126 chromosome, complete genome
2.keep your both fastq (R1 and R2 files) in same folder. such as: 
(snippy) irfan@User:~/snp/Salmonella_snp$ ls
ERR10359916_1.fastq.gz  ERR10359918_2.fastq.gz  ERR10359938_1.fastq.gz  ERR10359946_2.fastq.gz
ERR10359916_2.fastq.gz  ERR10359921_1.fastq.gz  ERR10359938_2.fastq.gz  ERR10359954_1.fastq.gz
ERR10359918_1.fastq.gz  ERR10359921_2.fastq.gz  ERR10359946_1.fastq.gz  ERR10359954_2.fastq.gz
3.(snippy) irfan@User:~/snp$ nano run_snippy_all.sh
[#!/bin/bash

# Paths
REF="Salmonella_enterica_ref.fasta"
READDIR="Salmonella_snp"
OUTDIR="mysnps_results"


# Make output dir
mkdir -p $OUTDIR

# List of sample prefixes (assuming _1.fastq.gz and _2.fastq.gz)
for r1 in ${READDIR}/*_1.fastq.gz; do
    sample=$(basename $r1 _1.fastq.gz)
    r2="${READDIR}/${sample}_2.fastq.gz"

    # Skip if pair not found
    [ -f "$r2" ] || { echo "Missing $r2"; continue; }

    # Run snippy
    snippy --cpus 4 \
        --outdir ${OUTDIR}/${sample} \
        --ref $REF \
        --R1 $r1 \
        --R2 $r2 \
        --force
done
note: here cpus 4/8 depends on the processor/ram of pc]
4.(snippy) irfan@User:~/snp$ chmod +x run_snippy_all.sh
5.(snippy) irfan@User:~/snp$ ./run_snippy_all.sh
6.(snippy) irfan@User:~/snp$ snippy-core --ref Salmonella_enterica_ref.fasta --prefix core mysnps_results/*
7.(phylogeny) irfan@User:~/snp$ fasttree -nt core.full.aln > core.full.tree (optional, if you to construct phylogenetic tree by snp)

output:
1.vcf file
2.core.aln file
3.tree file (visualization by itol web server)

# PCA analysis using snp
# install.packages("adegenet")
library("adegenet")
library("ggplot2")

snp <- fasta2genlight('core.full.aln', snpOnly=T)
meta <- read.table('meta.csv', sep=',', header = T)

snp$ploidy

snp$pop = as.factor(meta[match(snp$ind.names, meta$Strain),]$SamplingSite)

# Heat map of genotype
glPlot(snp)

# PCA_basic
pca <- glPca(snp, nf=10)
barplot(100*pca$eig/sum(pca$eig), main="Eigenvalues", col=heat.colors(length(pca$eig)))

# PCA_classic
pca_full <- glPca(snp, nf = nInd(snp))  # or use nf=500 to be safe
barplot(100 * pca_full$eig / sum(pca_full$eig), 
        main = "Eigenvalues", 
        col = heat.colors(length(pca_full$eig)),
        border = NA)

# fast plot
scatter(pca, psi='bottomright')

# Optimization of the plot

pca.dataset = as.data.frame(pca$scores)
pca.dataset$isolates = rownames(pca.dataset)
pca.dataset$pop = as.factor(meta[match(snp$ind.names, meta$Strain),]$SamplingSite)
pca.dataset$spp = as.factor(meta[match(snp$ind.names, meta$Strain),]$Species)

ggplot(pca.dataset, aes(PC1, PC2, fill=spp)) + geom_point(shape=21, size=3, alpha=0.7)
ggplot(pca.dataset, aes(PC1, PC3, fill=spp)) + geom_point(shape=21, size=3, alpha=0.7)
# Fst analysis using snp
# install.packages("hierfstat")
library("hierfstat")
library("adegenet")

# Lets load the fasta file by extract only the SNPs
seq.snp = fasta2DNAbin('core.full.aln', snpOnly=T)
obj = DNAbin2genind(seq.snp)

# Load the meta file that contains population structure 
meta <- read.table('meta.csv', sep=',', header = T)

# Add the structure I want to test in the GenInd object
obj$pop = as.factor(meta[match(rownames(obj$tab), meta$Strain),]$SamplingSite)

# Convert genind object into hierfstat object
obj.hf = genind2hierfstat(obj)

# Various estimation of Fst
genet.dist(obj.hf, method='Nei87', diploid=F)
genet.dist(obj.hf, method='Dch', diploid=F)





