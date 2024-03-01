# Comp-383-Pipeline-Project
We are comparing HCMV transcriptomes 2- and 6-days post-infection (dpi). I chose to do track 2 which focuses on genome assembly.

# Donor downloads: 
wget -O donor1_2dpi.fastq.gz "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030"
wget -O donor1_6dpi.fastq.gz "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033"
wget -O donor3_2dpi.fastq.gz "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044"
wget -O donor3_6dpi.fastq.gz "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"

# Build bowtie index:
bowtie2-build ncbi_dataset/data/genomic.DNA.fasta HCMV
bowtie2 --quiet -x HCMV -1 s-1.fastq -2 s-2.fastq -S Simap.sam --al-conc-gz s-mapped/_%fq.gz

# Track 2:
- Step 2: filter reads with Bowtie2
- Step 3: assemble transcriptomes with SPAdes
- Step 4: analysis the SPAdes assembly
- Step 5: Blast+ analysis

# Dependancies included:
- Bowtie2
- SPAdes
- BLAST+

# Quickly test pipeline using the small sample test data:
- modify the python scripts to use the sample test data
- run the python script
- verify the correct files output 
