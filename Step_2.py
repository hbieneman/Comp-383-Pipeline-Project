import os
genome_fasta = "/home/hbieneman/PipelineProject/NC_006273.2.fasta" # defining the path to fasta file
index_prefix = 'HCMV_index' # defining the index name

def count_reads_in_fastq(fastq_file):
    count = 0 # set count to zero
    with open(fastq_file, 'rb') as file: # open files in binary - looked this up
        for line in file:
            if line.startswith(b'@'): # allows us to read binary data
                count += 1 # counting each read in the fastq file
    return count

def create_bowtie2_index(index_path, genome_fasta): # function to create index for HCMV
    os.system(f"bowtie2-build {genome_fasta} {index_path}") # from slides

# creating index using the specific fasta file
create_bowtie2_index(index_prefix, genome_fasta)

# function to save  reads that map to the HCMV index
def map_reads_to_hcmv_index(reads_fastq, hcmv_index_prefix, output_sam, unmapped_fastq):
    command = f'bowtie2 -x {hcmv_index_prefix} -U {reads_fastq} -S {output_sam} --un {unmapped_fastq}' # slides
    os.system(command)
    # -x: Index we’re mapping to.
    # -U: Used for unpaired reads. Specifies input fastq file.
    # -S: Name of output file. Writes a *.sam format file.
    # --un: Specifies where unmapped reads will be written.

# for when files are already downloaded
# fastq_files = [
#     'Donor1_2dpi.fastq',
#     'Donor1_6dpi.fastq',
#     'Donor3_2dpi.fastq',
#     'Donor3_6dpi.fastq'
# ]

# making the log file
log_file_path = "PipelineProject.log"
# with open(log_file_path, "w") as log_file:

for sample in ['Donor1_2dpi', 'Donor1_6dpi', 'Donor3_2dpi', 'Donor3_6dpi']: # for each sample in data
#for fastq_file in fastq_files:
    # constructing the path
    reads_fastq = f'{sample}.fastq'
    #sample = fastq_file.split('.')[0]
    output_sam = f'{sample}_output.sam'
    unmapped_fastq = f'{sample}_unmapped.fastq'

    before_mapping_reads = count_reads_in_fastq(reads_fastq) # calling function to count reads before mapping

    map_reads_to_hcmv_index(reads_fastq, index_prefix, output_sam, unmapped_fastq) # mapping the reads

    after_mapping_reads = count_reads_in_fastq(unmapped_fastq) # calling function to count reads after mapping

    # Writing output to the log file
    with open(log_file_path, "a") as log_file:
        log_file.write(f"{sample} had {before_mapping_reads} read pairs before Bowtie2 filtering and {after_mapping_reads} read pairs after.\n")
