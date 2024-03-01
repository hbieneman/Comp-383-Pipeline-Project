import os
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

first_name = "Hailey"
last_name = "Bieneman"

output_directory = f"PipelineProject_{first_name}_{last_name}" # piece together name
os.system(f"mkdir {output_directory}") # making directory for output
os.chdir(output_directory)

# creating function to download fastq files using wget - slides
def download_fastq_file(url, output_fastq):
    command = f'wget {url} -O {output_fastq}'
    os.system(command)

# defining the URLs for each donor
# donor1_2dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030"
# donor1_6dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033"
# donor3_2dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044"
# donor3_6dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"
#
# # calling the function for each donor
# download_fastq_file(donor1_2dpi_url, 'Donor1_2dpi.fastq')
# download_fastq_file(donor1_6dpi_url, 'Donor1_6dpi.fastq')
# download_fastq_file(donor3_2dpi_url, 'Donor3_2dpi.fastq')
# download_fastq_file(donor3_6dpi_url, 'Donor3_6dpi.fastq')


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

#for when files are already downloaded
fastq_files = [
    'Donor1_2dpi.fastq',
    'Donor1_6dpi.fastq',
    'Donor3_2dpi.fastq',
    'Donor3_6dpi.fastq'
]

# making the log file
log_file_path = "PipelineProject.log"
# with open(log_file_path, "w") as log_file:

#for sample in ['Donor1_2dpi', 'Donor1_6dpi', 'Donor3_2dpi', 'Donor3_6dpi']: # for each sample in data
for fastq_file in fastq_files:
    # constructing the path
    #reads_fastq = f'{sample}.fastq'
    sample = fastq_file.split('.')[0]
    output_sam = f'{sample}_output.sam'
    unmapped_fastq = f'{sample}_unmapped.fastq'

    before_mapping_reads = count_reads_in_fastq(fastq_file) # calling function to count reads before mapping

    map_reads_to_hcmv_index(fastq_file, index_prefix, output_sam, unmapped_fastq) # mapping the reads

    after_mapping_reads = count_reads_in_fastq(unmapped_fastq) # calling function to count reads after mapping

    # Writing output to the log file
    with open(log_file_path, "a") as log_file:
        log_file.write(f"{sample} had {before_mapping_reads} read pairs before Bowtie2 filtering and {after_mapping_reads} read pairs after.\n")

#################### Part 3 ####################

# changing function to just output sam output files
def map_reads_to_hcmv_index(reads_fastq, hcmv_index_prefix, output_sam):
    command = f'bowtie2 -x {hcmv_index_prefix} -U {reads_fastq} -S {output_sam}'
    os.system(command)
    # -x: Index we’re mapping to.
    # -U: Used for unpaired reads. Specifies input fastq file.
    # -S: Name of output file. Writes a *.sam format file.

# extracts sample name and makes new output
for sam_file in ['Donor1_2dpi_output.sam', 'Donor1_6dpi_output.sam', 'Donor3_2dpi_output.sam', 'Donor3_6dpi_output.sam']:
    sample = sam_file.split('_output.sam')[0]
    output_sam = f'{sample}_output.sam' # output file path

    spades_assembly([sam_file], spades_output) #calls the spades_assembly function

    # Writing SPAdes command to the log file
    spades_command = f'spades.py --rna -s {sam_file} -o {spades_output}'
    # --rna: required for RNA-Seq data
    # - s: sam file format

    # appending the SPAdes command
    with open(log_file_path, "a") as log_file:
        log_file.write(f"SPAdes assembly command: {spades_command}\n")
    # -o: directory to store all the resulting files (required)

#################### Part 4 ####################

# calculate the number of contigs > 1000 bp and their total length
def calculate_contig_stats(fasta_file):
    contigs_gt_1000 = 0 # contigs > 1000 bp
    total_length_gt_1000 = 0 # total length of contigs > 1000 bp

    for record in SeqIO.parse(fasta_file, "fasta"): # for each
        contig_length = len(record.seq) # get the length of the contig

        # checking if the contig length is greater than 1000
        if contig_length > 1000:
            contigs_gt_1000 += 1
            total_length_gt_1000 += contig_length # add current length to total

    return contigs_gt_1000, total_length_gt_1000 # final counts

spades_assembly_fasta = 'SPAdes_Assembly/contigs.fasta' # path to SPAdes - fix this

contigs_gt_1000, total_length_gt_1000 = calculate_contig_stats(spades_assembly_fasta) # call function using SPAdes assembly

# write to the log file
with open(log_file_path, "a") as log_file:
    log_file.write(f"There are {contigs_gt_1000} contigs > 1000 bp in the assembly.\n")
    log_file.write(f"There are {total_length_gt_1000} bp in the assembly.\n")

#################### Part 5 ####################

def blast_analysis(assembly_dir):
    longest_contig = "" # define variables to store
    longest_contig_length = 0

    for filename in os.listdir(assembly_dir): # for each file in directory
        if filename.endswith(".fasta"): # for files ending in fasta
            with open(os.path.join(assembly_dir, filename), "r") as f:
                contig = "" # define variable
                contig_length = 0

                for line in f: # for each
                    if line.startswith(">"): # if line starts with >
                        # if longer than longest redefine
                        if contig_length > longest_contig_length:
                            longest_contig_length = contig_length
                            longest_contig = contig

                        # reset variables for loop
                        contig = ""
                        contig_length = 0
                    else:
                        contig += line.strip() # add line to contig and update length
                        contig_length += len(line.strip()) # remove white space

                # updates for longest contigs
                if contig_length > longest_contig_length:
                    longest_contig_length = contig_length
                    longest_contig = contig

                with open(log_file_path, "a") as log_file:
                    log_file.write(f"Longest Contig: {longest_contig}\n")
                    log_file.write(f"Longest Contig Length: {longest_contig_length}\n")

    return longest_contig, longest_contig_length # return longest contig and length

