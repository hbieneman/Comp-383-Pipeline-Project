import os

def download_fastq_file(url, output_fastq):
    command = f'wget {url} -O {output_fastq}'
    os.system(command)

# defining the URLs for each donor
donor1_2dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660030/SRR5660030"
donor1_6dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660033/SRR5660033"
donor3_2dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660044/SRR5660044"
donor3_6dpi_url = "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR5660045/SRR5660045"

# calling the function for each donor
download_fastq_file(donor1_2dpi_url, 'Donor1_2dpi.fastq')
download_fastq_file(donor1_6dpi_url, 'Donor1_6dpi.fastq')
download_fastq_file(donor3_2dpi_url, 'Donor3_2dpi.fastq')
download_fastq_file(donor3_6dpi_url, 'Donor3_6dpi.fastq')
