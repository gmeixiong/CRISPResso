# CRISPResso
Reflow workflow for running CRISPResso


install docker...



docker run -itv ~/Documents/GitHub/CRISPResso:/scratch gmeixiong/crispresso bash



if you want s3 upload/download, run 'aws configure' and input the proper credentials...

python crispresso.py fastqs flash trim crisp [--s3]

python crispresso.py <datafolder containing fastq files, excel sample sheet, and adapterloc fa file> 
<flash output directory> <trimmomatic output directory> <crispresso output directory> [s3 flag for upload]
