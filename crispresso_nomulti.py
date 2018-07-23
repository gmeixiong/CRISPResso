import pandas as pd
import argparse
import glob
import sys
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument('datafolder')
parser.add_argument('outdir')
args = parser.parse_args()

referenceseq = "96wp1sorted_sample_sheet.xlsx"
adapterloc = args.datafolder + "/TruSeq3-PE-2.fa"
readleng = "250"
fragleng = "300"
stdevfrag = "45"

foldertxt = args.datafolder + '/'
all_fq_names = glob.glob(foldertxt + '*R?_001.fastq.gz')



fq_df = pd.DataFrame({'all_fq_names' : all_fq_names})
fq_df['fq_fwd'] = None
fq_df['fq_rev'] = None

for i, row in fq_df.iterrows():
	if 'R1' in row.all_fq_names:
		row.fq_fwd = row.all_fq_names

for i, i_row in fq_df.iterrows():
	fq_name = i_row.all_fq_names
	if 'R2' in fq_name:
		match = fq_name.replace('R2', 'R1')
		for j, j_row in fq_df.iterrows():
			if j_row.all_fq_names == match:
				j_row.fq_rev = fq_name

fq_df = fq_df.iloc[:,1:3]
fq_df = fq_df.dropna()
# print(fq_df)

ref_import = pd.read_excel(foldertxt + referenceseq)

fq_df['reference'] = None
fq_df['HDR'] = None
fq_df['guideseq'] = None

for i, ref_row in ref_import.iterrows():
	sample_id = ref_row.Sample_ID
	# print(sample_id)
	for j, df_row in fq_df.iterrows():
		if sample_id in df_row.fq_fwd:
			df_row.reference = ref_row['WT amplicon']
			df_row.HDR = ref_row['HDR amplicon']
			df_row.guideseq = ref_row['Guide Sequence']

fq_names_only = pd.DataFrame()

fq_names = []
pattern = '^(.*)[_].*'
for i, row in fq_df.iterrows():
	filename = row.fq_fwd.replace(foldertxt, '')
	fq_name = re.findall(pattern, filename)[0]
	fq_names.append(fq_name)

# outdir = sys.argv[2]
if not os.path.exists(args.outdir):
	os.makedirs(args.outdir)

# print(fq_df)
# print(fq_names)

for i, row in fq_df.iterrows():
	# print i
	flash_cmd = '/root/CRISPResso_dependencies/bin/flash -r ' + str(readleng) + ' -f ' + str(fragleng) + ' -s ' + str(stdevfrag) + ' -o ' + fq_names[i-1] + ' -d ' + args.outdir + ' -z ' + row['fq_fwd'] + ' ' + row['fq_rev'] + " >> " + args.outdir + "/flash.txt"
	os.system(flash_cmd)


combined_fn = lambda row: os.path.basename(row.fq_fwd).replace('_001.fastq.gz', '.extendedFrags.fastq.gz')
qual_p_fn = lambda row: os.path.basename(row.fq_fwd).replace('_R1_001.fastq.gz', '_R1_qual_P.fastq.gz')
fq_df['combined'] = fq_df.apply(combined_fn, axis=1)
fq_df['qual_P'] = fq_df.apply(qual_p_fn, axis=1)

# filterdir = datafolder + '_filter'
# if not os.path.exists(filterdir):
# 	os.makedirs(filterdir)
combined_fn2 = lambda row: args.outdir + '/' + row.combined
fq_df['combined'] = fq_df.apply(combined_fn2, axis=1)
qual_p_fn2 = lambda row: args.outdir + '/' + row.qual_P
fq_df['qual_P'] = fq_df.apply(qual_p_fn2, axis=1)


###Removing adapters from reads, then passing them through quality filter, saving in new folder
##See trimmomatic for settings. A sliding window searches along for quality score at a minimum of 
##20 PHRED33 score, also trimming any adapter sequences.  
###############SET THESE VARIABLES###########################
leadingqual = 3
trailingqual = 3
slidewindsiz = 4
slidewindqual = 20
minlength = 50
############################################################
for i, row in fq_df.iterrows():
	trim_cmd = '/opt/conda/bin/trimmomatic SE ' + row.combined + ' ' + row.qual_P + ' ' + 'ILLUMINACLIP:' + adapterloc + ':2:30:10 LEADING:' + str(leadingqual) \
	+ ' TRAILING:' + str(trailingqual) + ' SLIDINGWINDOW:' + str(slidewindsiz) + ':' + str(slidewindqual) + ' MINLEN:' + str(minlength) + " >> " + args.outdir + "/trim.txt"
	os.system(trim_cmd)

# ##Clearing intermediate files
# remcmd <- paste0("find . -name '*extendedFrags.fastq.gz' -type f -delete")
# print(remcmd)
# message(remcmd, "\n"); system(remcmd)
# remcmd <- paste0("find . -name '*notCombined_[12].fastq.gz' -type f -delete")
# print(remcmd)
# message(remcmd, "\n"); system(remcmd)
# remcmd <- paste0("find . -name '*.hist' -type f -delete")
# print(remcmd)
# message(remcmd, "\n"); system(remcmd)
# remcmd <- paste0("find . -name '*.histogram' -type f -delete")
# print(remcmd)
# message(remcmd, "\n"); system(remcmd)

# os.system("ls -l " + datafolder + "_filter")

# ##Maping sorting and indexing bam files from sucessfully paired reads using bwa mem
# ##Only mapping paired reads, fwd, rev, and then placing them in the bm_fname folder location.

# print(fq_df)

for i, row in fq_df.iterrows():
	crispresso_cmd = '/opt/conda/bin/CRISPResso -r1 ' + row.qual_P + ' -a ' + row.reference + ' -e ' + row.HDR + ' -g ' + row.guideseq + ' -o ' + args.outdir + " >> " + args.outdir + "/crisp.txt"
	os.system(crispresso_cmd)
	# os.system("ls -l " + row.qual_P)
	# os.system("cd arg; ls -l ")




