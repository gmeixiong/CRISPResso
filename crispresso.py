from multiprocessing import Process
import pandas as pd
import argparse
import glob
import re
import os

def parseArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument('datafolder')
	parser.add_argument('flashdir')
	parser.add_argument('trimdir')
	parser.add_argument('crispdir')
	parser.add_argument('--s3', action="store_true")
	return parser.parse_args()

def initDataFrame(fastq_folder):
	fastq_names = glob.glob(fastq_folder + '*R?_001.fastq.gz')
	df = pd.DataFrame({'fq_name' : fastq_names})
	df['fq_fwd'] = None
	df['fq_rev'] = None
	pairReadFiles(df)
	df = df.iloc[:,1:3]
	df = df.dropna()
	return df

def pairReadFiles(df):
	for _, r1_row in df[df.fq_name.str.contains('R1')].iterrows():
		df.loc[df.fq_name == r1_row.fq_name, 'fq_fwd'] = r1_row.fq_name
			
	for _, r2_row in df[df.fq_name.str.contains('R2')].iterrows():
		fq_name = r2_row.fq_name
		match = fq_name.replace('R2', 'R1')
		df.loc[df.fq_name == match, 'fq_rev'] = fq_name

def importSampleSheet(df, datafolder):
	sample_sheet = glob.glob(datafolder + '*.xls*')[0]
	ref_import = pd.read_excel(sample_sheet)
	df['reference'] = None
	df['HDR'] = None
	df['guideseq'] = None

	for i, ref_row in ref_import.iterrows():
		sample_id = ref_row.Sample_ID
		df.loc[df.fq_fwd.str.contains(sample_id), 'reference'] = ref_row['WT amplicon']
		df.loc[df.fq_fwd.str.contains(sample_id), 'HDR'] = ref_row['HDR amplicon']
		df.loc[df.fq_fwd.str.contains(sample_id), 'guideseq'] = ref_row['Guide Sequence']

def getFastqName(df_row, foldertxt):
	pattern = '^(.*)_R.*'
	filename = df_row.fq_fwd.replace(foldertxt, '')
	fq_name = re.findall(pattern, filename)[0]
	return fq_name

def correctFilenames(fq_df):
	combined_fn = lambda row: os.path.basename(row.fq_fwd).replace('_R1_001.fastq.gz', '.extendedFrags.fastq.gz')
	qual_p_fn = lambda row: os.path.basename(row.fq_fwd).replace('_R1_001.fastq.gz', '_R1_qual_P.fastq.gz')
	fq_df['combined'] = fq_df.apply(combined_fn, axis=1)
	fq_df['qual_P'] = fq_df.apply(qual_p_fn, axis=1)

	combined_fn2 = lambda row: args.flashdir + '/' + row.combined
	qual_p_fn2 = lambda row: args.trimdir + '/' + row.qual_P
	fq_df['combined'] = fq_df.apply(combined_fn2, axis=1)
	fq_df['qual_P'] = fq_df.apply(qual_p_fn2, axis=1)

def crispresso_worker(row):
	fq_name = getFastqName(row, args.datafolder + '/')
	##FLASH run
	flash_cmd = '/root/CRISPResso_dependencies/bin/flash -r ' + str(readleng) + ' -f ' + str(fragleng) + ' -s ' + str(stdevfrag) + ' -o ' + fq_name + ' -d ' + args.flashdir + ' -z ' + row['fq_fwd'] + ' ' + row['fq_rev']

	##Removing adapters from reads, then passing them through quality filter, saving in new folder
	##See trimmomatic for settings. A sliding window searches along for quality score at a minimum of 
	##20 PHRED33 score, also trimming any adapter sequences. 
	trim_cmd = '/opt/conda/bin/trimmomatic SE ' + row.combined + ' ' + row.qual_P + ' ' + 'ILLUMINACLIP:' + adapterloc + ':2:30:10 LEADING:' + str(leadingqual) \
	+ ' TRAILING:' + str(trailingqual) + ' SLIDINGWINDOW:' + str(slidewindsiz) + ':' + str(slidewindqual) + ' MINLEN:' + str(minlength)	

	##CRISPResso run
	crispresso_cmd = '/opt/conda/bin/CRISPResso -r1 ' + row.qual_P + ' -a ' + row.reference + ' -e ' + row.HDR + ' -g ' + row.guideseq + ' -o ' + args.crispdir

	os.system(flash_cmd)
	os.system(trim_cmd)
	os.system(crispresso_cmd)

args = parseArgs()

# referenceseq = "96wp1sorted_sample_sheet.xlsx"
adapterloc = args.datafolder + "/TruSeq3-PE-2.fa"
readleng = "250"
fragleng = "300"
stdevfrag = "45"

################## CRISPRESSO VARIABLESE ####################
leadingqual = 3
trailingqual = 3
slidewindsiz = 4
slidewindqual = 20
minlength = 50
################## CRISPRESSO VARIABLESE ####################

def main():
	fq_df = initDataFrame(args.datafolder + '/')
	importSampleSheet(fq_df, args.datafolder + '/')
	correctFilenames(fq_df)

	if not os.path.exists(args.trimdir):
		os.makedirs(args.trimdir)

	procs = []
	for i, row in fq_df.iterrows():
		p = Process(target=crispresso_worker, args=(row,))
		procs.append(p)
		p.start()

	for p in procs:
		p.join()

	if args.s3:
		s3bucket = 'gmeixiong-bucket'
		destination = 'CRISPResso_OUT'
		s3_cmd = "aws s3 cp " + args.crispdir + ' s3://' + s3bucket + '/' + destination + ' --recursive'
		os.system(s3_cmd)
		
if __name__ == "__main__":
	main()
