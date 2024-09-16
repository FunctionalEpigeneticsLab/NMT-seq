import os
import glob
import pandas as pd
import argparse

wd="/staging/leuven/stg_00064/Xinran/project/snNMT_PD"
os.chdir(wd)

def get_combined_count_mtx(tsv_files, skiprows=None, split_pos=1):
	'''Iterate through each TSV file and combine the data'''
	combined_data = pd.DataFrame()
	for tsv_file in tsv_files:
		file_path = tsv_file
		#print('Reading '+file_path)
		current_data = pd.read_csv(file_path, sep='\t', skiprows=skiprows, index_col=0)
		current_data = current_data.drop(current_data.columns[0:5], axis=1)
		if 'NIST_ConsensusVector' in current_data.index:
			current_data.drop('NIST_ConsensusVector', inplace=True)
		combined_data = pd.concat([combined_data, current_data], axis=1)
	cells = [col.split('.')[split_pos].split('_S')[0] for col in combined_data.columns]
	combined_data.columns = cells
	print('Count matrix shape:',combined_data.shape)
	return combined_data

def write_combined_count_mtx(folder_pattern, formatted_date):
	ana_dir = os.path.join(wd, "RNA_"+formatted_date)
	counts_dir = os.path.join(ana_dir,"FEATURECOUNTS")
	if not os.path.exists(counts_dir):
		os.mkdir(counts_dir)

	tsv_files_exon = []
	for folder_path in glob.glob(folder_pattern):
		tsv_files_exon.extend([os.path.join(folder_path, f) for f in os.listdir(folder_path) if (f.endswith('featureCount.tsv'))])

	for t in tsv_files_exon:
		print(t)
	print('Number of exonic count matrices:',len(tsv_files_exon))

	mtx_exon = get_combined_count_mtx(tsv_files_exon)
	
	# Drop rows where the index starts with substring "ERCC"
	mtx_exon_sub = mtx_exon[~mtx_exon.index.str.startswith('ERCC')]
	mtx_exon_sub.shape
	
	tsv_files_intron = []
	for folder_path in glob.glob(folder_pattern):
		tsv_files_intron.extend([os.path.join(folder_path, f) for f in os.listdir(folder_path) if (f.endswith('featureCount.intron.tsv'))])

	for t in tsv_files_intron:
		print(t)
	print('Number of intronic count matrices:',len(tsv_files_intron))

	mtx_intron = get_combined_count_mtx(tsv_files_intron, skiprows=1, split_pos=2)

	if mtx_exon.shape[1]==mtx_intron.shape[1]:
		print('The number of cells match for the two count matrices. Good to go!')
	else:
		print('The number of cells for the two count matrices DO NOT MATCH! Please double check.')

	# Ensure they have the same shape and matching indices/columns
	aligned_exon, aligned_intron = mtx_exon.align(mtx_intron, join='outer', axis=0, fill_value=0)
	aligned_exon, aligned_intron = aligned_exon.align(aligned_intron, join='outer', axis=1, fill_value=0)
	# Merge
	merged_mtx = aligned_exon.add(aligned_intron, fill_value=0)
	print('Shape of merged count matrix:', merged_mtx.shape)

	merged_mtx.shape[0]==mtx_exon.shape[0] and merged_mtx.shape[1]==mtx_exon.shape[1]

	batches = merged_mtx.columns.to_list()
	print('Total number of batches sequenced: '+str(len(set([item.split('_')[0] for item in batches]))))
	print(set([item.split('_')[0] for item in batches]))

	merged_mtx.to_csv(os.path.join(counts_dir, "combined_featurecount.tsv"), sep='\t')
	aligned_exon.to_csv(os.path.join(counts_dir, "exon_featurecount.tsv"), sep='\t')
	aligned_intron.to_csv(os.path.join(counts_dir, "intron_featurecount.tsv"), sep='\t')

	df = pd.concat([mtx_exon_sub.sum(axis=0), mtx_intron.sum(axis=0)], axis=1)
	df.columns = ['exonic_count', 'intronic_count']
	df.to_csv(os.path.join(counts_dir, "count_exonic_intronic.tsv"), sep='\t')
	print('Merged count matrices were written to:', counts_dir)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Example program with a -d argument')
	parser.add_argument('-d', '--date', type=str, default="", help='Input date (prefix) for analysis (default: "")')
	args = parser.parse_args()
	folder_pattern = "./*/FEATURECOUNTS"
	write_combined_count_mtx(folder_pattern, args.date)
