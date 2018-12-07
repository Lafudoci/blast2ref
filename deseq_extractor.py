
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')
de_profile = config['DE_FILTER']['PROFILE']
max_pvalue = float(config['edgeR_PROFILE']['MAX_PValue'])
max_fdr = float(config['edgeR_PROFILE']['MAX_FDR'])
min_logcpm = float(config['edgeR_PROFILE']['MIN_logCPM'])
abs_min_logfc = float(config['edgeR_PROFILE']['ABS_MIN_logFC'])

trinity_fasta = {}

fasta_path =  r'D:\2017NGS\lbr2_trinity\Trinity.fasta'
de_result_path = r'D:\2017NGS\rsem\edgeR_results\RSEM.gene.counts.matrix.control_vs_early.edgeR.DE_results'

output_path = 'diff_fasta_all.fasta'

def read_fasta(fasta_path):
	logger.info('Loading '+ str(fasta_path))
	seq_id = ''
	with open(fasta_path, 'r') as f:
		for line in f.readlines():
			line = line.strip()
			if len(line) > 0:
				if line.startswith('>'):
					seq_id = line.split(' ')[0][1:]
					seq_len = line.split(' ')[1].split('=')[1]
					trinity_fasta[seq_id] = {'seq':''}
					trinity_fasta[seq_id]['length'] = seq_len
				else:
					trinity_fasta[seq_id]['seq']+=line

	# verify fasta seq length
	invalid_seq = []
	for seq_id, value  in trinity_fasta.items():
		if int(value['length']) != len(value['seq']):
			logger.warning('Seq length not match: '+ seq_id)
			invalid_seq.append(seq_id)
	
	# pop invalid seqs
	for invalid in invalid_seq:
		logger.warning('Dropping invalid seq: '+ seq_id)
		trinity_fasta.pop(seq_id)

	logger.info('Loaded seqs: '+ str(len(trinity_fasta)))
	return trinity_fasta

def read_de_edger():
	i = 0
	header = []
	de_dict = {}
	logger.info('Loading '+ str(de_result_path))
	with open (de_result_path, 'r') as f:
		for line in f.readlines():
			# get header
			if i == 0:
				line = line.strip()
				header = line.split('\t')
				header = ['seqid'] + header
			# build DE dict
			else:
				de = line.strip().split('\t')
				de_dict[de[0]] = {}
				j = 0
				for col in header:
					if col in ('logFC','logCPM','PValue','FDR'):
						de_dict[de[0]][col] = float(de[j])
					else:
						de_dict[de[0]][col] = de[j]
					j+=1
			i+=1
	logger.info('Loaded DE results: '+ str(len(de_dict)))
	return de_dict

def de_filter_edger(de_dict):
	filtered_de = []
	for de in de_dict.values():
		if de['logFC']>abs_min_logfc and de['logCPM']>min_logcpm and de['PValue']<max_pvalue and de['FDR']<max_fdr:
			filtered_de.append(de)
	logger.info('DE genes after filtering: '+ str(len(filtered_de)))
	return filtered_de

def de_seq_extract(fasta, de_dict, outpath):
	for de in de_dict:
		print(de)



def main():
	# read fasta
	fasta_dict = read_fasta(fasta_path)

	# read DE results
	if de_profile.lower() == 'edger':
		logger.info('DE profile: edgeR')
		de_dict = read_de_edger()
		filtered_de_dict = de_filter_edger(de_dict)
	else:
		logger.warning('Unknown DE profile: '+ de_profile)

	# extract DE seq
	de_seq_extract(fasta_dict, filtered_de_dict, output_path)

	# de_cache_edger()
	# de_tsv_edger()

if __name__ == '__main__':
	main()