import json

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')
de_profile = config['DE_FILTER']['PROFILE'].lower()
de_level = config['edgeR_PROFILE']['DE_Level'].lower()
max_pvalue = float(config['edgeR_PROFILE']['MAX_PValue'])
max_fdr = float(config['edgeR_PROFILE']['MAX_FDR'])
min_logcpm = float(config['edgeR_PROFILE']['MIN_logCPM'])
abs_min_logfc = float(config['edgeR_PROFILE']['ABS_MIN_logFC'])

trinity_fasta = {}

fasta_path =  r'D:\2017NGS\lbr2_trinity\Trinity.fasta'
de_result_path = r'D:\2017NGS\rsem\edgeR_results\RSEM.gene.counts.matrix.control_vs_early.edgeR.DE_results'

fasta_output_path = r'filtered_seq.fasta'

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

	# verify DE level
	global de_level
	for de in de_dict:
		sid = de.split('_')
		if de_level == 'gene' and sid[-1].startswith('i'):
			logger.warning('Isoform id was detected, DE level should be isoform. Overwritted config.')
			de_level = 'isoform'
			break
		if de_level == 'isoform' and sid[-1].startswith('g'):
			logger.warning('Gene id was detected, DE level should be gene. Overwritted config.')
			de_level = 'gene'
			break

	return de_dict

def de_filter_edger(de_dict):
	filtered_de_list = []
	for de in de_dict.values():
		if abs(de['logFC'])>abs_min_logfc and de['logCPM']>min_logcpm and de['PValue']<max_pvalue and de['FDR']<max_fdr:
			filtered_de_list.append(de['seqid'])
	logger.info('DE genes after filtering: '+ str(len(filtered_de_list)))
	return filtered_de_list

def de_seq_extract_edger(fasta, de_list):
	logger.info('Extracting DE seq from total fasta...')
	de_seq = {}
	if de_level == 'gene':
		for iso_id, value in fasta.items():
			gene_id = '_'.join(iso_id.split('_')[:4])
			if gene_id in de_list:
				de_seq[iso_id] = value['seq']
	elif de_level == 'isoform':
		for de in de_list:
			de_seq[de] = fasta[de]['seq']
	else:
		logger.warning('Unknown DE level: '+ de_level)
	
	logger.info('Extracted fasta seq: '+ str(len(de_seq)))

	return de_seq

def de_cache_edger(de_dict, de_list):
	cache = {}
	logger.debug('Writing DE seq cache...')
	for sid, value in de_dict.items():
		if sid in de_list:
			cache[sid] = value
	with open ('deseq.cache', 'w') as f:
		f.write(json.dumps(cache))



def write_fasta(seq_dict, path):
	logger.info('Writing DE seq to fasta file: '+ path)
	with open(path, 'w') as f:
		for sid, seq in seq_dict.items():
			f.write('>'+sid+'\n'+seq+'\n')
	logger.info('The fasta file was built.')

def deseq_extractor():
	# read fasta
	fasta_dict = read_fasta(fasta_path)

	# read DE results
	if de_profile == 'edger':
		logger.info('DE profile: edgeR')
		de_dict = read_de_edger()
		filtered_de_list = de_filter_edger(de_dict)
		de_cache_edger(de_dict, filtered_de_list)
		# de_tsv_edger()

		# extract DE seq
		de_seq = de_seq_extract_edger(fasta_dict, filtered_de_list)
		write_fasta(de_seq, fasta_output_path)

	else:
		logger.warning('Unknown DE profile: '+ de_profile)



if __name__ == '__main__':
	deseq_extractor()