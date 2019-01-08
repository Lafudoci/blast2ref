import json
import csv
import utils

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')

max_pvalue = 0.0
max_fdr = 0.0
min_logcpm = 0.0
abs_min_logfc = 0.0
de_level = ''


def read_trinity_fasta(fasta_path):
	trinity_fasta = {}
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

def read_de_edger(de_file):
	i = 0
	header = []
	de_dict = {}
	logger.info('Loading '+ str(de_file))
	with open (de_file, 'r') as f:
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
		if sid[-1].startswith('i'):
			de_level = 'isoform'
			break
		if sid[-1].startswith('g'):
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

def de_cache_edger(de_dict, de_list, output_prefix):
	cache = {}
	name = output_prefix + '_DEG'
	logger.debug('Writing DE seq cache: ' + name)
	for sid, value in de_dict.items():
		if sid in de_list:
			cache[sid] = value
	with open (name+'.cache', 'w') as f:
		f.write(json.dumps(cache))

	write_de_tsv(name, cache)

def write_de_tsv(name, de_dict):
	colnames = []
	for de in de_dict.values():
		for header in de.keys():
			colnames.append(header)
		break
	with open (name+'.tsv', 'w', newline="\n") as tsvfile:
		writer = csv.DictWriter(tsvfile, fieldnames=colnames, delimiter='\t')
		writer.writeheader()
		for de in de_dict.values():
			writer.writerow(de)


def write_fasta(seq_dict, prefix):
	file_name = prefix+'_filtered.fasta'
	logger.info('Writing DE seq to: ' + file_name)
	with open(file_name, 'w') as f:
		for sid, seq in seq_dict.items():
			f.write('>'+sid+'\n'+seq+'\n')
	logger.info('The filtered fasta file was built.')

def deseq_extractor(fasta_path, fasta_source, de_profile, de_file, fasta_output_prefix):
	# read fasta
	if fasta_source == 'trinity':
		fasta_dict = read_trinity_fasta(fasta_path)

	# read DE results
	if de_profile == 'edger':
		logger.info('DE profile: edgeR')
		global max_pvalue
		global max_fdr
		global min_logcpm
		global abs_min_logfc
		max_pvalue = float(config['edgeR_PROFILE']['MAX_PValue'])
		max_fdr = float(config['edgeR_PROFILE']['MAX_FDR'])
		min_logcpm = float(config['edgeR_PROFILE']['MIN_logCPM'])
		abs_min_logfc = float(config['edgeR_PROFILE']['ABS_MIN_logFC'])
		
		de_dict = read_de_edger(de_file)
		filtered_de_list = de_filter_edger(de_dict)
		de_cache_edger(de_dict, filtered_de_list, fasta_output_prefix)

		# extract DE seq
		de_seq = de_seq_extract_edger(fasta_dict, filtered_de_list)
		write_fasta(de_seq, fasta_output_prefix)

	else:
		logger.warning('Unknown DE profile: '+ de_profile)



if __name__ == '__main__':
	deseq_extractor('edger', 'GH2GH3.DE_results', 'filtered_seq_test')