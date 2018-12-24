import deseq_extractor
import eutility
import kegg_mapper
import lable_crawler
import utils

import argparse
import os
import sys
import subprocess

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.INFO)

version = '0.1.0'

argparse = argparse.ArgumentParser()
argparse.add_argument('-j','--job_type', help="Job types: default, acc2ref, cache2tsv, novel-first, refed-first.", default = 'default')
argparse.add_argument('-f','--fasta_file', help="Path to fasta file.")
argparse.add_argument('-fs','--fasta_source', help="The source of fasta file.", default = 'trinity')
# argparse.add_argument('-br','--blast_result', help="Path to BLAST+ result output.")
argparse.add_argument('-bf','--blast_format', help="Format of BLAST+ output.", default = 'fmt6')
argparse.add_argument('-db','--blast_db', help="Database option for BLAST+", default = 'nr')
argparse.add_argument('-de','--de_result', help="Path to differ expression result output.")
argparse.add_argument('-dp','--de_profile', help="Profile of differ expression result output.")
argparse.add_argument('-c','--cache_file', help="Prefix of Blast2Ref .cache file.")
argparse.add_argument('-out','--output_prefix', help="Prefix of Blast2Ref results output.")
argp = argparse.parse_args()

def deseq_ex():
	if argp.fasta_file and argp.de_result and argp.de_profile and argp.output_prefix:
		if  os.path.exists(argp.output_prefix+'_filtered.fasta'):
			logger.info('Filtered fasta exists, skipping deseq extraction.')
			return True
		if argp.de_profile == 'edger':
			deseq_extractor.deseq_extractor(argp.fasta_file, argp.fasta_source, 'edger', argp.de_result, argp.output_prefix)
		return True
	else:
		logger.warning('Error, lack of necessary args.')

def run_blast():
	py = r"blast_batch_helper\blast_batch_helper.py"
	query_string = argp.output_prefix+'_filtered.fasta'
	if os.path.exists(query_string):
		out_string = argp.output_prefix +'_' + argp.blast_db + '_filtered.fmt6'
		others_string = '-task blastx-fast -evalue 1e-3 -remote'
		cmd = ["blastx", "-db", argp.blast_db, "-query", query_string, "-out", out_string, "-others", others_string]
		subprocess.run([sys.executable, py] + cmd)
		if os.path.exists(out_string):
			logger.info('BLAST job was finished.')
			return True
		else:
			logger.warning('Error, blast result not found: '+ out_string)
	else:
		logger.warning('Error, fasta not found: ' + query_string)
	

def cache2tsv():
	if argp.cache_file and argp.output_prefix:
		try: 
			lable_crawler.cache2tsv(argp.cache_file, argp.output_prefix)
			logger.info('Cache2tsv job done.')
		except FileNotFoundError:
			logger.warning('Hits cache not found. Job skipping.')
		return True
	else:
		logger.warning('Error, lack of necessary args.')

def acc2ref():
	if argp.output_prefix:
		# Read blast result to hits dict
		if argp.blast_format == 'fmt6':
			logger.info('Reading blast output from fmt6 file.')
			hits_dict = lable_crawler.read_fmt6(argp.output_prefix +'_'+ argp.blast_db + '_filtered.fmt6')
		else:
			logger.warning('Unknow file format.')
		# Enrich hits dict with mesh
		enrich_hits_dict = lable_crawler.hits_enrich(hits_dict, 500, argp.output_prefix+'_'+ argp.blast_db)
		# Format hits dict into tsv
		logger.info('Writing result to file.')
		lable_crawler.write_tsv(argp.output_prefix+'_'+ argp.blast_db, 'w', enrich_hits_dict)

		gene_dict = kegg_mapper.mapper(argp.output_prefix +'_'+ argp.blast_db)

		temp_dict = enrich_hits_dict
		for qid, accs in temp_dict.items():
			for acc, values in accs.items():
				enrich_hits_dict[qid][acc]['geneid'] = gene_dict.get(acc,{}).get('gene',{}).get('geneid',{})
				enrich_hits_dict[qid][acc]['symbol'] = gene_dict.get(acc,{}).get('gene',{}).get('symbol',{})
				enrich_hits_dict[qid][acc]['des'] = gene_dict.get(acc,{}).get('gene',{}).get('des',{})
				enrich_hits_dict[qid][acc]['keggid'] = gene_dict.get(acc,{}).get('gene',{}).get('keggid',{})
				enrich_hits_dict[qid][acc]['keggko'] = gene_dict.get(acc,{}).get('gene',{}).get('keggko',{})
				enrich_hits_dict[qid][acc]['keggmap'] = gene_dict.get(acc,{}).get('gene',{}).get('keggmap',{})

		lable_crawler.write_tsv(argp.output_prefix+'_'+ argp.blast_db+'_kegg', 'w', enrich_hits_dict)

		return True
	else:
		logger.warning('Error, lack of necessary args.')

def default():
	if argp.fasta_file and argp.de_profile and argp.output_prefix:
		# Read blast result to hits dict
		if deseq_ex() and run_blast() and acc2ref():
			return True
	else:
		logger.warning('Error, lack of necessary args.')

def main():
	utils.print_logo(version)
	job = False
	# default pipeline
	if argp.job_type == 'default':
		logger.info('Job type: Blast2ref pipeline.')
		job = default()
	# deseq_ex pipeline
	elif argp.job_type == 'deseq_ex':
		logger.info('Job type: Extract DE seqs from fasta.')
		job = deseq_ex()
	# cache2tsv pipeline
	elif argp.job_type == 'cache2tsv':
		logger.info('Job type: Convert cache to tsv file.')
		job = cache2tsv()
	elif argp.job_type == 'run_blast':
		job = run_blast()
	else:
		logger.warning('Unknow job type: '+ argp.job_type)

	if job:
		logger.info('Blast2Ref job was successfully finished.')


if __name__ == '__main__':
	main()