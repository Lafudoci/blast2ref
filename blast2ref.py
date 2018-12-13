import eutility
import utils
import lable_crawler
import argparse
import os

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.INFO)

version = '0.1.0'

argparse = argparse.ArgumentParser()
argparse.add_argument('-j','--job_type', help="Job types: default, acc2ref, cache2tsv, novel-first, refed-first.", default = 'default')
argparse.add_argument('-br','--blast_result', help="Path to BLAST+ result output.")
argparse.add_argument('-bf','--blast_format', help="Format of BLAST+ output.", default = 'fmt6')
argparse.add_argument('-de','--de_result', help="Path to differ expression result output.")
argparse.add_argument('-df','--de_format', help="Format of differ expression result output.")
argparse.add_argument('-c','--cache_file', help="Path to Blast2Ref .cache file.")
argparse.add_argument('-out','--output_file', help="Path to Blast2Ref result output.")
argp = argparse.parse_args()

def main():
	utils.print_logo(version)

	if argp.job_type == 'default':
		logger.info('Job type: Blast2ref pipeline.')
		# check necessary args
		if argp.blast_result and argp.output_file:
			# Read blast result to hits dict
			if argp.blast_format == 'fmt6':
				logger.info('Reading blast output from fmt6 file.')
				hits_dict = lable_crawler.read_fmt6(argp.blast_result)
			else:
				logger.warning('Unknow file format.')
			# Enrich hits dict with mesh
			enrich_hits_dict = lable_crawler.hits_mesh_enrich(hits_dict, 500, argp.output_file)
			# Format hits dict into tsv
			logger.info('Writing result to file.')
			lable_crawler.write_tsv(argp.output_file, 'w', enrich_hits_dict)
			logger.info('Blast2Ref job done.')
		else:
			logger.warning('Error, lack of necessary args.')

	elif argp.job_type == 'acc2ref':
		logger.info('Job type: Collect Refs from BLAST fmt result.')
		# check necessary args
		# if argp.blast_result and argp.output_file:


	elif argp.job_type == 'cache2tsv':
		logger.info('Job type: Convert cache to tsv file.')
		# check necessary args
		if argp.cache_file and argp.output_file:
			try: 
				lable_crawler.cache2tsv(argp.cache_file, argp.output_file)
				logger.info('Cache2tsv job done.')
			except FileNotFoundError:
				logger.warning('Hits cache not found. Job skipping.')
		else:
			logger.warning('Error, lack of necessary args.')
	
	else:
		logger.warning('Unknow job type: '+ argp.job_type)


if __name__ == '__main__':
	main()