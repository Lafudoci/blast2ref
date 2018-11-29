import eutility
import utils
import argparse

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.INFO)

version = '0.1.0'

argparse = argparse.ArgumentParser()
argparse.add_argument('-j','--job_type', help="Job types: default, novel-first, refed-first.", default = 'default')
argparse.add_argument('-br','--blast_result', help="Path to BLAST+ result output.", required=True)
argparse.add_argument('-bf','--blast_format', help="Format of BLAST+ output.", default = 'fmt6')
argparse.add_argument('-de','--de_result', help="Path to differ expression result output.")
argparse.add_argument('-df','--de_format', help="Format of differ expression result output.")
argparse.add_argument('-out','--output_file', help="Path to Blast2Ref result output.")
argp = argparse.parse_args()



def main():
	utils.print_logo(version)
	# Read blast result to hits dict
	if argp.blast_format == 'fmt6':
		logger.info('Reading blast output from fmt6 file.')
		hits_dict = utils.read_fmt6(argp.blast_result)
	else:
		logger.warning('Unknow file format.')

	# Enrich hits dict with mesh
	enrich_hits_dict = utils.hits_mesh_enrich(hits_dict)

	# Dump hits dict in tsv
	logger.info('Writting result to file.')
	utils.hits2tsv('blast-test.tsv', enrich_hits_dict)

	logger.info('Blast2Ref job done.')


if __name__ == '__main__':
	main()