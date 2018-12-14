import json
import eutility
import utils

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

def read_cache(path):
	with open(path, 'r') as f:
		hits = json.load(f)
	return hits

def qid_represent(enriched_hits, target_taxid, out):
	tax_sim = tax_simility(enriched_hits, target_taxid)


	# write back to hits
	for qid, hit in enriched_hits.items():
		for acc, values in hit.items():
			enriched_hits[qid][acc]['tax_sim'] = tax_sim[acc]

	utils.dump_json_file(enriched_hits, out)


def tax_simility(enriched_hits, target_taxid):
	tax_sim = {}
	target_lineage = eutility.taxid2lineage(target_taxid)
	# logger.debug('Target tx lineage: '+target_lineage)
	logger.info('Calculating lineage simility..')
	for qid, hit in enriched_hits.items():
		for acc, values in hit.items():
			lineage_list = values['lineage'].strip().split(';')
			target_lineage_list = target_lineage.strip().split(';')
			match = len(set(lineage_list) & set(target_lineage_list))
			simility = float(match)/float(len(target_lineage_list))
			
			if not tax_sim.get(acc):
				tax_sim[acc] = simility
	print(tax_sim)
	return tax_sim



def tax_counter(enriched_hits):
	tax_set = set()
	all_tax = []
	tax_counts_dict = {}
	for hit in enriched_hits.values():
		for acc in hit.values():
			tax_set.add(acc['tax_name'])
			all_tax.append(acc['tax_name'])
	for tax in tax_set:
		tax_counts_dict[tax] = all_tax.count(tax)
	return tax_counts_dict




if __name__ == '__main__':
	enriched_hits = read_cache('blast2ref_test10.tsv.cache')
	# print(enriched_hits)
	qid_represent(enriched_hits, 162300, 'blast2ref_test10.tsv.tax.cache')
	