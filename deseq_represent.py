import json
import eutility
import utils

from pprint import pprint

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')

de_level = config['edgeR_PROFILE']['DE_Level'].lower()

def represent(out):
	enriched_hits = utils.load_json_file(out+'_enrich.cache')
	gene_dict = utils.load_json_file(out+'_kogene.cache')
	qid_final = {}
	score = {'tax_sim':{}}
	target_taxid = config['EVALUATION_PARAMETERS']['TAXID']
	min_pident = config['EVALUATION_PARAMETERS']['MIN_PIDENT']

	target_lineage = eutility.taxid2lineage(target_taxid)
	
	for qid, hit in enriched_hits.items():
		idents = []
		pubmed_all = []
		mesh_all = []
		geneid_all = []
		for acc, values in hit.items():
			# First filter by blast identity
			if values['pident'] > min_pident:
				# make idents list for caculating avg. identity
				idents.append(float(values['pident']))
				# caculate tax simility
				if target_taxid:
					simility = tax_simility(values['lineage'], target_lineage)
					if acc not in score['tax_sim'].keys():
						score['tax_sim'][acc] = simility
				# collect pubmed & mesh terms
				pubmed_all = pubmed_all + values['pubmed']
				mesh_all = mesh_all + values['mesh_all']
				# collect geneid & kegg info
				geneid = gene_dict[acc]['gene'].get('geneid','')
				geneid_all.append(geneid)
				gene_dict[geneid] = gene_dict[acc]['gene']

		# build final qid list
		qid_final[qid] = {}

		# identity
		if len(idents)>0:
			qid_final[qid]['pident'] = sum(idents)/len(idents)
			qid_final[qid]['hits'] = len(idents)
		else:
			qid_final[qid]['pident'] = 0
			qid_final[qid]['hits'] = 0
		# geneid
		top_geneid = ''
		if len(geneid_all)>0:
			top_geneid = max(set(geneid_all), key=geneid_all.count)
		qid_final[qid] = {**qid_final[qid], **gene_dict.get(top_geneid, {})}

		# mesh
		counts_dict = mesh_counter(mesh_all)
		qid_final[qid]['mesh'] = counts_dict
		
	utils.dump_json_file(qid_final, out+'_qid_final.cache')
	utils.write_final_tsv(out+'_represent_isoform', 'w', qid_final)
	return qid_final

def tax_simility(lineage, target_lineage):
	lineage_list = lineage.strip().split(';')
	target_lineage_list = target_lineage.strip().split(';')
	match = len(set(lineage_list) & set(target_lineage_list))
	simility = float(match)/float(len(target_lineage_list))
	return simility

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

def mesh_counter(mesh_lists):
	counts_dict = {}
	term_all = []
	for term in mesh_lists:
		if '/' in term:
			term_all = term_all + term.split('/')
		else:
			term_all.append(term)
	for term in set(term_all):
		counts_dict[term] = term_all.count(term)
	return counts_dict
		

if __name__ == '__main__':
	represent('blast2ref_test10_swissprot')
	