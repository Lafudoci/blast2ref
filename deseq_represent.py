import json
import eutility
import utils
from collections import Counter
from pprint import pprint

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')

def represent(prefix, db_name):
	out = prefix + '_' + db_name
	enriched_hits = utils.load_json_file(out+'_enrich.cache')
	gene_dict = utils.load_json_file(out+'_kogene.cache')
	represent_iso = {}
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
				if len(values['pubmed']) > 0:
					pubmed_all += values['pubmed']
				if len(values['mesh_all']) > 0:
					mesh_all += values['mesh_all']
				# collect geneid & kegg info
				geneid = gene_dict[acc]['gene'].get('geneid','')
				geneid_all.append(geneid)
				gene_dict[geneid] = gene_dict[acc]['gene']

		# build final qid list
		represent_iso[qid] = {}

		# identity
		if len(idents)>0:
			represent_iso[qid]['pident'] = sum(idents)/len(idents)
			represent_iso[qid]['hits'] = len(idents)
		else:
			represent_iso[qid]['pident'] = 0
			represent_iso[qid]['hits'] = 0
		# geneid
		top_geneid = ''
		if len(geneid_all)>0:
			top_geneid = max(set(geneid_all), key=geneid_all.count)
		represent_iso[qid] = {**represent_iso[qid], **gene_dict.get(top_geneid, {})}

		# mesh
		counts_dict = mesh_counter(mesh_all)
		represent_iso[qid]['mesh'] = counts_dict

	utils.dump_json_file(represent_iso, out+'_represent_isoform.cache')

	# Read DEG cache
	de_dict = utils.load_json_file(prefix+'_DEG.cache')
	
	# Verify DE_level
	for qid in de_dict.keys():
		if qid.split('_')[-1].startswith('g'):
			de_level = 'gene'
			break
		elif qid.split('_')[-1].startswith('i'):
			de_level = 'isoform'
			break
		else:
			logger.warning('Can not recognize DE level.')
			de_level = ''
			return False

	# if de level is gene then represent again, and append DEG data.
	if de_level == 'gene':
		# utils.write_final_tsv(out+'_represent_isoform', 'w', represent_iso)
		represent_gene = iso2gene_represent(represent_iso)
		represent_deg = append_deg(prefix, represent_gene, de_dict)
		utils.dump_json_file(represent_deg, out+'_represent_gene.cache')
		utils.write_final_tsv(out+'_represent_gene', 'w', represent_deg)
		return True
	elif de_level == 'isoform':
		represent_deg = append_deg(prefix, represent_iso, de_dict)
		utils.dump_json_file(represent_deg, out+'_represent_isoform.cache')
		utils.write_final_tsv(out+'_represent_isoform', 'w', represent_deg)
		return True
	else:
		logger.warning('Can not recognize DE level.')

def iso2gene_represent(represent_iso):
	all_geneid = {}
	qid_geneid_counts = {}
	represent_gene = {}
	
	# build geneid table
	for qid, gene in represent_iso.items():
		if gene.get('geneid') and (gene.get('geneid') not in all_geneid):
			all_geneid[gene['geneid']] = dict(gene)
			# all_geneid[gene['geneid']].pop('qid')
			all_geneid[gene['geneid']].pop('geneid')
			all_geneid[gene['geneid']].pop('hits')
			all_geneid[gene['geneid']].pop('pident')
			all_geneid[gene['geneid']].pop('mesh')

	# counts geneid hits
	for qid, gene in represent_iso.items():
		qid_gene = qid.split('i')[0][:-1]
		if qid_gene not in qid_geneid_counts:
			qid_geneid_counts[qid_gene] = {}
		if gene['pident'] > 0:
			if gene.get('geneid'):
				if gene['geneid'] not in qid_geneid_counts[qid_gene]:
					qid_geneid_counts[qid_gene][gene['geneid']] = 0
				qid_geneid_counts[qid_gene][gene['geneid']] += gene['hits']

	# pick top geneid by hits
	for qid, geneids in qid_geneid_counts.items():
		topid_hits = Counter(geneids).most_common(1)
		if len(topid_hits)>0:
			represent_gene[qid] = {}
			represent_gene[qid]['geneid'] = topid_hits[0][0]
			represent_gene[qid]['hits'] = topid_hits[0][1]
			represent_gene[qid] = {**represent_gene[qid], **all_geneid[topid_hits[0][0]]}

	# collesct mesh set
	for qid_gene, val_gene in represent_gene.items():
		if 'mesh' not in represent_gene[qid_gene]:
			represent_gene[qid_gene]['mesh'] = Counter()
		for val_iso in represent_iso.values():
			if val_iso.get('geneid') == val_gene['geneid']:
				represent_gene[qid_gene]['mesh'] =dict(Counter(represent_gene[qid_gene]['mesh']) + Counter(val_iso['mesh']))

	# caculate pident
	for qid_gene, val_gene in represent_gene.items():
		if 'pident' not in represent_gene[qid_gene]:
			represent_gene[qid_gene]['pident'] = 0
		idents = []
		for val_iso in represent_iso.values():
			if val_iso.get('geneid') == val_gene['geneid']:
				idents.append(val_iso['pident'])
		represent_gene[qid_gene]['pident'] = sum(idents)/len(idents)

	return represent_gene

def append_deg(prefix, represent_dict, de_dict):
	represent_deg = {}
	for qid, represent in represent_dict.items():
		if qid in de_dict:
			represent_deg[qid] = {**represent_dict[qid], **de_dict[qid]}
	return represent_deg


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
	represent('blast2ref_diff_fasta_cluster2', 'nr')
	