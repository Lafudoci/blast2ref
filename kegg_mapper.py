import csv
import eutility
import utils

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)


def name_clean(name_string):
	# remove tax and isoform
	name = name_string.split('[')[0].strip().split('isoform')[0].strip()
	# remove specific string
	name = name.replace('-like', '')
	name = name.replace(', partial', '')
	name = name.replace('(', '')
	name = name.replace(')', '').strip()
	name = name.replace('PREDICTED:', '').strip()
	name = name.replace('-', ' ')
	name = name.replace(':', ' ')
	name = name.replace('=', ' ')
	name = name.replace(' ', ' ')
	return name

def write_acc_gene_tsv(path, acc_gene):
	colnames = ['acc','geneid','symbol','des','keggid','keggko','keggmap']
	with open (path, 'w', newline="\n") as tsvfile:
		writer = csv.DictWriter(tsvfile, fieldnames=colnames, delimiter='\t')
		writer.writeheader()
		for acc, val in acc_gene.items():
			val['acc'] = acc
			writer.writerow(val)

def kegg_id_api(operation, identifer):
	if len(identifer) == 0:
		return ''
	if operation == 'conv':
		api_url = 'rest.kegg.jp/conv/genes/'
		query_url = 'ncbi-geneid:'+ str(identifer)
	elif operation == 'ko':
		api_url = 'rest.kegg.jp/link/ko/'
		query_url = str(identifer)
	else:
		logger.warning('operation'+operation+'not found.')
	
	resp = utils.http_get(api_url+query_url)
	if resp.text.startswith(query_url):
		kegg = resp.text.split('\t')[1].strip()
	else:
		kegg = ''
	return kegg

def kegg_map_api(ko):
	if len(ko) == 0:
		return []
	pathway_list = []
	api_url = 'rest.kegg.jp/link/pathway/'
	query_url = str(ko)
	resp = utils.http_get(api_url+query_url)
	if resp.text.startswith(query_url):
		for line in resp.text.splitlines():
			path_id = line.split('\tpath:')[1].strip()
			if path_id.startswith('ko'):
				pathway_list.append(path_id)
	else:
		kegg = ''
	print(pathway_list)
	return pathway_list

def gene2kegg(gene_id):
	kegg = {}
	kegg_id = kegg_id_api('conv', gene_id)
	kegg['keggid'] = kegg_id
	kegg_ko = kegg_id_api('ko', kegg_id)
	kegg['keggko'] = kegg_ko
	kegg_pathway = kegg_map_api(kegg_ko)
	kegg['keggmap'] = kegg_pathway
	logger.debug(kegg)
	return kegg
	

def mapper(hits):
	gene_dict = {}
	query_history = {}
	acc_gene = {}
	try:
		acc_gene = utils.load_json_file('acc_gene.table')
	except FileNotFoundError:
		logger.debug('Acc-gene table not found, creating new.')

	# prepare query for name2gene
	for qid, accs in hits.items():
		for acc, values in accs.items():
			if acc not in gene_dict.keys():
				gene_dict[acc] = {'origin': name_clean(values['title'])}

	# run name2gene
	for acc, values in gene_dict.items():
		if values['origin'] != '':
			if acc not in acc_gene.keys():
				if values['origin'] not in query_history.keys():
					gene_info = eutility.name2gene(values['origin'])
					kegg_info = gene2kegg(gene_info.get('geneid', ''))
					query_history[values['origin']] = [gene_info, kegg_info]
				else:
					logger.debug('Name was already Searched.')
					gene_info = query_history[values['origin']][0]
					kegg_info = query_history[values['origin']][1]
				# merge gene & kegg dict
				acc_gene[acc] = {**gene_info,**kegg_info}
				gene_dict[acc]['gene'] = acc_gene[acc]
			else:
				logger.debug('Acc was already in table.')
				gene_dict[acc]['gene'] = acc_gene[acc]
		else:
			gene_dict[acc]['gene'] = {}
	
	utils.dump_json_file(acc_gene, 'acc_gene.table')
	write_acc_gene_tsv('acc_gene.tsv',acc_gene)

	print(gene_dict)

	return gene_dict

if __name__ == '__main__':
	hits = utils.load_json_file('(test)blast2ref_cluster2_blastx_nr_seq1_1e-6.tsv.cache')
	# hits = utils.load_json_file('blast2ref_test10.tsv.cache')
	mapper(hits)