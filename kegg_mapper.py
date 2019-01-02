import csv
import eutility
import utils
import lable_crawler
import os

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)


def name_clean(name_string):
	# clean uniprot record name
	if name_string.startswith('RecName: Full='):
		name = name_string.replace('RecName: Full=', '').split(';')[0].strip()
		name = name.replace(',', '')
		name = name.replace('-', ' ')
	else:
		# remove tax and isoform
		name = name_string.split('[')[0].strip().split('isoform')[0].strip()
		# remove specific string
		name = name.replace('-like', '')
		name = name.replace(', partial', '')
		name = name.replace('(', '')
		name = name.replace(')', '').strip()
		name = name.replace('PREDICTED:', '').strip()
		name = name.replace('#', '')
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
		for acc, gene in acc_gene.items():
			val = {'acc':acc, **gene}
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
	

def mapper(name_prefix):

	if not os.path.exists(name_prefix+'_kogene.cache'):
		hits = utils.load_json_file(name_prefix+'_enrich.cache')

		gene_dict = {}
		query_history = {}
		acc_gene = {}
		try:
			logger.debug('Loading Acc-gene table')
			acc_gene = utils.load_json_file('acc_gene.table')
		except FileNotFoundError:
			logger.debug('Acc-gene table not found, creating new.')

		# temp workaroud to remove acc field in accgene table
		for acc, gene in acc_gene.items():
			if 'acc' in gene:
				print('popping')
				gene.pop('acc')

		# prepare query for name2gene
		for qid, accs in hits.items():
			for acc, values in accs.items():
				if acc not in gene_dict.keys():
					gene_dict[acc] = {'origin': name_clean(values['title'])}

		# run name2gene & gene2kegg
		i = 0
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
					i += 1
				else:
					# logger.debug('Acc was already in table.')
					gene_dict[acc]['gene'] = acc_gene[acc]
			else:
				gene_dict[acc]['gene'] = {}
			if i == 50:
				utils.dump_json_file(acc_gene, 'acc_gene.table')
				write_acc_gene_tsv('acc_gene.tsv',acc_gene)
				i = 0

		logger.debug('Saving accgene table...')
		utils.dump_json_file(acc_gene, 'acc_gene.table')
		write_acc_gene_tsv('acc_gene.tsv',acc_gene)
		utils.dump_json_file(gene_dict, name_prefix+'_kogene.cache')
	else:
		hits = utils.load_json_file(name_prefix+'_enrich.cache')
		gene_dict = utils.load_json_file(name_prefix+'_kogene.cache')

	logger.debug('Building kogene tsv...')
	# write gene & kegg info back to hits
	hits_ko = {}
	for qid, accs in hits.items():
		hits_ko[qid] = {}
		for acc, values in accs.items():
			hits_ko[qid][acc] = {**hits[qid][acc] , **gene_dict[acc]['gene']}

	lable_crawler.write_tsv(name_prefix+'_enrich_kogene','w',hits_ko)

if __name__ == '__main__':
	mapper('blast2ref_diff_fasta_cluster2_nr')

	# name = name_clean('RecName: Full=Hemoglobin subunit beta 2; AltName: Full=Beta-2-globin; AltName: Full=Hemoglobin beta-2 chain')
	# gene_info = eutility.name2gene(name)
	# print(gene_info)