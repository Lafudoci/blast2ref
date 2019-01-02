import difflib
import json
import time
import requests
import re
import xml.etree.ElementTree as ET

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')
ncbi_api_key = config['NCBI_APIKEY']['APIKEY']


last_call_time = 0

class EutilsAPI:
	def __init__(self):
		self.base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
		self.search_api = self.base_url + "esearch.fcgi?db=%s&term=%s&retmode=json&sort=relevance"
		self.summary_api = self.base_url + "esummary.fcgi?db=%s&id=%s&retmode=json"
		self.fetch_api = self.base_url + "efetch.fcgi?db=%s&id=%s&retmode=%s&rettype=%s"
		self.citmatch_api = self.base_url + "ecitmatch.cgi?retmode=json&db=pubmed&retmode=xml&bdata=%s"
		self.link_api = self.base_url + "elink.fcgi?retmode=json&dbfrom=%s&db=%s%s"
	def request(self, url, *args):
		global last_call_time
		# add api key
		url += "&api_key=" + ncbi_api_key
		i = 0
		while(True):
			# limit request to e-utils in 8 times/sec
			while((time.time() - last_call_time) < 0.125):
				time.sleep(0.1)
			try:
				last_call_time = time.time()
				resp = requests.get(url=url)
				if (resp.status_code == 200) or (i > 5):
					return resp
				elif resp.status_code == 400:
					logger.warning(str(resp)+', retry in 10 sec.')
					time.sleep(10)
					i += 1
				else:
					logger.warning(str(resp)+', retry in 5 sec.')
					time.sleep(5)
			except requests.exceptions.RequestException as err:
				logger.warning(err)
				time.sleep(5)
		
	def search(self, db, term):
		return self.request(self.search_api%(db, term))
	def summary(self, db, idlist):
		return self.request(self.summary_api%(db, idlist))
	def fetch(self, db, idlist, retmode, rettype, *args):
		return self.request(self.fetch_api%(db, idlist, retmode, rettype))
	def citmatch(self, bdata):
		return self.request(self.citmatch_api%(bdata))
	def link(self, dbfrom, db, idlist):
		return self.request(self.link_api%(dbfrom, db, idlist))

class SequenceRecord:
	def __init__(self):
		self.acc = ''
		self.uids = []
		self.resp = 400
		self.bioseq = {}
		self.source = {}
		self.pub = {}
		self.pmids = []
		self.mesh = {}
	# pass

def xmltree_bioclass_parser(acc, root):
	# identify bioseq class
	bioseq_class = 'unknown'
	for child in root.iter():
		if child.tag == 'Bioseq-set_class':
			bioseq_class = child.get('value')
			break
		# identify non-labeled pdb-entry
		if child.tag == 'PDB-mol-id':
			if acc.startswith(child.text):
				bioseq_class = 'pdb-entry'
				break
	logger.debug('Bioseq-set_class: '+ bioseq_class)
	return bioseq_class

def xmltree_bioseq_parser(acc, root, bioseq_class):
	# Build bioseq info
	bioseq = {}
	for bqs in root.iter('Seq-entry_seq'):
		for bq in bqs:
			for child in bq.iter():
				
				# logger.debug(child.tag)

				# Create new bioseq sub dict depends on record type
				if bioseq_class == 'pdb-entry':
					# name a new bioseq sub dict by acc
					bioseq_tmp = acc
					bioseq_acc = acc
					if bioseq_tmp not in bioseq:
						bioseq[bioseq_tmp] = {}
				else:
					# name a new bioseq sub dict by temp name
					bioseq_tmp = 'bioseq_sub' + str(list(bqs).index(bq))
					if bioseq_tmp not in bioseq:
						bioseq[bioseq_tmp] = {}
				
				if child.tag == 'Textseq-id_accession':
					if child.text == bioseq[bioseq_tmp].get('accession'):
						logger.warning('Second accession id was ignored in the same bioseq: '+ child.text)
					else:
						bioseq[bioseq_tmp]['accession'] = child.text
				if child.tag == 'Textseq-id_name':
					bioseq[bioseq_tmp]['id_name'] = child.text
				if child.tag == 'Textseq-id_version':
					bioseq[bioseq_tmp]['version'] = child.text
				if child.tag == 'Seqdesc_title':
					bioseq[bioseq_tmp]['title'] = child.text
				if child.tag == 'Seqdesc_comment':
					bioseq[bioseq_tmp]['seq_comment'] = child.text
				if child.tag == 'MolInfo_biomol':
					bioseq[bioseq_tmp]['mol'] = child.get('value')
				if child.tag == 'Seq-id_gi':
					bioseq[bioseq_tmp]['gi'] = child.text
				# logger.debug(bioseq)

			# == Check bioseq after each iter ==

			# Rename bioseq sub dict from temp name to acc. If no acc, use id name as acc.
			if 'accession' in bioseq[bioseq_tmp].keys():
				bioseq_acc = bioseq[bioseq_tmp]['accession']
				bioseq[bioseq_acc] = bioseq.pop(bioseq_tmp)
			elif 'id_name' in bioseq[bioseq_tmp].keys():
				bioseq_acc = bioseq[bioseq_tmp]['id_name']
				bioseq[bioseq_acc] = bioseq.pop(bioseq_tmp)
			else:
				bioseq_acc = bioseq_tmp
				logger.warning('Acc was not found in bioseq')

			# *** Drop non-target bioseq to avoid too many results cause low mem crash ***
			if not acc.startswith(bioseq_acc):
				bioseq.pop(bioseq_acc)
				logger.info('Dropping non-target bioseq')
				continue
			else:
				# add empty string if not found in final result
				if 'version' not in bioseq[bioseq_acc]:
					bioseq[bioseq_acc]['version'] = ''
				if 'seq_comment' not in bioseq[bioseq_acc]:
					bioseq[bioseq_acc]['seq_comment'] = ''
				if 'title' not in bioseq[bioseq_acc]:
					bioseq[bioseq_acc]['title'] = bioseq[bioseq_acc]['seq_comment']
				if 'seq_comment' not in bioseq[bioseq_acc]:
					bioseq[bioseq_acc]['seq_comment'] = ''
				if 'mol' not in bioseq[bioseq_acc]:
					bioseq[bioseq_acc]['mol'] = ''
				if 'gi' not in bioseq[bioseq_acc]:
					bioseq[bioseq_acc]['gi'] = ''
	# logger.debug(bioseq)
	return bioseq

def xmltree_source_parser(root):
	# Build source info
	source = {}
	for child in root.iter('BioSource'):
		for sc in child.iter():
			if sc.tag == 'Org-ref_taxname':
				source['tax_name'] = sc.text
			if sc.tag == 'Org-ref_common':
				source['common_name'] = sc.text
			if sc.tag == 'Object-id_id':
				source['tax_id'] = sc.text
			if sc.tag == 'OrgName_lineage':
				source['lineage'] = sc.text
	return source

def xmltree_pub_parser(root):
	# Build publication info
	citart = {}
	citgen = {}
	i = 1
	j = 1
	for child in root.iter('Seqdesc_pub'):
		# parse cit-art
		ca = {}
		for pub in child.iter('Cit-art'):
			for art in pub.iter():
				# logger.debug(art.tag)
				if art.tag == 'PubMedId':
					ca['pubmed_id'] = art.text
				if art.tag == 'MedlineUID':
					ca['medline_id'] = art.text
				if art.tag == 'Cit-art_title':
					for t in art.iter():
						if t.tag =='Title_E_name':
							ca['title'] = t.text
				if art.tag == 'PubStatus' and art.attrib == 'ppublish':
					ca['ppublish'] = True
				# logger.debug(citart)
		if len(ca)>0:
			citart[i] = ca
		i += 1

		# parse cit-gen
		cg = {}
		for pub in child.iter('Cit-gen'):
			for gen in pub.iter():
				if gen.tag == 'Cit-gen_title':
					cg['title'] = gen.text
				if gen.tag == 'Cit-gen_pmid':
					cg['pubmed_id'] = gen.text
				if gen.tag == 'Cit-gen_muid':
					cg['medline_id'] = gen.text
		if len(cg)>0:
			citgen[j] = cg
		j += 1
	return {'citart': citart, 'citgen': citgen}

def xmltree_fetch_n_builder(uids):
	# Use first uid (workaround)
	if len(uids) > 0:
		uid = uids[0]
	else:
		uid = 0
	
	i = 0
	while(True):
		logger.info('Fetching XML from UID: ' + str(uid))
		api = EutilsAPI()
		resp = api.fetch('protein', uid, 'xml', 'native')
		try:
			logger.debug('Loading XML.')
			root = ET.fromstring(resp.text)
			break
		except ET.ParseError as err:
			if i < 5:
				logger.warning(err)
				logger.warning('Refetch XML in 5 sec...')
				time.sleep(5)
				i += 1
			else:
				logger.warning('XML xmltree building failed.')
				root = ET.fromstring("<?xml version=\"1.0\" encoding=\"UTF-8\"?> <Data><ERROR>Parsing error: %s</ERROR></Data>" % err)
				break

	if resp.status_code != 200:
		err = root.find('ERROR')
		logger.warning(err.text)

	return root, resp.status_code

def prot_record_parser(acc, uids):
	'''
	Input acc id and uid list
	Return record object containing acc, source, bioseq, publication info
	'''
	sr = SequenceRecord()
	sr.acc = acc
	sr.uids = uids

	root, status_code = xmltree_fetch_n_builder(sr.uids)
	sr.resp = status_code

	logger.debug('Parsing XML.')
	bio_class = xmltree_bioclass_parser(acc, root)
	sr.bioseq = xmltree_bioseq_parser(acc, root, bio_class)
	sr.source = xmltree_source_parser(root)
	sr.pub = xmltree_pub_parser(root)
	return sr

def medline_fetch(pubmed_id):
	'''
	Input PubMed id.
	Return medline string.
	'''
	api = EutilsAPI()
	resp = api.fetch('pubmed', pubmed_id, 'text', 'medline')
	logger.debug(resp.text)
	return resp.text

def medline_parser(pubmed_id, medline):
	'''
	Input medline string.
	Return field-value dict.
	'''
	md_dict = {}
	field = ''

	# Check medline format
	if not medline[0:30].strip().startswith('PMID- '+pubmed_id+'\n'):
		logger.debug('Medline was not found.')
		return md_dict

	for line in medline.split('\n'):
		if len(line)>5:
			if line[4] == '-':
				new_filed = True
				field, value = line.split('- ',1)
				field = field.strip()
				value = value.strip()
				if field in md_dict:
					if type(md_dict[field]) == list:
						md_dict[field].append(value)
					else:
						# change str to list first if it wasn't then add value
						md_dict[field] = [md_dict[field]]
						md_dict[field].append(value)
				else:
					md_dict[field] = value
			else:
				# add string as value to current field
				line = line.strip()	
				if type(md_dict[field]) is list:
					# add string 
					md_dict[field][-1] = md_dict[field][-1] + ' ' + line
				else:
					md_dict[field] = md_dict[field] + ' ' + line
	return md_dict

def record_pmid_list(sr):
	pmid_list = []
	if len(sr.pub['citart']) > 0:
		for value in sr.pub['citart'].values():
			if 'pubmed_id' in value:
				pmid_list.append(value['pubmed_id'])

	if len(sr.pub['citgen']) > 0:
		for value in sr.pub['citgen'].values():
			if 'pubmed_id' in value:
				pmid_list.append(value['pubmed_id'])

	logger.debug('PMID found:' + str(pmid_list))
	return pmid_list

def mesh_collect(sr):
	'''
	Input record obj.
	Output major & all MeSH term sets, detail dict in a dict.
	'''
	MH_set_major = set()
	MH_set_all = set()
	detail_dict = {}
	if len(sr.pmids) > 0:
		for pmid in sr.pmids:
			md = medline_fetch(pmid)
			md_dict = medline_parser(pmid, md)
			# extract MH term from md_dict
			if 'MH' in md_dict:
				for value in md_dict['MH']:
					if '*' in value:
						MH_set_major.add(value)
						# remove '*' before add into MH_set_all
						value = value.replace('*','')
						MH_set_all.add(value)
					else:
						MH_set_all.add(value)

					# build detail dict of pmid-MH
					if detail_dict.get(pmid, '') == '':
						detail_dict[pmid] = []
						detail_dict[pmid].append(value)
					else:
						detail_dict[pmid].append(value)

	if len(MH_set_major)==0:
		MH_set_major = ''
	if len(MH_set_all)==0:
		MH_set_all = ''
	if len(detail_dict)==0:
		detail_dict = {}
	mh_set_dict = {'major':list(MH_set_major),'all':list(MH_set_all), 'detail': detail_dict}
	return mh_set_dict

def uid2record(acc, uids):
	record = prot_record_parser(acc, uids)
	if len(record.bioseq) > 0:
		record.pmids = record_pmid_list(record)
		record.mesh = mesh_collect(record)
		logger.debug(acc + ' record was succesfully built.')
	else:
		# if bioseg is empty then reset record
		logger.debug(acc + ' record not found.')
		record = SequenceRecord()
	return record


def acc2uid(acc):
	api = EutilsAPI()
	resp = api.search('protein', acc)
	if resp.status_code == 200:
		i = 0
		while(True):
			uids = json.loads(resp.text)['esearchresult'].get('idlist', '')
			# logger.debug(uids)
			if uids != '':
				break
			elif i > 5:
				logger.warning('UID not found.')
				uids = ['0']
				break
			else:
				logger.warning('UID not found. Retrying...')
				time.sleep(5)
				i += 1
	else:
		uids = ['0']
	return uids

def acclink_pubmed(accs):
	i = 0
	qstring = ''
	link_dict = {}
	for acc in accs:
		qstring += '&id='+acc
	api = EutilsAPI()
	resp = api.link('protein', 'pubmed', qstring)
	if resp.status_code == 200:
		linke_result = json.loads(resp.text)
		for linkset in linke_result.get('linksets'):
			gi = linkset.get('ids')
			if len(gi) > 0:
				linkdbs = linkset.get('linksetdbs')
				for link in linkdbs:
					if link.get('linkname') == 'protein_pubmed':
						link_dict[accs[i]] = {'gi':gi[0], 'pmids':link.get('links')}
						i += 1
			else:
				i += 1
	else:
		logger.warning('ELinl request failed.')
	logger.debug(link_dict)
	return link_dict

def taxid2lineage(taxid):
	while(True):
			logger.info('Fetching lineage of taxid: ' + str(taxid))
			api = EutilsAPI()
			resp = api.fetch('taxonomy', taxid, 'xml', 'native')
			try:
				# logger.debug('Loading XML...')
				root = ET.fromstring(resp.text)
				break
			except ET.ParseError as err:
				if i < 5:
					logger.warning(err)
					logger.warning('Refetch XML in 5 sec...')
					time.sleep(5)
					i += 1
				else:
					logger.warning('XML xmltree building failed.')
					root = ET.fromstring("<?xml version=\"1.0\" encoding=\"UTF-8\"?> <Data><ERROR>Parsing error: %s</ERROR></Data>" % err)
					break
	for el in root.iter():
		if el.tag == 'Lineage':
			return el.text

def pubmed_parser(pubmed_id):
	while(True):
			logger.info('Fetching pubmed_id: ' + str(pubmed_id))
			api = EutilsAPI()
			resp = api.fetch('pubmed', pubmed_id, 'xml', 'native')
			try:
				# logger.debug('Loading XML...')
				root = ET.fromstring(resp.text)
				break
			except ET.ParseError as err:
				if i < 5:
					logger.warning(err)
					logger.warning('Refetch XML in 5 sec...')
					time.sleep(5)
					i += 1
				else:
					logger.warning('XML xmltree building failed.')
					root = ET.fromstring("<?xml version=\"1.0\" encoding=\"UTF-8\"?> <Data><ERROR>Parsing error: %s</ERROR></Data>" % err)
					break
	pub_dict = {}
	for el in root.iter():
		if el.tag == 'ISOAbbreviation':
			pub_dict['iso_j_abb'] = el.text
		if el.tag == 'ArticleTitle':
			pub_dict['article_title'] = el.text
		if el.tag == 'AbstractText':
			pub_dict['abs'] = el.text.strip()
		if el.tag == 'PubMedPubDate' and el.attrib == {'PubStatus':"accepted"}:
			pub_dict['accepted_date'] = el.find('Year').text+el.find('Month').text+el.find('Day').text

	return pub_dict

def name2gene(name_string):
	gene_info = {}
	match_score = {}
	best_gene = {}
	api = EutilsAPI()
	while(True):
		name_string = name_string.lower().replace(' ', '+')
		resp = api.search('gene',name_string)
		try:
			search_results = json.loads(resp.text).get('esearchresult',-1).get('idlist',-1)
			if search_results == -1:
				logger.warning('Could not parse esearchresult. Retrying...')
				time.sleep(3)
			elif len(search_results) == 0:
				gene_info = {}
				break
			else:
				top_ids = ','.join(search_results)
				resp = api.summary('gene', top_ids)
				summary_results = json.loads(resp.text)['result']
				return_uids = summary_results.get('uids','')
				for uid in return_uids:
					gene_info[uid] = {}
					gene_info[uid]['geneid'] = summary_results[uid]['uid']
					gene_info[uid]['symbol'] = summary_results[uid]['name']
					gene_info[uid]['des'] = summary_results[uid]['description']
				break
		except json.decoder.JSONDecodeError as err:
			logger.warning(err)
	
	# Score the similarity between name_string and returned gene names
	for gene in gene_info.values():
		score = 0
		score += difflib.SequenceMatcher(None, name_string, gene['des'].lower()).ratio()
		score += difflib.SequenceMatcher(None, name_string, gene['symbol'].lower()).ratio()
		# print(gene['des'], gene['symbol'], score)
		if score > 0.6:
			match_score[gene['geneid']] = score
	if len(match_score) > 0:
		best_id = max(match_score.keys(), key=(lambda key: match_score[key]))
		best_gene = gene_info[best_id]
	return best_gene


if __name__ == '__main__':
	acc = 'XP_025759116.1'
	uids = acc2uid(acc)
	record = uid2record(acc, uids)
	logger.info('source: ' + str(record.source))
	logger.info('bioseq: ' + str(record.bioseq))
	logger.info('pub: ' + str(record.pub))
	logger.info('pmids: ' + str(record.pmids))
	logger.info('MeSh: ' + str(record.mesh))

	# accs = ['CAO01356.1','1713245A','1JC9_A', 'AAA19454']
	# acclink_pubmed(accs)

	# logger.info(taxid2lineage(162300))

	# logger.info(pubmed_parser(23954694))

	# logger.info(name2gene('signal+transducer+and+activator+of+transcription+1'))
	# logger.info(name2gene('stat1'))
	# logger.info(name2gene('cytochrome b mitochondrion'))
	# logger.info(name2gene('sodium channel subunit beta 4'))
