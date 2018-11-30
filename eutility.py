import json
import time
import requests
import xml.etree.ElementTree as ET

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

import configparser
config = configparser.ConfigParser()
config.read('config.ini')
ncbi_api_key = config['NCBI_APIKEY']['APIKEY']

eutility_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
last_call_time = 0

class EutilsAPI:
	def __init__(self):
		self.search_api = eutility_base_url + "esearch.fcgi?db=%s&term=%s&retmode=json" + "&api_key=" + ncbi_api_key
		self.summary_api = eutility_base_url + "esummary.fcgi?db=%s&id=%s&retmode=json" + "&api_key=" + ncbi_api_key
		self.fetch_api = eutility_base_url + "efetch.fcgi?db=%s&id=%s&retmode=%s&rettype=%s" + "&api_key=" + ncbi_api_key
		self.citmatch_api = eutility_base_url + "ecitmatch.cgi?retmode=json&db=pubmed&retmode=xml&bdata=%s" + "&api_key=" + ncbi_api_key
		self.link_api = eutility_base_url + "elink.fcgi?retmode=jsondbfrom=%s&db=%s&id=%s" + "&api_key=" + ncbi_api_key
	def request(delf, url, *args):
		global last_call_time
		while(True):
			# limit request to e-utils in 8 times/sec
			while((time.time() - last_call_time) < 0.125):
				time.sleep(0.1)
			try:
				last_call_time = time.time()
				resp = requests.get(url=url)
				if str(resp) == '<Response [200]>':
					return resp
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
	# def __init__(self):
	# 	self.acc = ''
	# 	self.bioseq = {}
	# 	self.source = {}
	# 	self.pub = {}
	# 	self.mesh = {}
	pass

def xmltree_bioclass_parser(acc, root):
	# identify bioseq class
	bioseq_class = ''
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
					if child.text in bioseq[bioseq_tmp].values():
						if child.text != bioseq[bioseq_tmp]['accession']:
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

def prot_record_parser(acc):
	'''
	Input acc id
	Return record object containing acc, source, bioseq, publication info
	'''
	logger.info('Fetching XML: ' + acc)
	api = EutilsAPI()
	resp = api.fetch('protein', acc, 'xml', 'native')
	logger.debug('Loading XML.')
	root = ET.fromstring(resp.text)
	logger.debug('Parsing XML.')
	sr = SequenceRecord()

	bio_class = xmltree_bioclass_parser(acc, root)

	sr.acc = acc
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
	# logger.debug(resp.text)
	return resp.text

def medline_parser(pubmed_id, medline):
	'''
	Input medline string.
	Return field-value dict.
	'''
	md_dict = {}
	field = ''

	# Check medline format
	if not medline[0:30].strip().startswith('PMID- '+pubmed_id):
		logger.debug('Medline was not found.')
		return md_dict

	for line in medline.split('\n'):
		if len(line)>0:
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
	Output major and all MeSH term sets in a dict.
	'''
	MH_set_major = set()
	MH_set_all = set()
	if len(sr.pmid_list) > 0:
		for pmid in sr.pmid_list:
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

	if len(MH_set_major)==0:
		MH_set_major = ''
	if len(MH_set_all)==0:
		MH_set_all = ''
	mh_set_dict = {'major':MH_set_major,'all':MH_set_all}
	return mh_set_dict

def acc2record(acc):
	record = prot_record_parser(acc)
	record.pmid_list = record_pmid_list(record)
	record.mesh = mesh_collect(record)
	logger.debug(acc + ' record was succesfully built.')
	return record


# def acc2uid(acc):
# 	api = EutilsAPI()
# 	resp = api.search('protein', acc)
# 	uids = json.loads(resp.text)['esearchresult']['idlist']
# 	return uids

if __name__ == '__main__':
	# record = acc2record('1713245A')
	record = acc2record('XP_017736661.1')
	logger.info('source: ' + str(record.source['tax_name']))
	logger.info('bioseq: ' + str(record.bioseq))
	logger.info('pmid: ' + str(record.pmid_list))
	logger.info('MeSh: ' + str(record.mesh['major']))