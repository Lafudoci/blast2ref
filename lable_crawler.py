import eutility
import csv
import json

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

def read_fmt6(path, start_qid=''):
	col = []
	hits = {}
	if start_qid == '':
		start = True
	else:
		start = False

	with open (path, 'r') as f:
		for line in f.readlines():
			col = line.strip().split('\t')
			# qid = col[0]
			# sid = col[1]
			pident = col[2]
			# evalue = col[10]
			# score = col[11]
			if (start == False) and (start_qid == col[0]):
				start = True

			if start == True:
				if (col[0] not in hits):
					hits[col[0]] = {}
					hits[col[0]][col[1]] = {}
					hits[col[0]][col[1]]['pident'] = col[2]
					hits[col[0]][col[1]]['evalue'] = col[10]
					hits[col[0]][col[1]]['score'] = col[11]
				else:
					if col[1] not in hits[col[0]]:
						hits[col[0]][col[1]] = {}
						hits[col[0]][col[1]]['pident'] = col[2]
						hits[col[0]][col[1]]['evalue'] = col[10]
						hits[col[0]][col[1]]['score'] = col[11]
					else:
						# only keep best score each acc in dict
						if col[11] > hits[col[0]][col[1]]['score']:
							hits[col[0]][col[1]] = {}
							hits[col[0]][col[1]]['pident'] = col[2]
							hits[col[0]][col[1]]['evalue'] = col[10]
							hits[col[0]][col[1]]['score'] = col[11]
	return hits

def hits_enrich(hits, cache_interval, out):
	# load cache
	cache = chache2hits(out)

	i = 1
	id_tmp = {}
	for qid, accs in hits.items():
		for acc, value in accs.items():
			if cache and cache.get(qid,{}).get(acc) and cache[qid][acc].get('status'):
				logger.debug('%s, %s already in cache.'% (qid,acc))
				hits[qid][acc] = cache[qid][acc]
				id_tmp[acc] = qid
			else:
				if acc not in id_tmp:
					r_acc = acc.split('.')[0]	# There is no version in record dict
					logger.info('Enriching '+ qid +' '+ acc)
					uids = eutility.acc2uid(acc) 
					record = eutility.uid2record(acc, uids)
					hits[qid][acc]['uid'] = record.uids
					hits[qid][acc]['title'] = record.bioseq.get(r_acc, {}).get('title','')
					hits[qid][acc]['gi'] = record.bioseq.get(r_acc, {}).get('gi','')
					hits[qid][acc]['tax_id'] = record.source.get('tax_id','')
					hits[qid][acc]['tax_name'] = record.source.get('tax_name','')
					hits[qid][acc]['lineage'] = record.source.get('lineage','')
					hits[qid][acc]['pubmed'] = record.pmids
					hits[qid][acc]['mesh_major'] = record.mesh.get('major',[])
					hits[qid][acc]['mesh_all'] = record.mesh.get('all',[])
					hits[qid][acc]['mesh_detail'] = record.mesh.get('detail',{})
					hits[qid][acc]['status'] = record.resp
					if qid not in cache:
						cache[qid] = {}
					else:
						cache[qid][acc] = hits[qid][acc]
					id_tmp[acc] = qid
				else:
					logger.info('Found %s from %s'%(acc,id_tmp[acc]))
					hits[qid][acc]['uid'] = hits[id_tmp[acc]][acc].get('uid',[])
					hits[qid][acc]['title'] = hits[id_tmp[acc]][acc].get('title','')
					hits[qid][acc]['gi'] = hits[id_tmp[acc]][acc].get('gi','')
					hits[qid][acc]['tax_id'] = hits[id_tmp[acc]][acc].get('tax_id','')
					hits[qid][acc]['tax_name'] = hits[id_tmp[acc]][acc].get('tax_name','')
					hits[qid][acc]['lineage'] = hits[id_tmp[acc]][acc].get('lineage','')
					hits[qid][acc]['pubmed'] = hits[id_tmp[acc]][acc].get('pubmed',[])
					hits[qid][acc]['mesh_major'] = hits[id_tmp[acc]][acc].get('mesh_major',[])
					hits[qid][acc]['mesh_all'] = hits[id_tmp[acc]][acc].get('mesh_all',[])
					hits[qid][acc]['mesh_detail'] = hits[id_tmp[acc]][acc].get('mesh_detail',{})
					hits[qid][acc]['status'] = hits[id_tmp[acc]][acc].get('status','')
					if qid not in cache:
						cache[qid] = {}
					else:
						cache[qid][acc] = hits[qid][acc]

				if i == cache_interval:
					hits2cache(cache, out)
					i = 1
				i+=1

	hits2cache(cache, out)
	write_tsv(out+'_enrich', 'w', hits)

	return True

def chache2hits(out):
	cache = {}
	try:
		with open(out+'_enrich.cache', 'r') as f:
			cache = json.load(f)
		logger.debug('Hits cache exists. Continuing previous job.')
	except FileNotFoundError:
		logger.debug('Hits cache not found. Creating new job.')
	return cache
	

def hits2cache(hits, out):
	with open(out+'_enrich.cache', 'w') as f:
		logger.info('Writing cache')
		json.dump(hits, f)

def cache2tsv(file_name, out):
	logger.info('Loading '+ file_name + '.cache')
	with open(file_name+'_enrich.cache', 'r') as f:
		hits = json.load(f)
	logger.info('Writing tsv file...')
	write_tsv(out, 'w', hits)

def write_tsv(out, mode, hits):
	with open (out+'.tsv', mode, newline="\n") as tsvfile:
		colnames = ['qid','sid',
					'pident','evalue','score',
					'uid','title','gi',
					'tax_id','tax_name','lineage',
					'pubmed','mesh_major','mesh_all', 'mesh_detail',
					'status','geneid','symbol','des',
					'keggid','keggko','keggmap'
		]
		writer = csv.DictWriter(tsvfile, fieldnames=colnames, delimiter='\t')
		writer.writeheader()
		for qid, subj in hits.items():
			for bioseq in subj.items():
				row_dict = {'qid':qid, 'sid':bioseq[0], **bioseq[1]}
				writer.writerow(row_dict)

if __name__ == '__main__':
	hits_dict = read_fmt6('blast2ref_test10_swissprot_filtered.fmt6')
	enrich_hits_dict = hits_enrich(hits_dict, 10, 'blast2ref_test10_swissprot')