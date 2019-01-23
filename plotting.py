import utils
import os



import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)


def kegg_map_coloring(name_prefix):
	deg_json = {}
	map2ko = {}
	ko2color = {}

	# make sub folder for map images
	if not os.path.exists(name_prefix+'_kegg_maps'):
		os.makedirs(name_prefix+'_kegg_maps')

	try:
		deg_json = utils.load_json_file(name_prefix+'_represent_gene.cache')
	except FileNotFoundError:
		try:
			deg_json = utils.load_json_file(name_prefix+'_represent_isoform.cache')
		except FileNotFoundError:
			logger.warning('represent_gene or _isoform cache not found.')

	if deg_json:
		# build map2ko dict
		for deg in deg_json.values():
			if deg.get('keggko'):
				# workaround
				if type(deg['keggko']) is not list:
					deg['keggko'] = [deg['keggko']]
				for keggko in deg['keggko']:
					ko = keggko.split(':')[1]
					for kmap in deg['keggmap']:
						if kmap not in map2ko:
							map2ko[kmap] = [ko]
						else:
							if ko not in map2ko[kmap]:
								map2ko[kmap].append(ko)
		# build ko2color dict
		for deg in deg_json.values():
			if deg.get('keggko'):
				for keggko in deg['keggko']:
					ko = keggko.split(':')[1]
					if ko not in ko2color:
						ko2color[ko] = [deg['hits'], deg['logFC']]
					else:
						if deg['hits'] > ko2color[ko][0]:
							ko2color[ko] = [deg['hits'], deg['logFC']]
		for ko, fcs in ko2color.items():
			ko2color[ko] = de_color_mapping(fcs[1])

		# build kegg map request string
		for kmap, kos in map2ko.items():
			req = '' + kmap
			for ko in kos:
				req += '/'+ ko +'%09'+ ko2color[ko]
			kegg_weblink_pathway(name_prefix+'_kegg_maps', req)
		return True

def de_color_mapping(fold_change):
	if fold_change > 0:
		return 'green'
	else:
		return 'red'

def kegg_weblink_pathway(path, dataset_string):
	url = 'www.kegg.jp/kegg-bin/show_pathway?' + dataset_string
	resp = utils.http_get(url, https=True)
	for line in resp.text.split("\n"):
		if line.startswith("<img src=\"/tmp/mark_pathway"):
			print(line)
			img_url = 'www.kegg.jp'+ line.split("\"")[1]
			img_name = img_url.split("/")[-1]
			break

	img_data = utils.http_get(img_url, https=True).content
	with open(path+'\\'+img_name, 'wb') as handler:
		handler.write(img_data)


if __name__ == '__main__':
	# kegg_map_coloring('blast2ref_diff_fasta_cluster2_nr')
	kegg_map_coloring('blast2ref_test10_swissprot')
