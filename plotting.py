import utils



import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)


def kegg_map_coloring(name_prefix):
	map2ko = {}
	ko2color = {}
	deg_list = utils.load_json_file(name_prefix+'_represent_gene.cache')

	# build map2ko dict
	for deg in deg_list.values():
		if len(deg['keggko'])>0:
			ko = deg['keggko'].split(':')[1]
			for kmap in deg['keggmap']:
				if kmap not in map2ko:
					map2ko[kmap] = [ko]
				else:
					if ko not in map2ko[kmap]:
						map2ko[kmap].append(ko)
	
	# build ko2color dict
	for deg in deg_list.values():
		if len(deg['keggko'])>0:
			ko = deg['keggko'].split(':')[1]
			if ko not in ko2color:
				ko2color[ko] = [deg['logFC']]
			else:
				ko2color[ko].append(deg['logFC'])
	for ko, fcs in ko2color.items():
		ko2color[ko] = de_color_mapping(sum(fcs)/len(fcs))

	# build kegg map request string
	for kmap, kos in map2ko.items():
		req = '' + kmap
		for ko in kos:
			req += '/'+ ko +'%09'+ ko2color[ko]
		kegg_weblink_pathway(req)
	return True

def de_color_mapping(fold_change):
	if fold_change > 0:
		return 'green'
	else:
		return 'red'

def kegg_weblink_pathway(dataset_string):
	url = 'www.kegg.jp/kegg-bin/show_pathway?' + dataset_string
	resp = utils.https_get(url)
	for line in resp.text.split("\n"):
		if line.startswith("<img src=\"/tmp/mark_pathway"):
			print(line)
			img_url = 'www.kegg.jp'+ line.split("\"")[1]
			img_name = img_url.split("/")[-1]
			break

	img_data = utils.https_get(img_url).content
	with open(img_name, 'wb') as handler:
		handler.write(img_data)


if __name__ == '__main__':
	kegg_map_coloring('blast2ref_diff_fasta_cluster2_nr')
	# kegg_map_coloring('blast2ref_test10_swissprot')
