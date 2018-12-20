import json
import time
import requests
import logging
logger = logging.getLogger(__name__)
logging.basicConfig(format='[%(asctime)s] [%(levelname)s] %(message)s', level=logging.DEBUG)

def print_logo(version):
	logo = "\
=============================================   \n\
  ____  _           _   ___  _____       __     \n\
 |  _ \\| |         | | |__ \\|  __ \\     / _| \n\
 | |_) | | __ _ ___| |_   ) | |__) |___| |_     \n\
 |  _ <| |/ _` / __| __| / /|  _  // _ \\  _|   \n\
 | |_) | | (_| \\__ \\ |_ / /_| | \\ \\  __/ |  \n\
 |____/|_|\\__,_|___/\\__|____|_|  \\_\\___|_|  \n\
                2018 NTOU AQUA MVIL Lafudoci    \n\
                                  Ver.%s\n\
============================================="%(version)
	print(logo)

def dump_json_file(dict, path):
	with open(path, 'w') as f:
		json.dump(dict, f)

def load_json_file(path):
	with open(path, 'r') as f:
		json_dict = json.load(f)
	return json_dict

def http_get(url):
	url = "http://"+url
	while(True):
		try:
			resp = requests.get(url)
			if resp.status_code == 200:
				return resp
			else:
				logger.warning('HTTP resp code:'+str(resp.status_code))
				time.sleep(3)
		except requests.exceptions.RequestException as err:
			logger.warning(err)

if __name__ == '__main__':
	print_logo('0.1.0')
	# http_get('rest.kegg.jp/conv/genes/ncbi-geneid:4519')