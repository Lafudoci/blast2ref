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



if __name__ == '__main__':
	print_logo('0.1.0')