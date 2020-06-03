import pickle
import os
import logging
from utils.util import download_file_dropbox, make_folder


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/resources/'
_enrich_dir = os.path.dirname(os.path.realpath(__file__)) + '/enrichment/data'
success = make_folder(_enrich_dir)
if not success:
    logger.info('Could not create {}.'.format(_enrich_dir))


def main():
    #download tf-target db
    tf_db_file = os.path.join(_resource_dir, 'TF_target_20191224.db')
    if not os.path.exists(tf_db_file):
        logger.info('Downloading TF_target db file...')
        url = 'https://www.dropbox.com/s/gjoal1xe420je6o/TF_target_20191224.db?dl=1'
        download_file_dropbox(url, tf_db_file)
        
    #download mirna-disease db file and related files
    mir_db_file = os.path.join(_resource_dir, 'mirnaDisease.db')
    if not os.path.exists(mir_db_file):
        logger.info('Downloading miRNA disease db file...')
        url = 'https://www.dropbox.com/s/gfi0gu9gry0hkxs/mirnaDisease.db?dl=1'
        download_file_dropbox(url, mir_db_file)
    
    fn = os.path.join(_resource_dir, 'mirna_precursor_dict.pickle')
    if not os.path.exists(fn):
        logger.info('Downloading miRNA precursor file...')
        url = 'https://www.dropbox.com/s/t783u9w6ugk1g9r/mirna_precursor_dict.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    fn = os.path.join(_resource_dir, 'precursor_mirna_dict.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading miRNA precursor file...')
        url = 'https://www.dropbox.com/s/s8fxni731i0a8on/precursor_mirna_dict.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    #hgnc symbol to id mapping file
    fn = os.path.join(_resource_dir, 'hgnc_symbol_id.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading miRNA precursor file...')
        url = 'https://www.dropbox.com/s/edm9mxq8wmbztpb/hgnc_symbol_id.pickle?dl=1'
        download_file_dropbox(url, fn)
        
    #download enrichment related files
    fn = os.path.join(_enrich_dir, 'go-basic.obo')
    if not os.path.exists(fn):
        logger.info('Downloading GO file...')
        url = 'https://www.dropbox.com/s/q7zlqej49exkk52/go-basic.obo?dl=1'
        download_file_dropbox(url, fn)
        
    fn = os.path.join(_enrich_dir, 'goa_human.gaf.gz')
    if not os.path.exists(fn):
        logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/yi7bce4l5qvcfrp/goa_human.gaf.gz?dl=1'
        download_file_dropbox(url, fn)
        
    fn = os.path.join(_enrich_dir, 'gaf_funcs.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/ebbhe8atailo3ya/gaf_funcs.pickle?dl=1'
        download_file_dropbox(url, fn)
        
    #gene disease data
    fn = os.path.join(_enrich_dir, 'ctd.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/tim3ywltd4eg7dx/ctd.pickle?dl=1'
        download_file_dropbox(url, fn)
        
    fn = os.path.join(_enrich_dir, 'disgenet.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/x7usg4qh3frthjw/disgenet.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    fn = os.path.join(_enrich_dir, 'hmdd.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/3elx75k1mes22oy/hmdd.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    #pathways
    fn = os.path.join(_enrich_dir, 'kegg.pickle')
    if not os.path.exists(fn):
        logger.info('Downloading pathway files...')
        url = 'https://www.dropbox.com/s/pzmfcheldabg1dx/kegg.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    fn = os.path.join(_enrich_dir, 'msigdb.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/l1zbm5gn2iszod0/msigdb.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    fn = os.path.join(_enrich_dir, 'reactome.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/tsetv9q83c6gr2h/reactome.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    fn = os.path.join(_enrich_dir, 'wikipathway.pickle')
    if not os.path.exists(fn):
        #logger.info('Downloading GOA file...')
        url = 'https://www.dropbox.com/s/gf8qzv4i275s9le/wikipathway.pickle?dl=1'
        download_file_dropbox(url, fn)
    
    #mirna
    
    
if __name__ == "__main__":
    main()