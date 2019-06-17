#import
import os
import logging
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
import wget
from ftplib import FTP
import gzip
import Bio.UniProt.GOA as GOA

from collections import defaultdict


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s', level=logging.INFO)
logger = logging.getLogger('GOEnrich')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/../enrichment/data'

class GOEnrich:
    def __init__(self, species=None):
        self.go = self.read_go()
        if self.go:
            logger.info('GOEnrich loaded GO data.')
        
        if not species:
            self.gaf_funcs = self.read_gaf()
        else:
            gaf_uri = '/pub/databases/GO/goa/' + species.upper() + '/goa_' + \
                       species.lower() + '.gaf.gz'
            self.gaf_funcs = self.read_gaf(gaf_uri=gaf_uri)
        if self.gaf_funcs:
            logger.info("GOEnrich loaded GAF file.")
        else:
            logger.error("GAF file doesn't exist or couldn't be loaded.")
        
        #background population
        self.pop = self.gaf_funcs.keys()
        
        self.assoc = self.get_assoc_gene_go()
    
    def read_go(self, go_obo_url=None, data_folder=None):
        if not go_obo_url:
            go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        if not data_folder:
            data_folder = _resource_dir
        
        #read GO data as dictionary
        go_obo = self.download_go(go_obo_url, data_folder)
        if go_obo:
            go = obo_parser.GODag(go_obo)
            return go
        else:
            logger.error("The GO data file doesn't exist.")
            return None
    
    def get_population(self):
        return self.population
    
    def get_annotations_keyword(self, keyword):
        """
        Generate a list of annotated genes which have the keyword in their names.
        
        parameter
        ------------
        keyword: str
        
        output
        ------------
        annots: dict
        """
        annots = {x: self.gaf_funcs[x]
               for x in self.gaf_funcs 
               if keyword.lower() in self.gaf_funcs[x]['DB_Object_Name'].lower()}
        return annots
        
    def get_annotations_genes(self, gene_list):
        annots = {x: self.gaf_funcs[x] for x in gene_list}
        return annots
    
    def read_gaf(self, gaf_uri=None, gaf_fn=None):
        if not gaf_uri:
            gaf_uri = '/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'
        if not gaf_fn:
            gaf_fn = gaf_uri.split('/')[-1]
        
        gaf_file = self.download_gaf(gaf_uri, gaf_fn)
        
        #read the annotations into a dictionary
        # File is a gunzip file, so we need to open it in this way
        with gzip.open(gaf_file, 'rt') as gaf_fp:
            gaf_funcs = {}
            # Iterate on each function using Bio.UniProt.GOA library.
            for entry in GOA.gafiterator(gaf_fp):
                gene_symbol = entry.pop('DB_Object_Symbol')
                gaf_funcs[gene_symbol] = entry
        return gaf_funcs
        
    def go_enrichment_analysis(self, study, pop=None, assoc=None, go=None, propagate_counts=True,
                               alpha=0.05, methods=None, p_bonferroni=0.01, limit=30):
        """
        GO enrichment analysis
        
        output
        ----------
        res: list of GOEnrichmentRecord
        """
        if not pop:
            pop = self.pop
        if not assoc:
            assoc = self.assoc
        if not go:
            go = self.go
        
        #methods = ["bonferroni", "sidak", "holm", "fdr"]
        #In order to reduce delay, only use bonferroni by default
        if not methods:
            methods = ['bonferroni']
        g = GOEnrichmentStudy(pop, assoc, go,
                         propagate_counts=propagate_counts,
                         alpha=alpha,
                         methods=methods)
        g_res = g.run_study(study)
        
        #only return results with pvalue < pvalue
        res = []
        num = 0
        for x in g_res:
            if num > limit:
                break
            if x.p_bonferroni <= p_bonferroni:
                res.append(x)
                num += 1
        return res
    
    def get_assoc_gene_go(self):
        assoc = defaultdict(set)
        for x in self.gaf_funcs:
            assoc[x].add(str(self.gaf_funcs[x]['GO_ID']))
        return assoc
    
    def download_go(self, go_obo_url, data_folder):
        # Check if we have the ./data directory already
        if(not os.path.isfile(data_folder)):
            # Emulate mkdir -p (no error if folder exists)
            try:
                os.mkdir(data_folder)
            except OSError as e:
                if(e.errno != 17):
                    logger.error('GOEnrich could not make the folder.')
                    return None
        else:
            logger.error('Data path (' + data_folder + ') exists as a file. '
                         'Please rename, remove or change the desired location of the data path.')
            return None
        
        # Check if the file exists already
        if(not os.path.isfile(data_folder+'/go-basic.obo')):
            go_obo = wget.download(go_obo_url, data_folder+'/go-basic.obo')
            
        else:
            go_obo = data_folder+'/go-basic.obo'
        return go_obo
        
    def download_gaf(self, gaf_uri, gaf_fn, data_folder=None):
        if not data_folder:
            data_folder = _resource_dir
        # Check if the file exists already
        gaf = os.path.join(data_folder, gaf_fn)
        if(not os.path.isfile(gaf)):
            # Login to FTP server
            ebi_ftp = FTP('ftp.ebi.ac.uk')
            ebi_ftp.login() # Logs in anonymously
            
            # Download
            with open(gaf,'wb') as fp:
                ebi_ftp.retrbinary('RETR {}'.format(gaf_uri), fp.write)
                
            # Logout from FTP server
            ebi_ftp.quit()
            
        return gaf
