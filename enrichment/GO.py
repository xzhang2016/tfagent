#import
import os
import time
import logging
from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
import wget
from ftplib import FTP
import gzip
import Bio.UniProt.GOA as GOA
from collections import defaultdict
import pickle
import sqlite3

logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s', level=logging.INFO)
logger = logging.getLogger('TFTA-GOEnrich')

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
            self.gaf_funcs = self.read_gaf2(gaf_uri=gaf_uri)
        if self.gaf_funcs:
            logger.info("GOEnrich loaded GAF file.")
        else:
            logger.error("GAF file doesn't exist or couldn't be loaded.")
        
        #background population
        self.pop = self.gaf_funcs.keys()
        
        self.assoc = self.get_assoc_gene_go()
        
        #load go_gene.db
        db_file = os.path.join(_resource_dir, 'go_gene1.db')
        if os.path.isfile(db_file):
            self.godb = sqlite3.connect(db_file, check_same_thread=False)
            logger.info('TFTA-GOEnrich loaded go_gene database.')
        else:
            num = self.generate_go2gene_db(db_file)
            if os.path.isfile(db_file):
                self.godb = sqlite3.connect(db_file, check_same_thread=False)
            else:
                logger.error('TFTA-GOEnrich could not load go_gene database.')
                self.godb = None;
    
    def read_go(self, go_obo_url=None, data_folder=None):
        if not go_obo_url:
            go_obo_url = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
        if not data_folder:
            data_folder = _resource_dir
        
        #read GO data as dictionary
        go_obo = self.download_file_http(go_obo_url, data_folder, 'go-basic.obo')
        if go_obo:
            go = obo_parser.GODag(go_obo)
            return go
        else:
            logger.error("The GO data file doesn't exist.")
            return None
            
    def read_gaf(self, gaf_uri=None, gaf_fn=None, data_folder=_resource_dir):
        """
        Download gaf file from http://current.geneontology.org/products/pages/downloads.html
        Annotations for 2019-07-01 release
        """
        if not gaf_uri:
            gaf_uri = "http://geneontology.org/gene-associations/goa_human.gaf.gz"
        if not gaf_fn:
            gaf_fn = gaf_uri.split('/')[-1]
        
        gaf_file = self.download_file_http(gaf_uri, data_folder, gaf_fn)
        
        #read the annotations into a dictionary
        # File is a gunzip file, so we need to open it in this way
        with gzip.open(gaf_file, 'rt') as gaf_fp:
            gaf_funcs = defaultdict(list)
            # Iterate on each function using Bio.UniProt.GOA library.
            for entry in GOA.gafiterator(gaf_fp):
                gene_symbol = entry.pop('DB_Object_Symbol')
                gaf_funcs[gene_symbol].append(entry)
        return gaf_funcs
    
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
        annots = dict()
        for g, entry in self.gaf_funcs.items():
            if keyword.lower() in entry[0]['DB_Object_Name'].lower():
                annots[g] = entry[0]
        return annots
        
    def get_annotations_genes(self, gene_list):
        annots = {x: self.gaf_funcs[x] for x in gene_list}
        return annots
    
    def get_go_keyword(self, keyword):
        """
        return GO terms whose names have keyword
        """
        #check if the corresponding file exists
        data_folder = _resource_dir + '/go_files'
        fn = os.path.join(data_folder, keyword + ".pickle")
        if os.path.isfile(fn):
            with open(fn, 'rb') as pickle_in:
                go_name = pickle.load(pickle_in)
        else:
            go_name = self._go_keyword(keyword)
        return go_name
    
    def _go_keyword(self, keyword):
        """
        Generate a dict, its key is GO_id, value is GO name which contains the keyword
        """
        data_folder = _resource_dir + '/go_files'
        if not os.path.isfile(data_folder):
            # Emulate mkdir -p (no error if folder exists)
            try:
                os.mkdir(data_folder)
            except OSError as e:
                if(e.errno != 17):
                    logger.error('GOEnrich could not make the go_files folder.')
                    return None
        else:
            logger.error('Data path (' + data_folder + ') exists as a file. '
                         'Please rename, remove or change the desired location of the data path.')
            return None
            
        go_name = dict()
        if self.go:
            go_name = {go_id:self.go[go_id].name for go_id in self.go if keyword in self.go[go_id].name}
            #save to file
            fn = os.path.join(data_folder, keyword + ".pickle")
            with open(fn, 'wb') as pickle_out:
                pickle.dump(go_name, pickle_out)
        return go_name
    
    
        
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
            for entry in self.gaf_funcs[x]:
                assoc[x].add(entry['GO_ID'])
        return assoc
    
    def download_file_http(self, url, data_folder, file_name):
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
        if not os.path.isfile(os.path.join(data_folder, file_name)):
            downloaded_file = wget.download(url, os.path.join(data_folder, file_name))
            
        else:
            downloaded_file = os.path.join(data_folder, file_name)
        return downloaded_file
        
    def download_file_ftp(self, ftp_site, gaf_uri, gaf_fn, data_folder=None):
        if not data_folder:
            data_folder = _resource_dir
        # Check if the file exists already
        gaf = os.path.join(data_folder, gaf_fn)
        if(not os.path.isfile(gaf)):
            # Login to FTP server
            ebi_ftp = FTP(ftp_site)
            ebi_ftp.login() # Logs in anonymously
            
            # Download
            with open(gaf,'wb') as fp:
                ebi_ftp.retrbinary('RETR {}'.format(gaf_uri), fp.write)
                
            # Logout from FTP server
            ebi_ftp.quit()
            
        return gaf
        
    def read_gaf2(self, ftp_site=None, gaf_uri=None, gaf_fn=None):
        """
        The gaf from ebi is 6/8/16 version. We will use the one from geneontology website.
        """
        if not ftp_site:
            ftp_site = 'ftp.ebi.ac.uk'
        if not gaf_uri:
            gaf_uri = '/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz'
        if not gaf_fn:
            gaf_fn = gaf_uri.split('/')[-1]
        
        gaf_file = self.download_file_ftp(ftp_site, gaf_uri, gaf_fn)
        
        #read the annotations into a dictionary
        # File is a gunzip file, so we need to open it in this way
        with gzip.open(gaf_file, 'rt') as gaf_fp:
            gaf_funcs = defaultdict(list)
            # Iterate on each function using Bio.UniProt.GOA library.
            for entry in GOA.gafiterator(gaf_fp):
                gene_symbol = entry.pop('DB_Object_Symbol')
                gaf_funcs[gene_symbol].append(entry)
        return gaf_funcs
        
    def generate_go2gene_db(self, go_file_name):
        """
        Generate a sqlite db file which contains GO id and its associated genes
        It's too slow.
        """
        t0 = time.perf_counter()
        db_file = os.path.join(_resource_dir, go_file_name)
        conn = sqlite3.connect(db_file)
        c = conn.cursor()
        c.execute('''CREATE TABLE go2genes
             (id integer, goid text, genesymbol text)''')
        t = []
        num = 1
        for gene,entry in self.gaf_funcs.items():
            for en in entry:
                go_terms = en['GO_ID']
                #c.execute("INSERT INTO go2genes VALUES (?,?,?)", (num, go_terms, gene))
                t.append((num, go_terms, gene))
                num += 1
        c.executemany('INSERT INTO go2genes VALUES (?,?,?)', t)
        conn.commit()
        logger.info('Created go_gene.db. go2genes table has {} items.'.format(num))
        time.sleep(0.1)
        c.close()
        conn.close()
        t1 = time.perf_counter()
        logger.info('Used {} seconds to generate go-gene db file.'.format(t1-t0))
        return num
        
    def generate_go2gene_file(self):
        num = 1
        fn = os.path.join(_resource_dir, 'go2genes_gaf.txt')
        with open(fn, 'w') as fw:
            for gene,entries in self.gaf_funcs.items():
                for entry in entries:
                    fw.write(str(num) + '\t' + str(entry['GO_ID']) + '\t' + gene + '\t')
                    name_space = self.go[entry['GO_ID']].namespace
                    #or
                    #name_space = entry['Aspect']
                    fw.write(name_space + '\n')
                    num += 1
        print('write {} rows to file.'.format(num))
        return num
            
