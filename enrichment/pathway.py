from scipy import stats
from statsmodels.sandbox.stats import multicomp
import os
import pickle
import operator
from collections import defaultdict
#from tfta.tfta import TFTA
from utils.util import download_file_dropbox
import logging
logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA-pathwayEnrich')

_enrich_dir = os.path.dirname(os.path.realpath(__file__)) + '/../enrichment/data'

class PathwayEnrich():
    pathway_db = ['kegg', 'reactome', 'wikipathway', 'msigdb']
    disease_db = ['ctd', 'disgenet']
    disease_link = {'ctd': 'https://www.dropbox.com/s/tim3ywltd4eg7dx/ctd.pickle?dl=1',
                    'disgenet': 'https://www.dropbox.com/s/x7usg4qh3frthjw/disgenet.pickle?dl=1'}
    
    def __init__(self, pop, tfta):
        self.pop = pop
        self.tfta = tfta
        
        self.gene_sets = self.get_pathway_geneset()
        if self.gene_sets:
            logger.info('PathwayEnrich loaded pathway-genesets.')
        else:
            logger.info('PathwayEnrich could not load pathway-genesets.')
        
        #load disease geneset when it's needed
        #self.disease2gene = defaultdict(dict)
        #In order to avoid long-time delay, just load the disease-genesets in advance
        self.disease2gene = self.get_disease_geneset()
        if self.disease2gene:
            logger.info('PathwayEnrich loaded disease-genesets.')
        else:
            logger.info('PathwayEnrich could not load disease-genesets.')
        
    
    def get_ora_pathway(self, study, db_str='kegg'):
        """
        Return the enriched pathways
        
        parameter
        ---------------
        study: list or set of genes
        db_str: str
        """
        #res = defaultdict(dict)
        try:
            res = self.ora(study, self.pop, self.gene_sets[db_str.lower()])
            return res
        except Exception as e:
            logger.error(e)
            return None
            
    def get_ora_disease(self, study, db_str='ctd'):
        """
        Return the enriched diseases
        
        parameter
        ---------------
        study: list or set of genes
        db_str: str
        """
        db_str = db_str.lower()
        if db_str in self.disease_db:
            if db_str in self.disease2gene:
                gene_set = self.disease2gene[db_str]
            else:
                fn = os.path.join(_enrich_dir, db_str + '.pickle')
                url = self.disease_link[db_str]
                gene_set = self.load_pickle_file(fn, url=url)
                if gene_set:
                    self.disease2gene[db_str] = gene_set
        else:
            return None
            
        try:
            res = self.ora(study, self.pop, gene_set)
            return res
        except Exception as e:
            logger.error(e)
            return None
        
    def get_pathway_geneset(self):
        gene_sets = dict()
        filenames = dict()
        is_exist = True
        for db in self.pathway_db:
            fn = os.path.join(_enrich_dir, db + '.pickle')
            filenames[db] = fn
            if not os.path.isfile(fn):
                is_exist = False
        
        if is_exist:
            for db in self.pathway_db:
                with open(filenames[db], 'rb') as pickle_in:
                    gene_sets[db] = pickle.load(pickle_in)
        else:
            for db in self.pathway_db:
                gset = self.tfta.get_pathway_genes(db)
                if len(gset):
                    gene_sets[db] = gset
                else:
                    return None
        return gene_sets
        
    def get_disease_geneset(self):
        disease2gene = dict()
        for db in self.disease_db:
            fn = os.path.join(_enrich_dir, db + '.pickle')
            url = self.disease_link[db]
            
            geneset = self.load_pickle_file(fn, url=url)
            if geneset:
                disease2gene[db] = geneset
        return disease2gene
        
    @staticmethod
    def ora(study, pop, gene_set, adjust='bonferroni', p_bonferroni=0.01, limit=30):
        """
        Over representation analysis based on hypergeometric test.
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.hypergeom.html
        
        parameters
        -------------
        study: list or set, the significant gene symbols
        pop: list or set, the background gene symbols
        gene_set: dict, the functional gene sets
        adjust: the adjust method in the multiple tests, see details at 
        https://www.statsmodels.org/0.8.0/generated/statsmodels.sandbox.stats.multicomp.multipletests.html
        
        return
        -------------
        dict, the ORA analysis result
        """
        pop_mapped = dict()
        s_mapped = dict()
        for k, v in gene_set.items():
            sm = set(v['gene']) & set(study)
            if len(sm) >= 1:
                pop_mapped[k] = set(v['gene']) & set(pop)
                s_mapped[k] = sm
        pvalues = {}
        for k in s_mapped.keys():
            pvalues[k] = stats.hypergeom.sf(len(s_mapped[k]) - 1, len(pop), len(pop_mapped[k]), len(study))
        
        #multiple testing correction
        _, pv, _, _ = multicomp.multipletests(list(pvalues.values()), method=adjust)
        res_cor = {list(pvalues.keys())[i]: pv[i] for i in range(len(pvalues))}
        #sorting in ascending order according to pvalues, list of tuples
        res_sorted = sorted(res_cor.items(), key=operator.itemgetter(1))
        
        # only consider those with pvalue < 0.01
#         sig_pvalues = dict()
#         for k,v in pvalues.items():
#             if v < p_bonferroni:
#                 sig_pvalues[k] = v
#         if sig_pvalues:
#             _, pv, _, _ = multicomp.multipletests(list(sig_pvalues.values()), method=adjust)
#             res_cor = {list(sig_pvalues.keys())[i]: pv[i] for i in range(len(sig_pvalues))}
#         
#             #sorting in ascending order according to pvalues, list of tuples
#             res_sorted = sorted(res_cor.items(), key=operator.itemgetter(1))
#         else:
#             res_sorted = []
        
        #test purpose
        save_res(res_sorted)
        
        num = 1
        res = defaultdict(dict)
        if res_sorted:
            for pn, pv in res_sorted:
                if num > limit:
                    break
                if pv < p_bonferroni:
                    res[pn]['p-bonferroni'] = pv
                    res[pn]['dblink'] = gene_set[pn]['dblink']
                    res[pn]['gene'] = s_mapped[pn]
                    num += 1
                else:
                    break
            #If there's no enriched ones, return top 5 terms
            if not res:
                nlimit = 5
                for pn, pv in res_sorted:
                    res[pn]['p-bonferroni'] = pv
                    res[pn]['dblink'] = gene_set[pn]['dblink']
                    res[pn]['gene'] = s_mapped[pn]
                    if len(res) >= nlimit:
                        break
        return res
        
    @staticmethod
    def load_pickle_file(fn, url=None):
        if os.path.isfile(fn):
            with open(fn, 'rb') as pickle_in:
                res = pickle.load(pickle_in)
            return res
        elif url:
            #download data
            logger.info('Downloading data...')
            download_file_dropbox(url, fn)
            if os.path.isfile(fn):
                with open(fn, 'rb') as pickle_in:
                    res = pickle.load(pickle_in)
                return res
            else:
                return None
        else:
            return None
        

def save_res(pvalue):
    with open('pathway_pvalue.txt', 'w') as fw:
        for pn, pv in pvalue:
            fw.write(pn + '\t')
            fw.write(str(pv))
            fw.write('\n') 
        