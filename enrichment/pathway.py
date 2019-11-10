from scipy import stats
from statsmodels.sandbox.stats import multicomp
import os
import pickle
import operator
from collections import defaultdict
from tfta.tfta import TFTA

import logging
logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA-pathwayEnrich')

_enrich_dir = os.path.dirname(os.path.realpath(__file__)) + '/../enrichment/data'
class PathwayEnrich():
    pathway_db = ['kegg', 'reactome', 'wikipathway']
    
    def __init__(self, pop, tfta):
        self.pop = pop
        self.tfta = tfta
        self.gene_sets = self.get_gene_sets()
        if self.gene_sets:
            logger.info('PathwayEnrich loaded gene sets.')
        else:
            logger.info('PathwayEnrich could not load gene sets.')
        
    
    def get_ora_enrich(self, study, db_str='kegg'):
        """
        Return the enriched kegg pathways
        
        parameter
        ---------------
        study: list or set of genes
        db_str: str
        """
        res = defaultdict(dict)
        try:
            res = self.ora(study, self.pop, self.gene_sets[db_str.lower()])
            return res
        except Exception as e:
            logger.error(e)
            return None
        
    def get_gene_sets(self):
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
    
    @staticmethod
    def ora(study, pop, gene_set, adjust='bonferroni', p_bonferroni=0.01, limit=30):
        """
        Over representation analysis.
        ------------------------------
        :param study: list or set, the significant gene symbols
        :param pop: list or set, the background gene symbols
        :param gene_set: dict, the functional gene sets
        :param adjust: the adjust method in the multiple tests, 
            details at https://www.statsmodels.org/0.8.0/generated/statsmodels.sandbox.stats.multicomp.multipletests.html
        :return: dict, the ORA analysis result
        """
        pop_mapped = dict()
        s_mapped = dict()
        for k, v in gene_set.items():
            sm = set(v['gene']) & set(study)
            if len(sm) > 1:
                pop_mapped[k] = set(v['gene']) & set(pop)
                s_mapped[k] = sm
        pvalues = {}
        for k in s_mapped.keys():
            pvalues[k] = stats.hypergeom.sf(len(s_mapped[k]) - 1, len(pop), len(pop_mapped[k]), len(study))
        
        #multiple testing correction
        #only consider those with pvalue < 0.01
        sig_pvalues = dict()
        for k,v in pvalues.items():
            if v < p_bonferroni:
                sig_pvalues[k] = v
        if sig_pvalues:
            _, pv, _, _ = multicomp.multipletests(list(sig_pvalues.values()), method=adjust)
            res_cor = {list(sig_pvalues.keys())[i]: pv[i] for i in range(len(list(sig_pvalues.keys())))}
        
            #sorting in ascending order according to pvalues, list of tuples
            res_sorted = sorted(res_cor.items(), key=operator.itemgetter(1))
        else:
            res_sorted = []
        #test purpose
        save_res(res_sorted)
        
        num = 1
        res = defaultdict(dict)
        if res_sorted:
            for pn, pv in res_sorted:
                if num > limit:
                    break
                if pv <= p_bonferroni:
                    res[pn]['p-bonferroni'] = pv
                    res[pn]['dblink'] = gene_set[pn]['dblink']
                    res[pn]['gene'] = s_mapped[pn]
                    num += 1
                else:
                    break
        return res

def save_res(pvalue):
    with open('pathway_pvalue.txt', 'w') as fw:
        for pn, pv in pvalue:
            fw.write(pn + '\t')
            fw.write(str(pv))
            fw.write('\n') 
        