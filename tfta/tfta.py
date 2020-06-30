#TFTA stands for TF-target agent whose task is to search for targets known to a TF and vice versa

import re
import os
import operator
import logging
import sqlite3
import numpy as np
from collections import defaultdict, Counter
import math
from indra import has_config
from indra.ontology.bio import bio_ontology
from indra.tools.expand_families import expand_agent
from utils.util import merge_dict_sum, merge_dict_list
from utils.util import download_file_dropbox
import pickle


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA')

if has_config('INDRA_DB_REST_URL') and has_config('INDRA_DB_REST_API_KEY'):
    from indra.sources.indra_db_rest import get_statements
    from indra.tools.assemble_corpus import filter_evidence_source
    CAN_CHECK_STATEMENTS = True
else:
    logger.warning("Harvard web api not specified. Cannot get evidence from literature.")
    CAN_CHECK_STATEMENTS = False
 
_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/../resources/'
_enrich_dir = os.path.dirname(os.path.realpath(__file__)) + '/../enrichment/data'

pmid_sublink = "https://www.ncbi.nlm.nih.gov/pubmed/"
 
class TFNotFoundException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
 
class TargetNotFoundException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class PathwayNotFoundException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class GONotFoundException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class miRNANotFoundException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)

class TissueNotFoundException(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)
        
class KinaseNotFoundException(Exception): pass

def _make_go_map():
    lines = open(_resource_dir + 'GO_mapping.txt', 'rt').readlines()
    go_map = defaultdict(list)
    for line in lines:
        kin, goid = line.strip().split('\t')
        go_map[kin].append(goid)
    return go_map
    
go_map = _make_go_map()

def _get_mirna_indra():
    lines = open(_resource_dir + 'mirna_indra_regulateAmount.txt', 'rt').readlines()
    mirna_indra = set()
    for line in lines:
        mirna_indra.add(line.strip())
    return mirna_indra

#mirna_indra_set = _get_mirna_indra()

#hgnc official symbol to id mapping
def _get_hgnc_genes():
    hgnc_genes = dict()
    #check if the pickle file exist
    fn = os.path.join(_resource_dir, 'hgnc_symbol_id.pickle')
    if os.path.isfile(fn):
        with open(fn, 'rb') as pickle_in:
            hgnc_genes = pickle.load(pickle_in)
    else:
        with open(_resource_dir + 'hgnc_symbol_id_20190225.txt', 'rt') as fr:
            lines = fr.readlines()
        for line in lines:
            s = line.strip().split('\t')
            hgnc_genes[s[1]] = s[0]
        #save to file
        with open(fn, 'wb') as pickle_out:
            pickle.dump(hgnc_genes, pickle_out)
    return hgnc_genes

hgnc_symbol_id = _get_hgnc_genes()
hgnc_genes_set = set(hgnc_symbol_id.keys())

#gene expression threshold
EXP_THR = 1.5

class TFTA:
    def __init__(self):
        #Load TF_target database
        self.tfdb = self.load_db()
        if self.tfdb:
            self.tfdb.row_factory = sqlite3.Row
        
        self.ldd = self.load_ldd_db()
        if self.ldd:
            self.ldd.row_factory = sqlite3.Row
        
        #for exclusive query
        self.tissue_gene_exclusive = defaultdict(set)
        self.trans_factor = self.tf_set()
        self.mirna = self.mirna_set()
        
    def __del__(self):
        self.tfdb.close()

    def Is_tf_target(self,tf_name,target_name):
        """
        Return True if the tf regulates the target, and False if not
        """
        dbname = ''
        if self.tfdb is not None:
            #check if target_name in the database
            t = (target_name,tf_name)
            res = self.tfdb.execute("SELECT DISTINCT dbnames FROM CombinedDB WHERE Target = ? "
                                    "AND TF = ?", t).fetchone()
            if res:
                dbname = res['dbnames']
                return True,dbname
        return False,dbname
        
    def Is_gene_onto(self, go_name, gene_name):
        """
        Return true if the gene is in the go category of go_name
        """
        if go_name in ['transcription factor', 'tf']:
            if not self.trans_factor:
                if self.tfdb is not None:
                    t = (gene_name, )
                    res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor "
                                            "WHERE tf = ?", t).fetchall()
                    if res:
                        return True
            else:
                if gene_name in self.trans_factor:
                    return True
        else:
            try:
                goids = go_map[go_name]
            except KeyError:
                raise GONotFoundException
            if self.tfdb is not None:
                for go in goids:
                    t = (go, gene_name)
                    res = self.tfdb.execute("SELECT * FROM go2Genes WHERE termId = ? AND geneSymbol = ? ", t).fetchall()
                    if res:
                        return True
        return False
        
    def find_gene_onto(self, go_name, gene_names):
        """
        Return the genes which are in the category of go_name
        """
        go_genes = []
        if go_name in ['transcription factor', 'tf']:
            if not self.trans_factor:
                if self.tfdb is not None:
                    res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                    self.trans_factor = set([r[0] for r in res])
            go_genes = list(self.trans_factor & set(gene_names))
        else:
            try:
                goids = go_map[go_name]
            except KeyError:
                raise GONotFoundException
            if self.tfdb is not None:
                for go in goids:
                    t = (go, )
                    res = self.tfdb.execute("SELECT geneSymbol FROM go2Genes WHERE termId = ? ", t).fetchall()
                    go_genes.extend([r[0] for r in res])
                go_genes = list(set(go_genes) & set(gene_names))
        return go_genes

    def Is_tf_target_tissue(self,tf_name,target_name,tissue_name):
        """
        Return True if the tf regulates the target in a given tissue, and False if not
        """
        if self.tfdb is not None:
            regstr = '%' + tissue_name + '%'
            t = (tf_name, target_name, regstr)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                    "WHERE TF = ? AND Target = ? AND Tissue LIKE ? ", t).fetchall()
            if res:
                return True
            
        return False

    def find_tfs(self,target_names):
        """
        Return TFs regulating all the given targets
        """
        #query
        dbname = dict()
        if self.tfdb is not None:
            t = (target_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT TF,dbnames FROM CombinedDB "
                                    "WHERE Target = ? ", t).fetchall()
            if res:
                tf_names = set([r[0] for r in res])
                for r in res:
                    dbname[(r[0],target_names[0])] = r[1]
            else:
                raise TargetNotFoundException 
            if len(target_names)>1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT TF,dbnames FROM CombinedDB "
                                            "WHERE Target = ? ", t).fetchall()
                    if res:
                        tf_names = tf_names & set([r[0] for r in res])
                        for r in res:
                            dbname[(r[0],target_names[i])] = r[1]
                    else:
                        raise TargetNotFoundException
            #make the results consistent in both db and TF list
            if len(tf_names):
                #tf_names.sort()
                tf_names = tf_names.intersection(self.trans_factor)
        else:
            tf_names = set()
        return tf_names,dbname
        
    def find_common_tfs(self,target_names):
        """
        Return TFs regulating the given target list as well as the regulated targets by each TF
        """
        tf_targets = defaultdict(list)
        counts = defaultdict(int)
        if self.tfdb is not None:
            target_names = list(set(target_names))
            thr = max(2, math.ceil(len(target_names)/2))
            for target_name in target_names:
                t = (target_name,)
                res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
                                        "WHERE Target = ? ", t).fetchall()
                if res:
                    tfs = [r[0] for r in res]
                    for tf in tfs:
                        tf_targets[tf].append(target_name)
                        counts[tf] += 1
            if len(tf_targets):
                max_count = max(counts.values())
                if max_count >= thr:
                    #only consider the TFs which regulate at least thr targets
                    for tf,c in counts.items():
                        if c < thr:
                            del tf_targets[tf]
                elif max_count >= 2:
                    #consider the TFs which regulate at least 2 targets
                    for tf, c in counts.items():
                        if c < 2:
                            del tf_targets[tf]
                else:
                    raise TFNotFoundException
            else:
                raise TFNotFoundException
        return tf_targets

    def find_tfs_tissue(self,target_names,tissue_name):
        """
        Return TFs regulating the targets in a given tissue
        """
        #query
        tf_names = []
        if self.tfdb is not None:
            regstr = '%' + tissue_name + '%'
            t = (target_names[0],regstr)
            res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                    "WHERE Target = ? AND Tissue LIKE ? ", t).fetchall()
            if res:
                tf_names = self.trans_factor.intersection(set([r[0] for r in res]))
                if not tf_names:
                    return []
            else:
                return []
            if len(target_names) > 1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],regstr)
                    res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                            "WHERE Target = ? AND Tissue LIKE ? ", t).fetchall()
                    if res:
                        tf_names = tf_names & set([r[0] for r in res])
                        if not tf_names:
                            return []
                    else:
                        return []
        return tf_names

    def find_pathway_db_gene(self, dbsource, gene_names, fmembers=None):
        """
        return pathway information for given dbsource and gene_names
        """
        pathwayName = dict()
        dblink = dict()
        fnum = 0
        if self.tfdb is not None:
            pathlist = []
            if gene_names:
                ids = self.find_pathwayID_genes(gene_names)
                if ids:
                    pathlist.extend(ids)
                else:
                    raise PathwayNotFoundException
                
            if fmembers:
                fnum = len(fmembers)
                ids = self.find_pathwayID_families(fmembers)
                if ids:
                    pathlist.extend(ids)
                else:
                    raise PathwayNotFoundException
            
            #Take intersection
            id_count = Counter(pathlist)
            for id in id_count:
                if id_count[id] == (len(gene_names) + fnum):
                    t = (id, dbsource)
                    res = self.tfdb.execute("SELECT pathwayName,dblink FROM pathwayInfo "
                                            "WHERE Id = ? AND source LIKE ?", t).fetchone()
                    if res:
                        pathwayName[id] = res['pathwayName']
                        dblink[id] = res['dblink']
            if not pathwayName:
                raise PathwayNotFoundException
        return pathwayName,dblink

    def find_gene_pathway(self, pathway_names):
        """
        return genes related to pathway_name
        """
        pathwayName = dict()
        genelist = dict()
        pw_link = dict()
        if self.tfdb is not None:
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        pathwayName[r[0]] = r[1]
                        pw_link[r[0]] = r[2]
            if pathwayName:
                for pthID in pathwayName:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                           "WHERE pathwayID = ? ORDER BY genesymbol", t).fetchall()
                    genes = [r[0] for r in res1]
                    genelist[pthID] = genes
            else:
                raise PathwayNotFoundException
        return pathwayName,genelist,pw_link

    def find_tf_pathway(self, pathway_names):
        """
        Return TFs within the given pathway
        """
        pathwayName = dict()
        dblink = dict()
        tflist = dict()
        if self.tfdb is not None:
            #tf_set = self.get_tf_set()
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                        "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        pathwayName[r[0]] = r[1]
                        dblink[r[0]] = r[2]
                
            if pathwayName:
                pathwayId = pathwayName.keys()
                for pthID in pathwayId:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    tfs = self.trans_factor.intersection(set([r[0] for r in res1]))
                    if tfs:
                        tflist[pthID] = tfs
            if not len(tflist):
                raise PathwayNotFoundException
        return pathwayName,tflist,dblink
        
    def find_tf_pathway_id(self, ids):
        """
        Return TFs within the given pathway
        """
        pathwayName = dict()
        dblink = dict()
        tflist = dict()
        if self.tfdb is not None:
            for id in ids:
                t = (id,)
                res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                tfs = self.trans_factor.intersection(set([r[0] for r in res1]))
                if tfs:
                    tflist[id] = tfs
                    res = self.tfdb.execute("SELECT pathwayName,dblink FROM pathwayInfo "
                                            "WHERE Id = ? ", t).fetchone()
                    pathwayName[id] = res['pathwayName']
                    dblink[id] = res['dblink']
            if not len(tflist):
                raise PathwayNotFoundException
        return pathwayName,tflist,dblink
        
    def find_kinase_pathway(self,pathway_names):
        """
        Return kinases within the given pathway
        """
        pathwayName = dict()
        dblink = dict()
        kinaselist = dict()
        if self.tfdb is not None:
            t = ('kinase',)
            res = self.tfdb.execute("SELECT DISTINCT geneSymbol FROM go2Genes "
                                        "WHERE termId = ? ", t).fetchall()
            kinases = set([r[0] for r in res])
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                        "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        pathwayName[r[0]] = r[1]
                        dblink[r[0]] = r[2]
            if pathwayName:
                pathwayId = pathwayName.keys()
                for pthID in pathwayId:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    kins = list(kinases & set([r[0] for r in res1]))
                    if kins:
                        kinaselist[pthID] = kins
            else:
                raise PathwayNotFoundException
            if not len(kinaselist):
                raise PathwayNotFoundException
        return pathwayName, kinaselist, dblink

    def find_pathway_genes(self,gene_names, fmembers=None):
        """
        return pathways having given genes
        
        parameter
        -----------
        gene_names: list of gene symbols
        fmembers: dict, family name as key and list of Agent for members as value
        """
        pathwayName = dict()
        dblink = dict()
        fnum = 0
        if self.tfdb is not None:
            pathlist=[]
            if gene_names:
                ids = self.find_pathwayID_genes(gene_names)
                if ids:
                    pathlist.extend(ids)
                else:
                    raise PathwayNotFoundException
            #Take OR operation on family members
            if fmembers:
                fnum = len(fmembers)
                ids = self.find_pathwayID_families(fmembers)
                if ids:
                    pathlist.extend(ids)
                else:
                    raise PathwayNotFoundException
            #Take intersection
            id_count = Counter(pathlist)
            for id in id_count:
                if id_count[id] == (len(gene_names) + fnum):
                    t = (id,)
                    res = self.tfdb.execute("SELECT pathwayName,dblink FROM pathwayInfo "
                                            "WHERE Id = ? ", t).fetchone()
                    pathwayName[id] = res['pathwayName']
                    dblink[id] = res['dblink']
            if not pathwayName:
                raise PathwayNotFoundException
        return pathwayName,dblink
    
    def find_pathwayID_genes(self, gene_names):
        ids = []
        for gene_name in gene_names:
            t = (gene_name,)
            res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                     "WHERE genesymbol = ? ", t).fetchall()
            if res1:
                ids.extend([r[0] for r in res1])
            else:
                return None
        return ids
    
    def find_pathwayID_families(self, fmembers):
        """
        Take OR operation for family members.
        
        parameter
        ----------
        fmembers: dict, family name as key and list of Agent for members as value
        """
        ids = []
        for f in fmembers:
            fid = set()
            for m in fmembers[f]:
                t = (m.name, )
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    fid = fid.union(set([r[0] for r in res1]))
            if fid:
                ids.extend(list(fid))
            else:
                return None
        return ids
        

    def Is_pathway_gene(self, pathway_names, gene_names, fmembers=None):
        """
        Return pathways which contain the given genes and whose name contain substring of pathway_name
        """
        fpname = []
        fdblink = []
        if self.tfdb is not None:
            pn = dict()
            dl = dict()
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        pn[r[0]] = r[1]
                        dl[r[0]] = r[2]
                    
            if pn:
                for pthID in pn:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    genes = [r[0] for r in res1]
                    #check if genes overlap with all gene_names and any member in each family
                    num = len(set(gene_names) & set(genes))
                    if fmembers:
                        if all([num == len(gene_names), _check_overlap(set(genes), fmembers)]):
                            fpname.append(pn[pthID])
                            fdblink.append(dl[pthID])
                    elif num == len(gene_names):
                        fpname.append(pn[pthID])
                        fdblink.append(dl[pthID])
    
            if not len(fpname):
                raise PathwayNotFoundException
        return fpname, fdblink

    def find_pathway_gene_keyword(self,gene_names, keyword, fmembers=None):
        """
        return pathways having given genes and some information in pathway name
        """
        pathwayName = dict()
        dblink = dict()
        fnum = 0
        if self.tfdb is not None:
            pathlist = []
            if gene_names:
                ids = self.find_pathwayID_genes(gene_names)
                if ids:
                    pathlist.extend(ids)
                else:
                    raise PathwayNotFoundException
                
            #Take OR operation on family members
            if fmembers:
                fnum = len(fmembers)
                ids = self.find_pathwayID_families(fmembers)
                if ids:
                    pathlist.extend(ids)
                else:
                    raise PathwayNotFoundException
                        
            #intersection
            id_count = Counter(pathlist)
            regstr = '%' + keyword + '%'
            for id in id_count:
                if id_count[id] == (len(gene_names) + fnum):
                    t = (id, regstr)
                    res = self.tfdb.execute("SELECT pathwayName,dblink FROM pathwayInfo "
                                            "WHERE Id = ? AND pathwayName LIKE ? ", t).fetchone()
                    if res:
                        pathwayName[id] = res['pathwayName']
                        dblink[id] = res['dblink']
            
            if not len(dblink):
                raise PathwayNotFoundException
        return pathwayName, dblink

    def find_pathway_keyword(self, keyword):
        """
        return pathways containing the given keyword
        """
        pathwayName = dict()
        dblink = dict()
        if self.tfdb is not None:
            regstr = '%' + keyword + '%'
            t = (regstr,)
            res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ORDER BY pathwayName", t).fetchall()
            if res:
                for r in res:
                    pathwayName[r[0]] = r[1]
                    dblink[r[0]] = r[2]
            else:
                raise PathwayNotFoundException		
        
        return pathwayName,dblink

    def find_tf_keyword(self, keyword_name):
        """
        return pathways containing the given keyword and related tfs 
        """
        tflist = dict()
        newpathwayName = dict()
        newdblink = dict()
        if self.tfdb is not None:
            regstr = '%' + keyword_name + '%'
            t = (regstr,)
            res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
            if res:
                pathwayId = [r[0] for r in res]
                pathwayName = [r[1] for r in res]
                dblink = [r[2] for r in res]
                #search tfs
                for pn,pthID,lk in zip(pathwayName,pathwayId,dblink):
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    genes = [r[0] for r in res1]
                    #intersection with TF list
                    genes = set(genes).intersection(self.trans_factor)
                    if genes:
                        tflist[pthID] = genes
                        newpathwayName[pthID] = pn
                        newdblink[pthID] = lk
            else:
                raise PathwayNotFoundException
       
        if not len(newdblink):
            raise PathwayNotFoundException
        return newpathwayName,tflist,newdblink

    def find_common_pathway_genes(self, gene_names, fmembers=None, limit=30):
        """
        For a given gene list, find the pathways containing some of the genes,
        and the corresponding genes contained in each of the pathways
        
        parameter
        -----------
        gene_names: list of genes
        fmembers: dict, family name as key, the list of Agent for members as value
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        genes = defaultdict(list)
        fgenes = defaultdict(list)
        counts = defaultdict(int)
        if self.tfdb is not None:
            if gene_names:
                gene_names = list(set(gene_names))
                for gene_name in gene_names:
                    t = (gene_name,)
                    res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                             "WHERE genesymbol = ? ", t).fetchall()
                    if res1:
                        for r in res1:
                            genes[r[0]].append(gene_name)
                            counts[r[0]] += 1

            if fmembers:
                #OR operation on members of one family
                gene1,count1 = self.get_pathwayID_families(fmembers)
                #logger.info('tfta:find_common_pathway_genes:gene1=' + str(gene1))
                if gene1:
                    genes = merge_dict_sum(gene1, genes)
                if count1:
                    counts = merge_dict_sum(count1, counts)
            #logger.info('tfta:find_common_pathway_genes:after: genes=' + str(genes))
            
            if len(genes):
                max_count = max(counts.values())
                sorted_counts = sorted(counts.items(),key=lambda kv: kv[1], reverse=True)
                if max_count >= 2:
                    for pth, ct in sorted_counts:
                        if ct >= 2:
                            num += 1
                            if num > limit:
                                break
                            fgenes[pth] = genes[pth]
                        else:
                            break
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
            if len(fgenes):
                for pth in fgenes:
                    t = (pth,)
                    res = self.tfdb.execute("SELECT pathwayName, dblink FROM pathwayInfo "
                                            "WHERE Id = ? ", t).fetchone()
                    pathwayName[pth] = res['pathwayName']
                    dblink[pth] = res['dblink']
            else:
                raise PathwayNotFoundException
        return pathwayName,dblink,fgenes
    
    def get_pathwayID_families(self, fmembers, pathIDs=None):
        gene = defaultdict(list)
        count = defaultdict(int)
        for f in fmembers:
            fid = set()
            for m in fmembers[f]:
                t = (m.name,)
                #logger.info('get_pathwayID_families: m.name=' + m.name)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    fid.union(set([r[0] for r in res1]))
            if fid:
                if pathIDs:
                    fid = fid.intersection(pathIDs)
                for id in fid:
                    gene[id].append(f)
                    count[id] += 1
        return gene,count
    
    def find_common_pathway_genes_keyword(self, gene_names, keyword, fmembers=None, limit=30):
        """
        For a given gene list and keyword, find the pathways containing some of the genes,
        and return the corresponding given genes in each of the pathways 
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        counts = defaultdict(int)
        genes = defaultdict(list)
        fgenes = dict()
        if self.tfdb is not None:
            if gene_names:
                gene_names = list(set(gene_names))
                for gene_name in gene_names:
                    t = (gene_name,)
                    res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                             "WHERE genesymbol = ? ", t).fetchall()
                    if res1:
                        for r in res1:
                            genes[r[0]].append(gene_name)
                            counts[r[0]] += 1
            if fmembers:
                #OR operation on members of one family
                gene1,count1 = self.get_pathwayID_families(fmembers)
                #logger.info('tfta:find_common_pathway_genes_keyword:gene1=' + str(gene1))
                if gene1:
                    genes = merge_dict_sum(gene1, genes)
                if count1:
                    counts = merge_dict_sum(count1, counts)
            
            if len(genes):
                max_count = max(counts.values())
                sorted_counts = sorted(counts.items(),key=lambda kv: kv[1], reverse=True)
                if max_count >= 2:
                    regstr = '%' + keyword + '%'
                    for pth, ct in sorted_counts:
                        if ct >= 2:
                            t = (pth, regstr)
                            res = self.tfdb.execute("SELECT pathwayName, dblink FROM pathwayInfo "
                                            "WHERE Id = ? AND pathwayName LIKE ?", t).fetchone()
                            if res:
                                pathwayName[pth] = res['pathwayName']
                                dblink[pth] = res['dblink']
                                fgenes[pth] = genes[pth]
                                num += 1
                                if num > limit:
                                    break
                        else:
                            break
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
        return pathwayName,dblink,fgenes
        
    def find_common_pathway_genes_keyword2(self, gene_names, keyword, fmembers=None, limit=30):
        """
        For a given gene list and keyword, find the pathways containing some of the genes,
        and return the corresponding given genes in each of the pathways
        Check pathway name first, then check genes. This method aims to handle large gene list.
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        counts = defaultdict(int)
        genes = dict()
        fgenes = dict()
        temp = defaultdict(dict)
        if self.tfdb is not None:
            #find pathways whose names contain the keyword
            regstr = '%' + keyword + '%'
            t = (regstr, )
            res = self.tfdb.execute("SELECT Id, pathwayName, dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ?", t).fetchall()
            if res:
                for r in res:
                    temp[r[0]]['name'] = r[1]
                    temp[r[0]]['dblink'] = r[2]
            
            #find the pathways containing the given genes in the selected pathways 
            if gene_names:
                gene_names = set(gene_names)
                for pth in temp.keys():
                    t = (pth,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    ovg = gene_names.intersection(set([r[0] for r in res1]))
                    if len(ovg) > 1:
                        genes[pth] = ovg
                        counts[pth] = len(ovg)
            if fmembers:
                #OR operation on members of one family
                gene1,count1 = self.get_pathwayID_families(fmembers, pathIDs=set(temp.keys()))
                #logger.info('tfta:find_common_pathway_genes_keyword:gene1=' + str(gene1))
                if gene1:
                    genes = merge_dict_sum(gene1, genes)
                if count1:
                    counts = merge_dict_sum(count1, counts)
            
            if len(genes):
                max_count = max(counts.values())
                sorted_counts = sorted(counts.items(),key=lambda kv: kv[1], reverse=True)
                if max_count >= 2:
                    for pth, ct in sorted_counts:
                        if ct >= 2:
                            pathwayName[pth] = temp[pth]['name']
                            dblink[pth] = temp[pth]['dblink']
                            fgenes[pth] = genes[pth]
                            num += 1
                            if num > limit:
                                break
                        else:
                            break
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
        return pathwayName,dblink,fgenes
        
    def find_common_pathway_genes_db(self, gene_names, db_name, fmembers=None, limit=30):
        """
        For a given gene list and db name, find the pathways containing at least two of 
        the genes, and return the corresponding given genes in each of the pathways 
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        counts = defaultdict(int)
        genes = defaultdict(list)
        fgenes = dict()
        if self.tfdb is not None:
            if gene_names:
                gene_names = list(set(gene_names))
                for gene_name in gene_names:
                    t = (gene_name,)
                    res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                    if res1:
                        for r in res1:
                            genes[r[0]].append(gene_name)
                            counts[r[0]] += 1
            if fmembers:
                #OR operation on members of one family
                gene1,count1 = self.get_pathwayID_families(fmembers)
                #logger.info('tfta:find_common_pathway_genes_db:gene1=' + str(gene1))
                if gene1:
                    genes = merge_dict_sum(gene1, genes)
                if count1:
                    counts = merge_dict_sum(count1, counts)
            
            if len(genes):
                max_count = max(counts.values())
                sorted_counts = sorted(counts.items(),key=lambda kv: kv[1], reverse=True)
                if max_count >= 2:
                    for pth, ct in sorted_counts:
                        if ct >= 2:
                            t = (pth, db_name)
                            res = self.tfdb.execute("SELECT pathwayName, dblink FROM pathwayInfo "
                                            "WHERE Id = ? AND source LIKE ?", t).fetchone()
                            if res:
                                pathwayName[pth] = res['pathwayName']
                                dblink[pth] = res['dblink']
                                fgenes[pth] = genes[pth]
                                num += 1
                                if num > limit:
                                    break
                        else:
                            break
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
        return pathwayName, dblink, fgenes
        
    def find_common_pathway_genes_db2(self, gene_names, db_name, fmembers=None, limit=30):
        """
        For a given gene list and db name, find the pathways containing at least two of 
        the genes, and return the corresponding given genes in each of the pathways.
        Check pathway source first, then check genes. This method aims to handle large gene list.
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        counts = defaultdict(int)
        genes = defaultdict(list)
        fgenes = dict()
        temp = defaultdict(dict)
        if self.tfdb is not None:
            t = (db_name,)
            res = self.tfdb.execute("SELECT Id, pathwayName, dblink FROM pathwayInfo "
                                    "WHERE source LIKE ?", t).fetchall()
            if res:
                for r in res:
                    temp[r[0]]['name'] = r[1]
                    temp[r[0]]['dblink'] = r[2]
            
            if gene_names:
                gene_names = set(gene_names)
                for pth in temp.keys():
                    t = (pth,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                         "WHERE pathwayID= ? ", t).fetchall()
                    ovg = gene_names.intersection(set([r[0] for r in res1]))
                    if len(ovg) >= 2:
                        genes[pth] = ovg
                        counts[pth] = len(ovg)
            if fmembers:
                #OR operation on members of one family
                gene1,count1 = self.get_pathwayID_families(fmembers, pathIDs=set(temp.keys()))
                #logger.info('tfta:find_common_pathway_genes_db:gene1=' + str(gene1))
                if gene1:
                    genes = merge_dict_sum(gene1, genes)
                if count1:
                    counts = merge_dict_sum(count1, counts)
            
            if len(genes):
                max_count = max(counts.values())
                sorted_counts = sorted(counts.items(),key=lambda kv: kv[1], reverse=True)
                if max_count >= 2:
                    for pth, ct in sorted_counts:
                        if ct >= 2:
                            pathwayName[pth] = temp[pth]['name']
                            dblink[pth] = temp[pth]['dblink']
                            fgenes[pth] = genes[pth]
                            num += 1
                            if num > limit:
                                break
                        else:
                            break
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
        return pathwayName, dblink, fgenes

    def find_common_pathway_genes_db3(self, gene_names, db_name, fmembers=None, limit=30):
        """
        For a given gene list and db name, find the pathways containing at least two of 
        the genes, and return the corresponding given genes in each of the pathways.
        Check pathway source first, ignore familiy, then check genes and when it find the limit number of
        pathways it will return instead of checking all the pathways. 
        This method aims to handle large gene list.
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        genes = defaultdict(list)
        temp = defaultdict(dict)
        numgene = dict()
        if self.tfdb is not None:
            t = (db_name,)
            res = self.tfdb.execute("SELECT Id, pathwayName, dblink, numGenes FROM pathwayInfo "
                                    "WHERE source LIKE ?", t).fetchall()
            if res:
                for r in res:
                    temp[r[0]]['name'] = r[1]
                    temp[r[0]]['dblink'] = r[2]
                    numgene[r[0]] = r[3]
            if numgene:
                gene_names = set(gene_names)
                sorted_numgene = sorted(numgene.items(),key=lambda kv: kv[1], reverse=True)
                for pth,numg in sorted_numgene:
                    t = (pth,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID= ? ", t).fetchall()
                    ovg = gene_names.intersection(set([r[0] for r in res1]))
                    if len(ovg) > 2:
                        genes[pth] = ovg
                        #counts[pth] = len(ovg)
                        pathwayName[pth] = temp[pth]['name']
                        dblink[pth] = temp[pth]['dblink']
                        num += 1
                        if num > limit:
                            break
            else:
                raise PathwayNotFoundException
        return pathwayName, dblink, genes

    def find_common_pathway_genes_db4(self, gene_names, db_name, fmembers=None, limit=30):
        """
        For a given gene list and db name, find the pathways containing at least two of 
        the genes, and return the corresponding given genes in each of the pathways.
        Check pathway source first, ignore familiy, then check genes and when it find the limit number of
        pathways it will return instead of checking all the pathways. 
        This method aims to handle large gene list.
        """
        num = 0
        pathwayName = dict()
        dblink = dict()
        genes = defaultdict(list)
        gene_names = set(gene_names)
        if self.tfdb is not None:
            t = (db_name,)
            res = self.tfdb.execute("SELECT Id, pathwayName, dblink FROM pathwayInfo "
                                    "WHERE source LIKE ? ORDER BY numGenes DESC", t).fetchall()
            if res:
                for r in res:
                    t = (r[0],)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID= ? ", t).fetchall()
                    ovg = gene_names.intersection(set([r1[0] for r1 in res1]))
                    if len(ovg) > 2:
                        genes[r[0]] = ovg
                        pathwayName[r[0]] = r[1]
                        dblink[r[0]] = r[2]
                        num += 1
                        if num > limit:
                            break
            else:
                raise PathwayNotFoundException
        return pathwayName, dblink, genes

    def find_targets(self,tf_names, fmembers=None):
        """
        Return Targets regulated by the tf or tf list
        """
        dbname = dict()
        target_names = []
        if self.tfdb is not None:
            if tf_names:
                t = (tf_names[0],)
                res = self.tfdb.execute("SELECT DISTINCT Target,dbnames FROM CombinedDB "
                                        "WHERE TF = ? ", t).fetchall()
                if res:
                    target_names = [r[0] for r in res]
                    for r in res:
                        dbname[(tf_names[0],r[0])] = r[1]
                else:
                    raise TFNotFoundException
                if len(tf_names) > 1:
                    for i in range(1,len(tf_names)):
                        t = (tf_names[i],)
                        res = self.tfdb.execute("SELECT DISTINCT Target,dbnames FROM CombinedDB "
                                                "WHERE TF = ? ", t).fetchall()
                        if res:
                            target_names = list(set(target_names) & set([r[0] for r in res]))
                            for r in res:
                                dbname[(tf_names[i],r[0])] = r[1]
                        else:
                            raise TFNotFoundException
            #For families, not consider dbname for now
            if fmembers:
                for f in fmembers:
                    ftarget = set()
                    for m in fmembers[f]:
                        t = (m.name,)
                        res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                                "WHERE TF = ? ", t).fetchall()
                        if res:
                            ftarget = ftarget.union(set([r[0] for r in res]))
                    if ftarget:
                        target_names = list(set(target_names) & ftarget)
                    else:
                        raise TFNotFoundException
        return target_names,dbname

    def find_targets_tissue(self,tf_names, tissue_name):
        """
        Return Targets regulated by the tf list in a given tissue
        """
        target_names = set()
        if self.tfdb is not None:
            regstr = '%' + tissue_name + '%'
            t = (tf_names[0], regstr)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                    "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
            if res:
                target_names = set([r[0] for r in res])
            else:
                return []
                
            if len(tf_names) > 1:
                for i in range(1,len(tf_names)):
                    t = (tf_names[i],regstr)
                    res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                            "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
                    if res:
                        target_names = set(target_names) & set([r[0] for r in res])
                    else:
                        return []
        return target_names

    def find_genes_GO_tf(self, go_name, tf_names):
        """
        Return GO terms and genes regulated by tf_names
        """
        res_go_ids = []
        res_go_types = []
        res_go_names = []
        res_go_genes = dict()
        if self.tfdb is not None:
            t = (tf_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                    "WHERE TF = ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                raise TFNotFoundException
                
            if (len(tf_names)>1):
                for i in range(1,len(tf_names)):
                    t = (tf_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                            "WHERE TF = ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        raise TFNotFoundException
                                    
            regstr = '%' + go_name + '%'
            t = (regstr,)
            res = self.tfdb.execute("SELECT * FROM goInfo "
                                    "WHERE goName LIKE ? ", t).fetchall()
            #print('res=%s' % res)
            if not res:
                raise GONotFoundException
            else:
                record_ids = [r[0] for r in res]
                go_ids = [r[1] for r in res]
                go_types = [r[3] for r in res]
                go_names = [r[2] for r in res]
                
                #search genes
                for i in range(len(record_ids)):
                    t = (record_ids[i],)
                    res1 = self.tfdb.execute("SELECT DISTINCT geneSymbol FROM go2Genes "
                                    "WHERE termId = ? ", t).fetchall()
                    tmp = list(set(target_names) & set([r[0] for r in res1]))
                    if len(tmp):
                        res_go_ids.append(go_ids[i])
                        res_go_types.append(go_types[i])
                        res_go_names.append(go_names[i])
                        tmp.sort()
                        res_go_genes[go_ids[i]] = tmp        
        return res_go_ids,res_go_types,res_go_names,res_go_genes

    def find_genes_GO_tf2(self, go_id, tf_names):
        """
        Return GO terms and genes regulated by tf_names
        """
        res_go_ids = []
        res_go_types = []
        res_go_names = []
        res_go_genes = dict()
        if self.tfdb is not None:
            t = (tf_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                    "WHERE TF = ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                raise TFNotFoundException
                
            if (len(tf_names)>1):
                for i in range(1,len(tf_names)):
                    t = (tf_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                            "WHERE TF = ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        raise TFNotFoundException
                        
            #regstr = '%' + go_name + '%'
            t = (go_id,)
            res = self.tfdb.execute("SELECT * FROM goInfo "
                                    "WHERE goId = ? ", t).fetchall()
            if not res:
                raise GONotFoundException
            else:
                record_id = [r[0] for r in res]
                go_ids = [r[1] for r in res]
                go_names = [r[2] for r in res]
                go_types = [r[3] for r in res]
                #search genes
                for i in range(len(record_id)):
                    t = (record_id[i],)
                    res1 = self.tfdb.execute("SELECT DISTINCT geneSymbol FROM go2Genes "
                                             "WHERE termId = ? ", t).fetchall()
                    tmp = list(set(target_names) & set([r[0] for r in res1]))
                    if len(tmp):
                        res_go_ids.append(go_ids[i])
                        res_go_types.append(go_types[i])
                        res_go_names.append(go_names[i])
                        tmp.sort()
                        res_go_genes[go_ids[i]] = tmp
        return res_go_ids,res_go_types,res_go_names,res_go_genes

    def Is_miRNA_target(self, miRNA_name_dict, target_name):
        """
        Return True if the miRNA regulates the target, and False if not;
        also return evidence for provenance support
        query example: Does miR-20b-5p target STAT3?
        """
        expr = defaultdict(list)
        supt = defaultdict(list)
        pmid = defaultdict(list)
        miRNA_mis = dict()
        miRNA_name = list(miRNA_name_dict.keys())[0]
        if self.tfdb is not None:
            t = (miRNA_name, target_name)
            res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ? "
                                     "AND target = ? ", t).fetchall()
            if res:
                for r in res:
                    expr[target_name].append(r[3])
                    supt[target_name].append(r[4])
                    pmid[target_name].append(str(r[5]))
                return True, expr, supt, pmid, miRNA_mis
            else:
                #check if miRNA_name in the database
                t = (miRNA_name,)
                res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ? ", t).fetchall()
                if not res:
                    miRNA_mis[miRNA_name] = miRNA_name_dict[miRNA_name]
        return False, expr, supt, pmid, miRNA_mis
        
    def Is_miRNA_target_strength(self, miRNA_name_dict, target_name, strength):
        """
        Return True if the miRNA regulates the target, and False if not;
        also return evidence for provenance support
        query example: Does miR-20b-5p target STAT3?
        """
        expr = defaultdict(list)
        supt = defaultdict(list)
        pmid = defaultdict(list)
        miRNA_mis = dict()
        miRNA_name = list(miRNA_name_dict.keys())[0]
        if self.tfdb is not None:
            t = (miRNA_name, target_name, '%Weak%')
            if strength == 'strong':
                res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ? "
                                     "AND target = ? AND supportType NOT LIKE ?", t).fetchall()
            else:
                res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ? "
                                     "AND target = ? AND supportType LIKE ?", t).fetchall()
            if res:
                for r in res:
                    expr[target_name].append(r[3])
                    supt[target_name].append(r[4])
                    pmid[target_name].append(str(r[5]))
                return True, expr, supt, pmid, miRNA_mis
            else:
                #check if miRNA_name in the database
                t = (miRNA_name,)
                res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ? ", t).fetchone()
                if not res and (not miRNA_name[-1] in ['P', 'p']):
                    miRNA_mis[miRNA_name] = miRNA_name_dict[miRNA_name]
        return False, expr, supt, pmid, miRNA_mis
        
    def find_miRNA_target(self, target_names):
        """
        Return miRNAs regulating all the given targets
        example: What microRNAs target STAT3?
        
        parameter
        ----------
        target_names: list
        """
        miRNAs = set()
        expr = defaultdict(list)
        supt = defaultdict(list)
        pmid = defaultdict(list)
        if self.tfdb is not None:
            t = (target_names[0],)
            res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE target = ? ", t).fetchall()
            if res:
                for r in res:
                    miRNAs.add(r[1])
                    expr[(r[1],target_names[0])].append(r[3])
                    supt[(r[1],target_names[0])].append(r[4])
                    pmid[(r[1],target_names[0])].append(str(r[5]))
            else:
                raise TargetNotFoundException
            
            if len(target_names)>1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],)
                    res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE target = ? ", t).fetchall()
                    if res:
                        miRNAs = miRNAs & set([r[1] for r in res])
                        for r in res:
                            expr[(r[1],target_names[i])].append(r[3])
                            supt[(r[1],target_names[i])].append(r[4])
                            pmid[(r[1],target_names[i])].append(str(r[5]))
                    else:
                        raise TargetNotFoundException
        return miRNAs,expr,supt,pmid
    
    def find_miRNA_target_strength(self, target_names, evidence_strength):
        """
        Return miRNAs regulating all the given targets
        example: What microRNAs target STAT3?
        
        parameter
        ----------
        target_names: list
        """
        miRNAs = set()
        expr = defaultdict(list)
        supt = defaultdict(list)
        pmid = defaultdict(list)
        if self.tfdb is not None:
            if evidence_strength == 'strong':
                t = (target_names[0],'%Weak%')
                res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE target = ? AND supportType NOT LIKE ?", t).fetchall()
            else:
                t = (target_names[0], '%Weak%')
                res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE target = ? AND supportType LIKE ?", t).fetchall()
            if res:
                for r in res:
                    miRNAs.add(r[1])
                    expr[(r[1],target_names[0])].append(r[3])
                    supt[(r[1],target_names[0])].append(r[4])
                    pmid[(r[1],target_names[0])].append(str(r[5]))
            else:
                raise TargetNotFoundException
            
            if len(target_names)>1:
                if evidence_strength == 'strong':
                    for i in range(1,len(target_names)):
                        t = (target_names[i],'%Weak%')
                        res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE target = ? AND supportType NOT LIKE ?", t).fetchall()
                        if res:
                            miRNAs = miRNAs & set([r[1] for r in res])
                            for r in res:
                                expr[(r[1],target_names[i])].append(r[3])
                                supt[(r[1],target_names[i])].append(r[4])
                                pmid[(r[1],target_names[i])].append(str(r[5]))
                        else:
                            raise TargetNotFoundException
                else:
                    for i in range(1,len(target_names)):
                        t = (target_names[i],'%Weak%')
                        res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE target = ? AND supportType LIKE ?", t).fetchall()
                        if res:
                            miRNAs = miRNAs & set([r[1] for r in res])
                            for r in res:
                                expr[(r[1],target_names[i])].append(r[3])
                                supt[(r[1],target_names[i])].append(r[4])
                                pmid[(r[1],target_names[i])].append(str(r[5]))
                        else:
                            raise TargetNotFoundException
        return miRNAs,expr,supt,pmid
    
    def find_target_miRNA(self, miRNA_name_dict):
        """
        Return Targets regulated by the given miRNAs
        example: What genes does miR-20b-5p target?
        
        parameter
        ------------
        miRNA_name_dict: dict
        """
        miRNA_mis = dict()
        target_names = set()
        expr = defaultdict(list)
        supt = defaultdict(list)
        pmid = defaultdict(list)
        miRNA_names = list(miRNA_name_dict.keys())
        if self.tfdb is not None:
            t = (miRNA_names[0],)
            res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
            if res:
                for r in res:
                    target_names.add(r[2])
                    expr[(miRNA_names[0],r[2])].append(r[3])
                    supt[(miRNA_names[0],r[2])].append(r[4])
                    pmid[(miRNA_names[0],r[2])].append(str(r[5]))
            else:
                #raise miRNANotFoundException
                miRNA_mis[miRNA_names[0]] = miRNA_name_dict[miRNA_names[0]]
                return [],miRNA_mis,[],[],[]
            if len(miRNA_names) > 1:
                for i in range(1, len(miRNA_names)):
                    t = (miRNA_names[i],)
                    res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
                    if res:
                        target_names = target_names & set([r[2] for r in res])
                        for r in res:
                            expr[(miRNA_names[i],r[2])].append(r[3])
                            supt[(miRNA_names[i],r[2])].append(r[4])
                            pmid[(miRNA_names[i],r[2])].append(str(r[5]))
                    else:
                        #raise miRNANotFoundException
                        miRNA_mis[miRNA_names[i]] = miRNA_name_dict[miRNA_names[i]]
                        return [],miRNA_mis,[],[],[]
                    if not target_names:
                        break
        return target_names,miRNA_mis,expr,supt,pmid
        
    def find_target_miRNA_strength(self, miRNA_name_dict, evidence_strength):
        """
        Return Targets regulated by the given miRNAs
        example: What genes does miR-20b-5p target with strong evidence?
        
        parameters
        ------------
        miRNA_name_dict: dict
        evidence_strength: string
        """
        miRNA_mis = dict()
        target_names = set()
        expr = defaultdict(list)
        supt = defaultdict(list)
        pmid = defaultdict(list)
        miRNA_names = list(miRNA_name_dict.keys())
        if self.tfdb is not None:
            if evidence_strength == 'strong':
                t = (miRNA_names[0], '%Weak%')
                res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType NOT LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        target_names.add(r[2])
                        expr[(miRNA_names[0],r[2])].append(r[3])
                        supt[(miRNA_names[0],r[2])].append(r[4])
                        pmid[(miRNA_names[0],r[2])].append(str(r[5]))
                else:
                    #raise miRNANotFoundException
                    if self.check_mirna_not_in_db(miRNA_names[0]):
                        miRNA_mis[miRNA_names[0]] = miRNA_name_dict[miRNA_names[0]]
                    return [],miRNA_mis,[],[],[]
                    
                if len(miRNA_names)>1:
                    for i in range(1, len(miRNA_names)):
                        t = (miRNA_names[i], '%Weak%')
                        res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType NOT LIKE ?", t).fetchall()
                        if res:
                            target_names = target_names.intersection(set([r[2] for r in res]))
                            for r in res:
                                expr[(miRNA_names[i],r[2])].append(r[3])
                                supt[(miRNA_names[i],r[2])].append(r[4])
                                pmid[(miRNA_names[i],r[2])].append(str(r[5]))
                        else:
                            #raise miRNANotFoundException
                            if self.check_mirna_not_in_db(miRNA_names[i]):
                                miRNA_mis[miRNA_names[i]] = miRNA_name_dict[miRNA_names[i]]
                            return [],miRNA_mis,[],[],[]
                        if not target_names:
                            break
            else:
                t = (miRNA_names[0], '%Weak%')
                res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        target_names.add(r[2])
                        expr[(miRNA_names[0],r[2])].append(r[3])
                        supt[(miRNA_names[0],r[2])].append(r[4])
                        pmid[(miRNA_names[0],r[2])].append(str(r[5]))
                else:
                    #raise miRNANotFoundException
                    if self.check_mirna_not_in_db(miRNA_names[0]):
                        miRNA_mis[miRNA_names[0]] = miRNA_name_dict[miRNA_names[0]]
                    return [],miRNA_mis,[],[],[]
                if len(miRNA_names)>1:
                    for i in range(1, len(miRNA_names)):
                        t = (miRNA_names[i], '%Weak%')
                        res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType LIKE ?", t).fetchall()
                        if res:
                            target_names = target_names.intersection(set([r[2] for r in res]))
                            for r in res:
                                expr[(miRNA_names[i],r[2])].append(r[3])
                                supt[(miRNA_names[i],r[2])].append(r[4])
                                pmid[(miRNA_names[i],r[2])].append(str(r[5]))
                        else:
                            #raise miRNANotFoundException
                            if self.check_mirna_not_in_db(miRNA_names[i]):
                                miRNA_mis[miRNA_names[i]] = miRNA_name_dict[miRNA_names[i]]
                            return [],miRNA_mis,[],[],[]
                        if not target_names:
                            break
        return target_names,miRNA_mis,expr,supt,pmid
    
    def check_mirna_not_in_db(self, miRNA_name):
        if self.tfdb is not None:
            t = (miRNA_name,)
            res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
            if res:
                return False
            else:
                return True
    
    def find_tf_miRNA(self, miRNA_names):
        """
        Return TFs regulated by the given miRNAs
        example: what transcription factors does miR-124-3p regulate?
        miRNA_names: dict
        """
        target_names,miRNA_mis = self.find_target_miRNA(miRNA_names)
        tf_names = []
        if not miRNA_mis:
            if target_names:
                if not self.trans_factor:
                    if self.tfdb is not None:
                        res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                        self.trans_factor = set([r[0] for r in res])
                tf_names = list(set(target_names) & self.trans_factor)
        return tf_names,miRNA_mis
    
    def find_evidence_strength(self, strength, limit=30):
        """
        
        Return experimental evidence for certain evidence
        example: What constitutes strong evidence?
        
        parameter
        --------------
        strength: str
        """
        experiments = []
        if self.tfdb is not None:
            t = ('%Weak%',)
            if strength == 'strong':
                res = self.tfdb.execute("SELECT DISTINCT experiments FROM mirnaInfo "
                                    "WHERE supportType NOT LIKE ? ", t).fetchall()
            elif strength == 'weak':
                res = self.tfdb.execute("SELECT DISTINCT experiments FROM mirnaInfo "
                                    "WHERE supportType LIKE ? ", t).fetchall()
            else:
                res = self.tfdb.execute("SELECT DISTINCT experiments FROM mirnaInfo ").fetchall()
            if res:
                for r in res:
                    experiments.append(r[0])
                    if len(experiments) > limit:
                        break
        return experiments
    
    def find_evidence_miRNA_target(self, miRNA_name_dict, target_name):
        """
        
        Return experimental evidence that the given miRNA regulates the target
        example: What is the evidence that miR-148a-3p targets DNMT1?
        
        parameter
        --------------
        miRNA_name_dict: dict
        target_name: str
        """
        miRNA_mis = {}
        experiments = []
        support_types = []
        pmid_link = []
        miRNA_name = list(miRNA_name_dict.keys())[0]
        if self.tfdb is not None:
            t = (miRNA_name, target_name)
            res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND target = ? ", t).fetchall()
            if res:
                experiments = [r[3] for r in res]
                support_types = [r[4] for r in res]
                pmid_link = [pmid_sublink+str(r[5]) for r in res]
            else:
                #check if miRNA_name in the database
                t = (miRNA_name,)
                res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                        "WHERE mirna LIKE ? ", t).fetchall()
                if not res:
                    miRNA_mis[miRNA_name] = miRNA_name_dict[miRNA_name]
        return experiments, support_types, pmid_link, miRNA_mis
        
    def find_evidence_miRNA_target_strength(self, miRNA_name_dict, target_name, strength):
        """
        
        Return experimental evidence that the given miRNA regulates the target
        example: What is the strong (or weak) evidence that miR-148a-3p targets DNMT1?
        
        parameter
        --------------
        miRNA_name_dict: dict
        target_name: str
        strength: str
        """
        miRNA_mis = {}
        experiments = []
        support_types = []
        pmid_link = []
        miRNA_name = list(miRNA_name_dict.keys())[0]
        if self.tfdb is not None:
            t = (miRNA_name, target_name, '%Weak%')
            if strength == 'strong':
                res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ?"
                                    " AND target = ? AND supportType NOT LIKE ?", t).fetchall()
            else:
                res = self.tfdb.execute("SELECT * FROM mirnaInfo WHERE mirna LIKE ?"
                                    " AND target = ? AND supportType LIKE ?", t).fetchall()
            if res:
                experiments = [r[3] for r in res]
                support_types = [r[4] for r in res]
                pmid_link = [pmid_sublink+str(r[5]) for r in res]
            else:
                #check if miRNA_name in the database
                t = (miRNA_name,)
                res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                        "WHERE mirna LIKE ? ", t).fetchall()
                if not res and (not miRNA_name[-1] in ['P', 'p']):
                    miRNA_mis[miRNA_name] = miRNA_name_dict[miRNA_name]
        return experiments, support_types, pmid_link, miRNA_mis

    def find_miRNA_count_gene(self, gene_names, of_those=None, limit=30):
        """
        For a given gene list, find the miRNAs regulate one of the genes, and the
        frequency of the miRNAs
        What(which of those) mirs most frequently or commonly regulate a list of genes
        """
        mirna_count = dict()
        temp = []
        mir_targets = defaultdict(list)
        if self.tfdb is not None:
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                         "WHERE target = ? ", t).fetchall()
                if res1:
                    temp.extend([r[0].upper() for r in res1])
                    for mir in [r[0].upper() for r in res1]:
                        mir_targets[mir].append(gene_name)
            if temp:
                mcount = Counter(temp)
                sorted_c = sorted(mcount.items(), key=operator.itemgetter(1), reverse=True)
                num = 0
                if of_those:
                    for mir,count in sorted_c:
                        if num < limit:
                            if mir in of_those and count > 1:
                                mirna_count[mir] = count
                                num += 1
                        else:
                            break
                else:
                    for mir,count in sorted_c:
                        if num < limit:
                            if count > 1:
                                mirna_count[mir] = count
                                num += 1
                        else:
                            break
        return mirna_count,mir_targets
                         
    def find_gene_count_miRNA(self, miRNA_name_dict, of_those=None, limit=30):
        """
        For a given miRNA list, find the genes regulated by one of the miRNAs, and the
        frequency of the genes
        What genes are most frequently (or commonly) regulated by a list of microRNAs?
        
        parameters
        ----------------
        miRNA_name_dict:dict
        of_those: list of gene names or None
        limit: int
        """
        gene_count = dict()
        temp = []
        target_mirna = defaultdict(list)
        miRNA_mis = {}
        miRNA_names = list(miRNA_name_dict.keys())
        if self.tfdb is not None:
            for mir in miRNA_names:
                t = (mir,)
                res1 = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                         "WHERE mirna LIKE ? ", t).fetchall()
                if res1:
                    temp.extend([r[0] for r in res1])
                    for target in [r[0] for r in res1]:
                        target_mirna[target].append(mir)
                else:
                     miRNA_mis[mir] = miRNA_name_dict[mir]
            if temp:
                c = Counter(temp)
                sorted_c = sorted(c.items(), key=operator.itemgetter(1), reverse=True)
                num = 0
                if of_those:
                    for gene,count in sorted_c:
                        if num < limit:
                            if gene in of_those and count > 1:
                                gene_count[gene] = count
                                num += 1
                        else:
                            break
                else:
                    for gene,count in sorted_c:
                        if num < limit:
                            if count > 1:
                                gene_count[gene] = count
                                num += 1
                        else:
                            break
        return gene_count,target_mirna,miRNA_mis

    def find_pathway_db_keyword(self, db_source, pathway_names):
        """
        Return a list of pathways which come from given db and whose name containing the given keywords
        Which reactome pathways involve immune signaling?
        """
        pathwayName = dict()
        dblink = dict()
        if self.tfdb is not None:
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr, db_source)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                        "WHERE pathwayName LIKE ? AND source LIKE ?", t).fetchall()
                if res:
                    for r in res:
                        pathwayName[r[0]] = r[1]
                        dblink[r[0]] = r[2]
            if not dblink:
                raise PathwayNotFoundException
        return pathwayName, dblink

    def gets_similar_miRNAs(self, miRNA_names):
        """
        return a list of miRNAs beginning with the given miRNA which 
        is not in the database, for user clarification purpose
        miRNA_names is a list
        """
        clari_miRNA = dict()
        if self.tfdb is not None:
            for miRNA_name in miRNA_names:
                t = (miRNA_name,)
                res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                         "WHERE mirna LIKE ? ", t).fetchall()
                if not res:
                    regstr = miRNA_name + '%'
                    t = (regstr,)
                    res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                            "WHERE mirna LIKE ? ", t).fetchall()
                    if res:
                        clari_miRNA[miRNA_name] = [r[0] for r in res]
        return clari_miRNA
        
    def get_similar_miRNAs(self, miRNA_name):
        """
        For each miRNA, return a list of miRNAs beginning with the given miRNA which 
        is not in the database, for user clarification purpose
        miRNA_name: str, miRNA name
        """
        clari_miRNA = []
        if self.tfdb is not None:
            regstr = miRNA_name + '-%'
            t = (regstr,)
            res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
            if res:
                clari_miRNA = [r[0] for r in res]
            else:
                regstr = miRNA_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
                if res:
                    clari_miRNA = [r[0] for r in res]
        return clari_miRNA

    def find_tissue_gene(self, gene_name):
        """
        For a given gene, return a list of tissue names it's expressed in
        What tissues is STAT3 expressed in?
        """
        tissue_names = []
        if self.tfdb is not None:
            t = (gene_name, EXP_THR)
            res = self.tfdb.execute("SELECT DISTINCT tissue FROM geneTissue "
                                    "WHERE genesymbol = ? AND enrichment > ? ", t).fetchall()
            if res:
                tissue_names = [r[0] for r in res]
            else:
                raise TissueNotFoundException
        return tissue_names
    
    def find_target_kinase(self, kinase_names):
        """
        For given kinases, return the genes regulated by them
        """
        gene_names = set()
        if self.tfdb is not None:
            t = (kinase_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT target FROM kinaseReg "
                                    "WHERE kinase = ? ", t).fetchall()
            if res:
                gene_names = set([r[0] for r in res])
            else:
                raise TargetNotFoundException 
            if len(kinase_names)>1:
                for i in range(1,len(kinase_names)):
                    t = (kinase_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT target FROM kinaseReg "
                                            "WHERE kinase = ? ", t).fetchall()
                    if res:
                        gene_names = set(gene_names) & set([r[0] for r in res])
                    else:
                        raise TargetNotFoundException
        if not gene_names:
            raise TargetNotFoundException
        
        return gene_names
        
    def find_target_kinase_keyword(self, kinase_names, keyword):
        """
        For given kinases, return the genes regulated by them
        """
        gene_names = set()
        if self.tfdb is not None:
            t = (kinase_names[0], keyword)
            res = self.tfdb.execute("SELECT DISTINCT target FROM kinaseReg "
                                    "WHERE kinase = ? AND direction LIKE ? ", t).fetchall()
            if res:
                gene_names = set([r[0] for r in res])
            else:
                raise TargetNotFoundException 
            if len(kinase_names)>1:
                for i in range(1,len(kinase_names)):
                    t = (kinase_names[i], keyword)
                    res = self.tfdb.execute("SELECT DISTINCT target FROM kinaseReg "
                                            "WHERE kinase = ? AND direction LIKE ? ", t).fetchall()
                    if res:
                        gene_names = set(gene_names) & set([r[0] for r in res])
                    else:
                        raise TargetNotFoundException
        if not gene_names:
            raise TargetNotFoundException
        return gene_names
    
    def find_kinase_target(self, target_names):
        """
        For given genes, return the kinases that regulate them
        """
        kinase_names = []
        if self.tfdb is not None:
            t = (target_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT kinase FROM kinaseReg "
                                    "WHERE target = ? ", t).fetchall()
            if res:
                kinase_names = [r[0] for r in res]
            else:
                raise KinaseNotFoundException 
            if len(target_names)>1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT kinase FROM kinaseReg "
                                            "WHERE target = ? ", t).fetchall()
                    if res:
                        kinase_names = list(set(kinase_names) & set([r[0] for r in res]))
                    else:
                        raise KinaseNotFoundException
        return kinase_names
        
    def find_kinase_target_keyword(self, target_names, keyword_name):
        """
        For given genes and keyword, return kinases that regulate them with keyword (up or down) mode
        """
        kinase_names = []
        if self.tfdb is not None:
            t = (target_names[0], keyword_name)
            res = self.tfdb.execute("SELECT DISTINCT kinase FROM kinaseReg "
                                    "WHERE target = ? AND direction LIKE ? ", t).fetchall()
            if res:
                kinase_names = [r[0] for r in res]
            else:
                raise KinaseNotFoundException 
            if len(target_names)>1:
                for i in range(1,len(target_names)):
                    t = (target_names[i], keyword_name)
                    res = self.tfdb.execute("SELECT DISTINCT kinase FROM kinaseReg "
                                            "WHERE target = ? AND direction LIKE ? ", t).fetchall()
                    if res:
                        kinase_names = list(set(kinase_names) & set([r[0] for r in res]))
                    else:
                        raise KinaseNotFoundException
        return kinase_names
        
    def find_gene_tissue(self, tissue_name):
        """
        For a given tissue, return genes expressed in this tissue
        """
        gene_names = []
        if self.tfdb is not None:
            regstr = '%' + tissue_name + '%'
            t = (regstr, EXP_THR)
            res = self.tfdb.execute("SELECT DISTINCT genesymbol FROM geneTissue "
                                    "WHERE tissue LIKE ? AND enrichment > ? ", t).fetchall()
            if res:
                gene_names = [r[0] for r in res]
            else:
                raise TissueNotFoundException
        return gene_names
        
    def find_gene_tissue_exclusive(self, tissue_name):
        """
        For a given tissue, return genes exclusively expressed in this tissue
        """
        gene_names = []
        if not len(self.tissue_gene_exclusive):
            self.tissue_gene_exclusive = self.map_exclusive_tissue_gene()
        if tissue_name in self.tissue_gene_exclusive:
            gene_names = list(self.tissue_gene_exclusive[tissue_name])
        return gene_names
        
    def is_tissue_gene(self, tissue_name, gene_name):
        """
        For a given gene and a tissue, return if this gene is expressed in this tissue
        """
        if self.tfdb is not None:
            t = (tissue_name, gene_name, EXP_THR)
            res = self.tfdb.execute("SELECT DISTINCT genesymbol FROM geneTissue "
                                    "WHERE tissue LIKE ? AND genesymbol = ? "
                                    "AND enrichment > ? ", t).fetchall()
            if res:
                return True
            else:
                return False
    
    def is_tissue_gene_exclusive(self, tissue_name, gene_name):
        """
        For a given gene and a tissue, return if this gene is exclusively expressed in this tissue
        """
        if self.tfdb is not None:
            t = (gene_name, EXP_THR)
            res = self.tfdb.execute("SELECT DISTINCT genesymbol, tissue FROM geneTissue "
                                    "WHERE genesymbol = ? AND enrichment > ? ", t).fetchall()
            if res:
                tissues = set([r[1] for r in res])
                if len(tissues) == 1 and tissue_name in tissues:
                    return True
                else:
                    return False
            else:
                return False
                
    def map_exclusive_tissue_gene(self):
        """
        return gene expressions exclusively in each tissue
        """
        gene_exp_all = defaultdict(set)
        gene_exp_exclusive = defaultdict(set)
        if self.tfdb is not None:
            t = (EXP_THR,)
            res = self.tfdb.execute("SELECT DISTINCT tissue FROM geneTissue "
                                    "WHERE enrichment > ? ", t).fetchall()
            tissues = [r[0] for r in res]
            for tiss in tissues:
                t = (tiss, EXP_THR)
                res = self.tfdb.execute("SELECT DISTINCT genesymbol FROM geneTissue "
                                    "WHERE tissue LIKE ? AND enrichment > ? ", t).fetchall()
                gene_exp_all[tiss] = set([r[0] for r in res])
            #get exclusively expressed genes
            for i in range(len(tissues)):
                temp = gene_exp_all[tissues[i]]
                for j in range(len(tissues)):
                    if j != i:
                        temp = temp - gene_exp_all[tissues[j]]
                #print('tissue ' + tissues[i] + ' have %d exclusively expressively expressed genes' % len(temp))
                if len(temp):
                    gene_exp_exclusive[tissues[i]] = temp
        return gene_exp_exclusive
    
    def find_evidence_dbname(self, tf_name, target_name):
        """
        Return the dbnames supporting the regulation
        """
        if self.tfdb is not None:
            t = (tf_name, target_name)
            res = self.tfdb.execute("SELECT DISTINCT dbnames FROM CombinedDB "
                                    "WHERE TF = ? AND Target = ? ", t).fetchall()
            if res:
                dbs = [r[0] for r in res]
                db_names = set(','.join(dbs).split(','))
            else:
                db_names = set()
        return db_names
    
    def is_gene_disease(self, gene, disease, keyword):
        """
        Return true if the gene is regulated in the disease in the keyword direction, otherwise false
        """
        if self.ldd is not None:
            reg1 = '%' + disease + '%'
            if keyword.lower() == 'regulate':
                t = (gene, reg1)
                res = self.ldd.execute("SELECT * FROM diseaseGene WHERE gene = ? AND diseaseId IN "
                      "(SELECT Id FROM diseaseName WHERE disease LIKE ?)", t).fetchall()
            else:
                t = (gene, reg1, keyword)
                res = self.ldd.execute("SELECT * FROM diseaseGene WHERE gene = ? AND diseaseId IN "
                     "(SELECT Id FROM diseaseName WHERE disease LIKE ? AND direction LIKE ?)", t).fetchall()
            if res:
                return True
            else:
                return False
    
    def find_gene_disease(self, disease, keyword, of_those):
        """
        Return genes perturbated in the disease
        For 'regulate', combine the results from 'increase' and 'decrease'
        """
        dname = dict()
        genes = defaultdict(set)
        if self.ldd is not None:
            reg1 = '%' + disease + '%'
            if keyword.lower() == 'regulate':
                t = (reg1, )
                res = self.ldd.execute("SELECT Id, disease FROM diseaseName "
                                    "WHERE disease LIKE ? ", t).fetchall()
            else:
                t = (reg1, keyword)
                res = self.ldd.execute("SELECT Id, disease FROM diseaseName "
                                    "WHERE disease LIKE ? AND direction LIKE ? ", t).fetchall()
            if res:
                for r in res:
                    dname[r[0]] = r[1]
            if dname:
                for id in dname:
                    t = (id,)
                    res = self.ldd.execute("SELECT DISTINCT gene FROM diseaseGene "
                                    "WHERE diseaseId = ? ", t).fetchall()
                    if of_those:
                        temp = set(of_those).intersection(set([r[0] for r in res]))
                        if temp:
                            genes[dname[id]] = genes[dname[id]].union(temp)
                    else:
                        genes[dname[id]] = genes[dname[id]].union(set([r[0] for r in res]))
        return genes
        
    def find_gene_ligand(self, ligand_name, keyword, of_those, limit=10):
        """
        Return genes perturbated by the ligand
        """
        lname = dict()
        genes = dict()
        if self.ldd is not None:
            ligands = [ligand_name, ligand_name.replace(' ', '-')]
            for ligand in ligands:
                reg1 = '%' + ligand + '%'
                t = (reg1, keyword)
                res = self.ldd.execute("SELECT Id, ligand FROM ligandName "
                                    "WHERE ligand LIKE ? AND direction LIKE ? ", t).fetchall()
                if res:
                    for r in res:
                        lname[r[0]] = r[1]
            if lname:
                num = 1
                for id in lname:
                    t = (id,)
                    res = self.ldd.execute("SELECT DISTINCT gene FROM ligandGene "
                                    "WHERE ligandId = ? ORDER BY gene", t).fetchall()
                    if of_those:
                        temp = set(of_those).intersection(set([r[0] for r in res]))
                        if temp:
                            genes[id] = temp
                            num += 1
                    else:
                        genes[id] = [r[0] for r in res]
                        num += 1
                    if num > limit:
                        break
        return lname, genes
        
    def find_gene_drug(self, drug, keyword, of_those, limit=10):
        """
        Return genes perturbated by the drug
        """
        dname = dict()
        genes = dict()
        if self.ldd is not None:
            reg1 = '%' + drug + '%'
            t = (reg1, keyword)
            res = self.ldd.execute("SELECT Id, drug FROM drugName "
                                    "WHERE drug LIKE ? AND direction LIKE ? ", t).fetchall()
            if res:
                for r in res:
                    dname[r[0]] = r[1]
            if dname:
                num = 1
                for id in dname:
                    t = (id,)
                    res = self.ldd.execute("SELECT DISTINCT gene FROM drugGene "
                                    "WHERE drugId = ? ORDER BY gene", t).fetchall()
                    if of_those:
                        temp = set(of_those).intersection(set([r[0] for r in res]))
                        if temp:
                            genes[id] = temp
                            num += 1
                    else:
                        genes[id] = [r[0] for r in res]
                        num += 1
                    num += 1
                    if num > limit:
                        break
        return dname, genes
    
    def get_onto_set(self, go_name):
        """
        Return the genes which are in the category of go_name
        """
        go_genes = []
        go_name = go_name.lower()
        if go_name in ['tf', 'transcription factor']:
            if not self.trans_factor:
                if self.tfdb is not None:
                    res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                    self.trans_factor = set([r[0] for r in res])
            return self.trans_factor
        else:
            try:
                goids = go_map[go_name]
            except KeyError:
                raise GONotFoundException
            if self.tfdb is not None:
                for go in goids:
                    t = (go, )
                    res = self.tfdb.execute("SELECT geneSymbol FROM go2Genes WHERE termId = ? ", t).fetchall()
                    go_genes.extend([r[0] for r in res])
                go_genes = set(go_genes)
        return go_genes
        
    def get_hgnc_mapping(self):
        return hgnc_symbol_id
    
    def get_hgnc_symbols(self):
        return hgnc_genes_set
    
    def get_tf_set(self):
        tf_set = set()
        if not self.trans_factor:
            if self.tfdb is not None:
                res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                self.trans_factor = set([r[0] for r in res])
                tf_set = self.trans_factor
        else:
            tf_set = self.trans_factor
        return tf_set
        
    def tf_set(self):
        tfs = set()
        if self.tfdb is not None:
            res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
            tfs = set([r[0] for r in res])
        return tfs
        
    def mirna_set(self):
        mirnas = set()
        if self.tfdb is not None:
            res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo").fetchall()
            mirnas = set([r[0] for r in res])
        return mirnas
    
    def get_pathway_genes(self, db_source):
        """
        Return pathway-gene dict. 
        """
        p_genes = defaultdict(dict)
        fn = os.path.join(_enrich_dir, db_source.lower() + '.pickle')
        if os.path.isfile(fn):
            with open(fn, 'rb') as pickle_in:
                p_genes = pickle.load(pickle_in)
        else:
            if self.tfdb is not None:
                pathw = defaultdict(dict)
                reg_str = db_source + '%'
                t = (reg_str,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                        "WHERE source LIKE ?", t).fetchall()
                for r in res:
                    pathw[r[0]]['name'] = r[1]
                    pathw[r[0]]['dblink'] = r[2]
                
                for id in pathw.keys():
                    t = (id,)
                    res = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    p_genes[pathw[id]['name']]['gene'] = [r[0] for r in res]
                    p_genes[pathw[id]['name']]['dblink'] = pathw[id]['dblink']
                #save to file
                with open(fn, 'wb') as pickle_out:
                    pickle.dump(p_genes, pickle_out)
        return p_genes
    
    
        
    
    #---------------------------------------------#
    #--------methods for querying indra database--#
    #---------------------------------------------#
    def find_statement_indraDB(self, subj=None, obj=None, stmt_types=None):
        """
        subj: str, HGNC symbol
        obj: str, HGNC symbol
        stmt_type: list[str], list of statement type
        """
        if not CAN_CHECK_STATEMENTS:
            return [],False
        statements = []
        try:
            for stype in stmt_types:
                stmts = get_statements(subject=subj, object=obj, stmt_type=stype, simple_response=True, ev_limit=5, max_stmts=100, persist=False)
                if len(stmts):
                    stmts_filtered = filter_evidence_source(stmts, ['reach'], policy='one')
                    if len(stmts_filtered):
                        statements.extend(stmts_filtered)
        except Exception as e:
            #print(e)
            return statements,False
        return statements,True
            
    def find_regulator_indra(self, stmts):
        """
        stmts: indra statements
        return the list of TFs from the subject of the stmts, as well as other subjects
        """
        subjects = set()
        for stmt in stmts:
            subj = stmt.subj
            if subj is not None:
                subjects.add(subj.name)
        #get all the tfs in the db
        if not self.trans_factor:
            if self.tfdb is not None:
                res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                self.trans_factor = set([r[0] for r in res])
       
        tfs = subjects.intersection(self.trans_factor)
        nontfs = subjects - tfs
        #mirnas = nontfs.intersection(mirna_indra_set)
        #others = nontfs - mirnas
        genes = nontfs.intersection(hgnc_genes_set)
        #other = others - genes
        return tfs, genes#, mirnas, other
        
    def find_regulator_indra_gene(self, stmts, of_those=None):
        """
        stmts: indra statements
        of_those: set or None
        return the list of genes from the object of the stmts
        """
        subjects = set()
        for stmt in stmts:
            subj = stmt.subj
            if subj is not None:
                subjects.add(subj.name)
        #genes = subjects.intersection(hgnc_genes_set)
        if of_those:
            genes = subjects.intersection(of_those)
        else:
            genes = subjects
        return genes
        
    def find_regulators_indra(self, stmts_d):
        """
        stmts_d: dict
        """
        targets = list(stmts_d.keys())
        stmt_list = []
        for t in targets:
            stmt_list.extend(stmts_d[t])
        genes = self.find_regulator_indra_gene(stmts_d[targets[0]])
        if len(targets) > 1:
            for i in range(1,len(targets)):
                temp = self.find_regulator_indra_gene(stmts_d[targets[i]])
                genes = genes.intersection(temp)
        #filter statements
        stmt_f = []
        for stmt in stmt_list:
            subj = stmt.subj
            if subj:
                if subj.name in genes:
                    stmt_f.append(stmt)
        
        ##get all the tfs in the db
        if not self.trans_factor:
            if self.tfdb is not None:
                res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                self.trans_factor = set([r[0] for r in res])
       
        tfs = genes.intersection(self.trans_factor)
        nontfs = genes - tfs
        mirnas = nontfs.intersection(self.mirna)
        others = nontfs - mirnas
        gene = others.intersection(hgnc_genes_set)
        #other = others - gene
        return tfs, gene, mirnas, stmt_f
        
    def find_regulator_indra_targets(self, stmts_d, of_those=None, target_type=None):
        """
        stmts_d: dict
        return the list of genes, as well as the filtered statements
        """
        targets = list(stmts_d.keys())
        stmt_list = []
        for t in targets:
            stmt_list.extend(stmts_d[t])
        genes = self.find_regulator_indra_gene(stmts_d[targets[0]], of_those=of_those)
        
        if target_type:
            try:
                target_type_set = self.get_onto_set(target_type)
            except Exception as e:
                target_type_set = []
            if target_type_set:
                genes = genes.intersection(target_type_set)
        
        if len(targets) > 1:
            for i in range(1,len(targets)):
                temp = self.find_regulator_indra_gene(stmts_d[targets[i]])
                genes = genes.intersection(temp)
                
        #filter statements
        stmt_f = []
        for stmt in stmt_list:
            subj = stmt.subj
            if subj:
                if subj.name in genes:
                    stmt_f.append(stmt)
        return genes, stmt_f
        
    def find_target_indra(self, stmts):
        """
        stmts: indra statements
        return the list of genes from the object of the stmts
        """
        objs = set()
        for stmt in stmts:
            obj = stmt.obj
            if obj is not None:
                objs.add(obj.name)
        genes = objs.intersection(hgnc_genes_set)
        return genes
        
    def find_target_indra_regulators(self, stmts_d, of_those=None, target_type=None):
        """
        stmts_d: dict
        return the list of genes, as well as the filtered statements
        """
        regulators = list(stmts_d.keys())
        stmt_list = []
        for r in regulators:
            stmt_list.extend(stmts_d[r])
        genes = self.find_target_indra(stmts_d[regulators[0]])
        
        if of_those:
            genes = genes.intersection(of_those)
            
        if target_type and isinstance(target_type, set):
            genes = genes.intersection(target_type)
        else:
            try:
                target_type_set = self.get_onto_set(target_type)
            except Exception as e:
                target_type_set = []
            if target_type_set:
                genes = genes.intersection(target_type_set)
            
        if len(regulators) > 1:
            for i in range(1, len(regulators)):
                temp = self.find_target_indra(stmts_d[regulators[i]])
                genes = genes.intersection(temp)
        #filter statements
        stmt_f = []
        for stmt in stmt_list:
            obj = stmt.obj
            if obj:
                if obj.name in genes:
                    stmt_f.append(stmt)
        return genes, stmt_f
        
    def find_tf_indra(self, stmts):
        """
        stmts: indra statements
        return the list of tfs from the subject of the stmts
        """
        subjs = set()
        for stmt in stmts:
            subj = stmt.subj
            if subj is not None:
                subjs.add(subj.name)
                
        #get all the tfs in the db
        if not self.trans_factor:
            if self.tfdb is not None:
                res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                self.trans_factor = set([r[0] for r in res])
            
        tfs = subjs.intersection(self.trans_factor)
        return tfs
        
    def find_tfs_indra(self, stmts_d):
        """
        stmts_d: dict
        return the list of tfs, as well as the filtered statements
        """
        targets = list(stmts_d.keys())
        stmt_list = []
        for t in targets:
            stmt_list.extend(stmts_d[t])
        tfs = self.find_tf_indra(stmts_d[targets[0]])
        if len(targets) > 1:
            for i in range(1, len(targets)):
                temp = self.find_tf_indra(stmts_d[targets[i]])
                tfs = tfs.intersection(temp)
        #filter statements
        stmt_f = []
        for stmt in stmt_list:
            subj = stmt.subj
            if subj:
                if subj.name in tfs:
                    stmt_f.append(stmt)
        return tfs, stmt_f
        
    def find_evidence_indra(self, stmts):
        """
        stmts: indra statements
        return a set of tuple(source_api, pmid, text)
        """
        evidences = set()
        for st in stmts:
            evis = st.evidence
            source_api = ''
            pmid = ''
            text = ''
            for ev in evis:
                source_api = ev.source_api
                pmid = ev.pmid
                text = ev.text
                evidences.add((source_api, pmid, text))
        return evidences
    
    def find_member(self, agent):
        """
        Find members for a collection
        
        Parameters
        -------------
        agent: Agent

        Returns
        ------------
        magents: list of Agent
        """
        magents = _get_members(agent)
        return magents
        
    def find_members(self, agents):
        """
        Find members for a collection
        
        Parameters
        -------------
        agents: list of Agent
        Value is indra.statements.Agent
        
        Returns
        ------------
        members: dict
        Its key is family name
        Value is list[indra.statements.Agent]
        """
        members = dict()
        for agent in agents:
            magents = _get_members(agent)
            if magents:
                members[agent.name] = magents
        return members
    
    #-------------------------------------------------------------------------------------
    @staticmethod
    def load_db():
        logger.debug('Using resource folder: %s' % _resource_dir)
        #Load TF_target database
        tf_db_file = _resource_dir + 'TF_target_20191224.db'
        if not os.path.exists(tf_db_file):
            logger.info('Downloading TF_target db file...')
            url = 'https://www.dropbox.com/s/gjoal1xe420je6o/TF_target_20191224.db?dl=1'
            download_file_dropbox(url, tf_db_file)
            
        if os.path.isfile(tf_db_file):
            tfdb = sqlite3.connect(tf_db_file, check_same_thread=False)
            logger.info('TFTA loaded TF-target database')
        else:
            logger.error('TFTA could not load TF-target database.')
            tfdb = None;
        return tfdb
        
    @staticmethod
    def load_ldd_db():
        logger.info('Loading ldd.db...')
        ldd_file = os.path.join(_resource_dir, 'ldd.db')
        if not os.path.exists(ldd_file):
            logger.info('Downloading ldd.db file...')
            url = 'https://www.dropbox.com/s/zihor8mhpozle1e/ldd.db?dl=1'
            download_file_dropbox(url, ldd_file)
            
        if os.path.isfile(ldd_file):
            ldd = sqlite3.connect(ldd_file, check_same_thread=False)
            logger.info('TFTA loaded ldd database')
        else:
            logger.error('TFTA could not load ldd database.')
            ldd = None
        return ldd

def _get_members(agent):
    return expand_agent(agent, bio_ontology, ns_filter=['HGNC'])


def _check_overlap(genes, fmembers):
    """
    If genes overlap with any member of each family in fmembers, return True, otherwise False
        
    parameter
    --------------
    genes: set
    fmembers: defaultdict(list)
    """
    ind = []
    for f in fmembers:
        fg = set()
        for m in fmembers[f]:
            fg.add(m.name)
        ind.append(len(genes.intersection(fg)))
    return all(ind)
        
#test functions
#if __name__ == "__main__":
    #a = TFTA()
    
