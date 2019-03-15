#TFTA stands for TF-target agent whose task is to search for targets known to a TF and vice versa

import re
import os
import operator
import logging
import sqlite3
import numpy as np
from collections import defaultdict
import math
from indra import has_config

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

mirna_indra_set = _get_mirna_indra()

#hgnc official symbol to id mapping
def _get_hgnc_genes():
    lines = open(_resource_dir + 'hgnc_symbol_id_20190225.txt', 'rt').readlines()
    hgnc_genes = dict()
    for line in lines:
        s = line.strip().split('\t')
        hgnc_genes[s[1]] = s[0]
    return hgnc_genes

hgnc_symbol_id = _get_hgnc_genes()
hgnc_genes_set = set(hgnc_symbol_id.keys())

#gene expression threshold
EXP_THR = 1.5

class TFTA:
    #transcription factor list, will set when it's first used
    trans_factor = set()
    
    def __init__(self):
        logger.debug('Using resource folder: %s' % _resource_dir)
        #Load TF_target database
        tf_db_file = _resource_dir + 'TF_target_v7_1.db'
        if os.path.isfile(tf_db_file):
            self.tfdb = sqlite3.connect(tf_db_file, check_same_thread=False)
            logger.info('Loaded TF-target database')
            self.tfdb.row_factory = sqlite3.Row
        else:
            logger.error('TFTA could not load TF-target database.')
            self.tfdb = None;
        #for exclusive query
        self.tissue_gene_exclusive = defaultdict(set)
            
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
        if go_name == 'transcription factor':
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
        if go_name == 'transcription factor':
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
                    go_genes += [r[0] for r in res]
                go_genes = list(set(go_genes) & set(gene_names))
        if len(go_genes):
            go_genes.sort()
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
                tf_names = [r[0] for r in res]
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
                        tf_names = list(set(tf_names) & set([r[0] for r in res]))
                        for r in res:
                            dbname[(r[0],target_names[i])] = r[1]
                    else:
                        raise TargetNotFoundException
            if len(tf_names):
                tf_names.sort()
        else:
            tf_names = []
        return tf_names,dbname

    def find_tfs_count(self,target_names):
        """
        Return TFs regulating the given target list as well as the frequency (count)
        of each TF
        """
        #query
        tf_counts = dict()
        tf_names = []
        targets = defaultdict(list)
        if self.tfdb is not None:
            target_names = list(set(target_names))
            for target_name in target_names:
                t = (target_name,)
                res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
                                        "WHERE Target = ? ", t).fetchall()
                if res:
                    tfs = [r[0] for r in res]
                    tf_names = tf_names + tfs
                    for tf in tfs:
                        targets[tf].append(target_name)
            #count for each TF
            if len(tf_names):
                unique_tfs = list(set(tf_names))
                for i in range(len(unique_tfs)):
                    ct = tf_names.count(unique_tfs[i])
                    if len(target_names) > 1:
                        if ct > 1:
                            tf_counts[unique_tfs[i]] = ct
                    else:
                        tf_counts[unique_tfs[i]] = ct
            if len(tf_counts):
                tf_counts = sorted(tf_counts.items(), key=operator.itemgetter(1))
                tf_counts.reverse()
                #tf_counts, unique_tfs = list(zip(*sorted(zip(tf_counts, unique_tfs),reverse=True)))
            else:
                raise TFNotFoundException	
        return tf_counts
        
    def find_common_tfs(self,target_names):
        """
        Return TFs regulating the given target list as well as the regulated targets by each TF
        """
        #query
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
                tf_names = [r[0] for r in res]
            else:
                #check if the target is in the database
                t = (target_names[0],)
                res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                        "WHERE Target = ? ", t).fetchall()
                if res:
                    return tf_names
                else:
                    raise TargetNotFoundException	
            if len(target_names) > 1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],regstr)
                    res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                            "WHERE Target = ? AND Tissue LIKE ? ", t).fetchall()
                    if res:
                        tf_names = list(set(tf_names) & set([r[0] for r in res]))
                    else:
                        #check if the target is in the database
                        t = (target_names[i],)
                        res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                                "WHERE Target = ? ", t).fetchall()
                        if res:
                            return tf_names
                        else:
                            raise TargetNotFoundException
            if len(tf_names):
                tf_names.sort()
        return tf_names

    def find_pathways_from_name(self, pathway_name):
        """
        return pathway information related to pathway_name
        """
        #query
        if self.tfdb is not None:
            regstr = '%' + pathway_name + ' %'
            t = (regstr,)
            #get pathwayId
            res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ORDER BY pathwayName", t).fetchall()
            if res:
                pathwayId = [r[0] for r in res]
                pathwayName = [r[1] for r in res]
                externalId = [r[2] for r in res]
                source = [r[3] for r in res]
                dblink = [r[4] for r in res]
            else:
                raise PathwayNotFoundException
        else:
            pathwayId = []
            pathwayName = []
            externalId = []
            source = []
            dblink = []
        return pathwayId,pathwayName,externalId,source,dblink

    def find_pathways_from_dbsource_geneName1(self, dbsource,gene_name):
        """
        return pathway information for given dbsource and gene_name
        """
        pathwayId = []
        pathwayName = []
        externalId = []
        psource = []
        dblink = []
        if self.tfdb is not None:
            #regstr='%'+pathway_name+'%'
            t = (gene_name[0],)
            res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                    "WHERE genesymbol = ? ", t).fetchall()
            if res1:
                pids = [r[0] for r in res1]
            else:
                raise PathwayNotFoundException
            
            for pid in pids:  
                t = (pid, dbsource)
                res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                    "WHERE Id = ? AND source LIKE ? ", t).fetchall()
            
                if res:
                    pathwayId = pathwayId + [r[0] for r in res]
                    pathwayName = pathwayName + [r[1] for r in res]
                    externalId = externalId + [r[2] for r in res]
                    psource = psource + [r[3] for r in res]
                    dblink = dblink + [r[4] for r in res]
            #sort
            if len(pathwayName):
                pathwayName,pathwayId,externalId,psource,dblink = \
                     list(zip(*sorted(zip(pathwayName,pathwayId,externalId,psource,dblink))))
        return pathwayId,pathwayName,externalId,psource,dblink

    def find_pathways_from_dbsource_geneName(self, dbsource,gene_names):
        """
        return pathway information for given dbsource and gene_names
        """
        #query
        pathwayId = []
        pathwayName = []
        externalId = []
        psource = []
        dblink = []
        if self.tfdb is not None:
            pathlist = []
            for gene_name in gene_names:
                t = (gene_name,dbsource)
                print(t)
                res = self.tfdb.execute("SELECT Id FROM pathwayInfo "
                                   "WHERE Id in (SELECT DISTINCT pathwayID FROM pathway2Genes "
                                   "WHERE genesymbol = ?) AND source LIKE ? ", t).fetchall()
                if res:
                    pathlist = pathlist + [r[0] for r in res]
                else:
                    raise PathwayNotFoundException
            #interaction
            pathIDs = []
            for pth in set(pathlist):
                if pathlist.count(pth) == len(gene_names):
                    pathIDs.append(pth)
            if len(pathIDs):
                for pth in pathIDs:
                    t = (pth,)
                    res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                   "WHERE Id = ? ", t).fetchall()
                    pathwayId = pathwayId + [r[0] for r in res]
                    pathwayName = pathwayName + [r[1] for r in res]
                    externalId = externalId + [r[2] for r in res]
                    psource = psource + [r[3] for r in res]
                    dblink = dblink + [r[4] for r in res]
            else:
                raise PathwayNotFoundException    
            #sort
            pathwayName,pathwayId,externalId,psource,dblink = \
                    list(zip(*sorted(zip(pathwayName,pathwayId,externalId,psource,dblink))))
        return pathwayId,pathwayName,externalId,psource,dblink

    def find_genes_from_pathwayName(self, pathway_names):
        """
        return genes related to pathway_name
        """
        pathwayId = []
        pathwayName = dict()
        genelist = dict()
        pw_link = dict()
        #query
        if self.tfdb is not None:
            pn = []
            pids = []
            plink = []
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                #get pathwayId
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    pids = pids + [r[0] for r in res]
                    pn = pn + [r[1] for r in res]
                    plink = plink + [r[2] for r in res]    
            if len(pids):
                pathwayId = list(set(pids))
                for i in range(len(pids)):
                    pathwayName[pids[i]] = pn[i]
                    pw_link[pids[i]] = plink[i]
                for pthID in pathwayId:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                           "WHERE pathwayID = ? ORDER BY genesymbol", t).fetchall()
                    genes = [r[0] for r in res1]
                    genelist[pthID] = genes
            else:
                raise PathwayNotFoundException
        return pathwayId,pathwayName,genelist,pw_link

    def find_tf_pathway(self,pathway_names):
        """
        Return TFs within the given pathway
        """
        pathwayName = dict()
        dblink = dict()
        tflist = dict()
        if self.tfdb is not None:
            pids = []
            pn = []
            dl = []
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                        "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    pids += [r[0] for r in res]
                    pn += [r[1] for r in res]
                    dl += [r[2] for r in res]
            if len(pids):
                pathwayId = list(set(pids))
                for i in range(len(pids)):
                    pathwayName[pids[i]] = pn[i]
                    dblink[pids[i]] = dl[i]
                for pthID in pathwayId:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? AND isTF = 1 ORDER BY genesymbol", t).fetchall()
                    if res1:
                        tfs = [r[0] for r in res1]
                        tflist[pthID] = tfs
                    else:
                        del pathwayName[pthID]
                        del dblink[pthID]
                        #pathwayId.remove(pthID)
            else:
                raise PathwayNotFoundException
        return pathwayName,tflist,dblink
        
    def find_kinase_pathway(self,pathway_names):
        """
        Return kinases within the given pathway
        """
        pathwayId = []
        pathwayName = dict()
        dblink = dict()
        kinaselist = dict()
        if self.tfdb is not None:
            #get kinases
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
            if len(pathwayName):
                pathwayId = pathwayName.keys()
                for pthID in pathwayId:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? ", t).fetchall()
                    kins = list(kinases & set([r[0] for r in res1]))
                    if len(kins):
                        kins.sort()
                        kinaselist[pthID] = kins
            else:
                raise PathwayNotFoundException
            if not len(kinaselist):
                raise PathwayNotFoundException
        return pathwayName, kinaselist, dblink

    def find_pathways_from_genelist(self,gene_names):
        """
        return pathways having given genes
        """
        pathwayId=[]
        pathwayName=[]
        dblink=[]
        if self.tfdb is not None:
            pathlist=[]
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pathlist = pathlist+[r[0] for r in res1]
                else:
                    raise PathwayNotFoundException
             #intersection
            if len(gene_names)>1:
                pathIDs=[]
                for pth in set(pathlist):
                    if pathlist.count(pth) == len(gene_names):
                        pathIDs.append(pth)
            else:
                pathIDs = pathlist
            if len(pathIDs):
                for pth in pathIDs:
                    t = (pth,)
                    res = self.tfdb.execute("SELECT pathwayName,dblink FROM pathwayInfo "
                                            "WHERE Id = ? ", t).fetchall()
                    #pathwayId = pathwayId + [r[0] for r in res]
                    pathwayName = pathwayName + [r[0] for r in res]
                    #externalId = externalId + [r[2] for r in res]
                    #source = source + [r[3] for r in res]
                    dblink = dblink + [r[1] for r in res]
                #sort
                pathwayName,dblink = \
                    list(zip(*sorted(zip(pathwayName,dblink))))
            else:
                raise PathwayNotFoundException	
        return pathwayName,dblink

    def Is_pathway_gene(self, pathway_names, gene_names):
        """
        Return pathways which contain the given genes and whose name contain substring of pathway_name
        """
        pname = dict()
        plink = dict()
        upid = []
        if self.tfdb is not None:
            pids = []
            pn = []
            dl = []
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
                if res:
                    pids = pids + [r[0] for r in res]
                    pn = pn + [r[1] for r in res]
                    dl = dl + [r[2] for r in res]
                    
            if len(pids):
                for id,n,l in zip(pids,pn,dl):
                    pname[id] = n
                    plink[id] = l
            else:
                raise PathwayNotFoundException
                    
            upid = list(set(pids))
            fpid = []
            for pthID in upid:
                t = (pthID,)
                res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                        "WHERE pathwayID = ? ORDER BY genesymbol", t).fetchall()
                genes = [r[0] for r in res1]
                overlap_genes = list(set(gene_names) & set(genes))
                if len(overlap_genes) == len(gene_names):
                    fpid.append(pthID)
            if not len(fpid):
                raise PathwayNotFoundException
            return fpid, pname, plink

    def find_pathway_gene_keyword(self,gene_names, keyword):
        """
        return pathways having given genes and some information in pathway name
        """
        pathwayId = []
        pathwayName1 = []
        externalId = []
        source = []
        dblink = []
        if self.tfdb is not None:
            pathlist = []
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pathlist = pathlist + [r[0] for r in res1]
                else:
                    raise PathwayNotFoundException
            #intersection
            pathIDs = []
            for pth in set(pathlist):
                if pathlist.count(pth) == len(gene_names):
                    pathIDs.append(pth)
            if len(pathIDs)>0:
                regstr = '%' + keyword + '%'
                for pth in pathIDs:
                    t = (pth, regstr)
                    res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                            "WHERE Id = ? AND pathwayName LIKE ? ", t).fetchall()                    
                    if res:
                        pathwayId = pathwayId + [r[0] for r in res]
                        pathwayName1 = pathwayName1 + [r[1] for r in res]
                        externalId = externalId + [r[2] for r in res]
                        source = source + [r[3] for r in res]
                        dblink = dblink + [r[4] for r in res]
                #sort
                if len(pathwayId):
                    pathwayName1,pathwayId,externalId,source,dblink = \
                        list(zip(*sorted(zip(pathwayName1,pathwayId,externalId,source,dblink))))
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException    
        return pathwayName1, dblink

    def find_pathway_keyword(self, keyword):
        """
        return pathways containing the given keyword
        """
        if self.tfdb is not None:
            regstr = '%' + keyword + '%'
            t = (regstr,)
            res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ORDER BY pathwayName", t).fetchall()
            if res:
                pathwayName = [r[1] for r in res]
                dblink = [r[4] for r in res]
            else:
                raise PathwayNotFoundException		
        else:
            pathwayName = []
            dblink = []
        return pathwayName,dblink

    def find_tf_keyword(self, keyword_name):
        """
        return pathways containing the given keyword and related tfs 
        """
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
                tflist = dict()
                newpathwayName = []
                newpathwayId = []
                newdblink = []
                for pn,pthID,lk in zip(pathwayName,pathwayId,dblink):
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? AND isTF = 1 ORDER BY genesymbol", t).fetchall()
                    if res1:
                        tfs = [r[0] for r in res1]
                        tflist[pthID] = tfs
                        newpathwayId.append(pthID)
                        newpathwayName.append(pn)
                        newdblink.append(lk)
            else:
                raise PathwayNotFoundException
        else:
            newpathwayId = []
            newpathwayName = []
            tflist = dict()
            newdblink = []
        if not len(newdblink):
            raise PathwayNotFoundException
        return newpathwayId,newpathwayName,tflist,newdblink

    def find_pathway_count_genes(self, gene_names):
        """
        For a given gene list, find the pathways containing one of the genes,
        and the frequency of the pathways 
        """
        pathwayName = []
        externalId = []
        source = []
        dblink = []
        counts = []
        if self.tfdb is not None:
            pathlist = []
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pathlist = pathlist + [r[0] for r in res1]
            #pathway frequency
            uniq_path = list(set(pathlist))
            if len(uniq_path) and len(gene_names)>1:
                uniq_path2 = []
                for pth in uniq_path:
                    ct = pathlist.count(pth)
                    if ct > 1:
                        counts.append(ct)
                        uniq_path2.append(pth)
            else:
                uniq_path2 = uniq_path
                counts = np.ones(len(uniq_path2), dtype=np.int)

            if len(uniq_path2):
                for pth in uniq_path2:
                    t = (pth,)
                    res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                            "WHERE Id = ? ", t).fetchall()
                    pathwayName = pathwayName + [r[1] for r in res]
                    externalId = externalId + [r[2] for r in res]
                    source = source + [r[3] for r in res]
                    dblink = dblink + [r[4] for r in res]
                #sort
                counts,pathwayName,externalId,source,dblink = \
                    list(zip(*sorted(zip(counts,pathwayName,externalId,source,dblink),reverse=True)))
            else:
                raise PathwayNotFoundException
        return pathwayName,externalId,source,dblink,counts

    def find_common_pathway_genes(self, gene_names):
        """
        For a given gene list, find the pathways containing some of the genes,
        and the corresponding genes contained in each of the pathways 
        """
        pathwayName = dict()
        dblink = dict()
        genes = defaultdict(list)
        counts = defaultdict(int)
        if self.tfdb is not None:
            gene_names = list(set(gene_names))
            thr1 = max(2, math.ceil(len(gene_names)/2))
            thr2 = max(2, math.ceil(math.sqrt(len(gene_names))))
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pids = [r[0] for r in res1]
                    for pid in pids:
                        genes[pid].append(gene_name)
                        counts[pid] += 1
            if len(genes):
                max_count = max(counts.values())
                if max_count >= thr1:
                    #only consider pathways which contain at least thr given genes
                    for pth, ct in counts.items():
                        if ct < thr1:
                            del genes[pth]
                elif max_count >= thr2:
                    for pth, ct in counts.items():
                        if ct < thr2:
                            del genes[pth]
                elif max_count >= 2:
                    #consider pathways which contain at least 2 given genes
                    for pth, ct in counts.items():
                        if ct < 2:
                            del genes[pth]
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
            if len(genes):
                for pth in genes.keys():
                    t = (pth,)
                    res = self.tfdb.execute("SELECT pathwayName, dblink FROM pathwayInfo "
                                            "WHERE Id = ? ", t).fetchone()
                    pathwayName[pth] = res['pathwayName']
                    dblink[pth] = res['dblink']
            else:
                raise PathwayNotFoundException
        return pathwayName, dblink, genes
        
    def find_common_pathway_genes_keyword(self, gene_names, keyword):
        """
        For a given gene list and keyword, find the pathways containing some of the genes,
        and return the corresponding given genes in each of the pathways 
        """
        pathwayName = dict()
        dblink = dict()
        counts = defaultdict(int)
        genes = defaultdict(list)
        if self.tfdb is not None:
            gene_names = list(set(gene_names))
            thr1 = max(2, math.ceil(len(gene_names)/2))
            thr2 = max(2, math.ceil(math.sqrt(len(gene_names))))
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pids = [r[0] for r in res1]
                    for pid in pids:
                        genes[pid].append(gene_name)
                        counts[pid] += 1
            if len(genes):
                max_count = max(counts.values())
                if max_count >= thr1:
                    #only consider pathways which contain at least thr given genes
                    for pth, ct in counts.items():
                        if ct < thr1:
                            del genes[pth]
                elif max_count >= thr2:
                    for pth, ct in counts.items():
                        if ct < thr2:
                            del genes[pth]
                elif max_count >= 2:
                    #consider pathways which contain at least 2 given genes
                    for pth, ct in counts.items():
                        if ct < 2:
                            del genes[pth]
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
            if len(genes):
                regstr = '%' + keyword + '%'
                #pths = genes.keys()
                for pth in list(genes.keys()):
                    t = (pth, regstr)
                    res = self.tfdb.execute("SELECT pathwayName, dblink FROM pathwayInfo "
                                            "WHERE Id = ? AND pathwayName LIKE ?", t).fetchone()
                    if res:
                        pathwayName[pth] = res['pathwayName']
                        dblink[pth] = res['dblink']
                    else:
                        del genes[pth]
            else:
                raise PathwayNotFoundException
            if not len(genes):
                raise PathwayNotFoundException
        return pathwayName, dblink, genes
        
    def find_common_pathway_genes_db(self, gene_names, db_name):
        """
        For a given gene list and db name, find the pathways containing at least two of 
        the genes, and return the corresponding given genes in each of the pathways 
        """
        pathwayName = dict()
        dblink = dict()
        counts = defaultdict(int)
        genes = defaultdict(list)
        if self.tfdb is not None:
            gene_names = list(set(gene_names))
            thr1 = max(2, math.ceil(len(gene_names)/2))
            thr2 = max(2, math.ceil(math.sqrt(len(gene_names))))
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pids = [r[0] for r in res1]
                    for pid in pids:
                        genes[pid].append(gene_name)
                        counts[pid] += 1
            if len(genes):
                max_count = max(counts.values())
                if max_count >= thr1:
                    #only consider pathways which contain at least thr1 given genes
                    for pth, ct in counts.items():
                        if ct < thr1:
                            del genes[pth]
                elif max_count >= thr2:
                    for pth, ct in counts.items():
                        if ct < thr2:
                            del genes[pth]
                elif max_count >= 2:
                    #consider pathways which contain at least 2 given genes
                    for pth, ct in counts.items():
                        if ct < 2:
                            del genes[pth]
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
            if len(genes):
                for pth in list(genes.keys()):
                    t = (pth, db_name)
                    res = self.tfdb.execute("SELECT pathwayName, dblink FROM pathwayInfo "
                                            "WHERE Id = ? AND source LIKE ?", t).fetchone()
                    if res:
                        pathwayName[pth] = res['pathwayName']
                        dblink[pth] = res['dblink']
                    else:
                        del genes[pth]
            else:
                raise PathwayNotFoundException
            if not len(genes):
                raise PathwayNotFoundException
        return pathwayName, dblink, genes

    def find_pathway_count_genes_keyword(self, gene_names, keyword):
        """
        For a given gene list and keyword, find the pathways containing one of the genes,
        and the frequency of the pathways 
        """
        pathwayName = []
        externalId = []
        source = []
        dblink = []
        counts = []
        if self.tfdb is not None:
            pathlist = []
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT pathwayID FROM pathway2Genes "
                                         "WHERE genesymbol = ? ", t).fetchall()
                if res1:
                    pathlist = pathlist + [r[0] for r in res1]
            #pathway frequency
            uniq_path = list(set(pathlist))
            if len(uniq_path):
                uniq_path2 = []
                for pth in uniq_path:
                    ct = pathlist.count(pth)
                    if ct > 1:
                        #counts.append(ct)
                        uniq_path2.append(pth)
            if len(uniq_path2):
                regstr = '%' + keyword + ' %'
                for pth in uniq_path2:
                    t = (pth, regstr)
                    res = self.tfdb.execute("SELECT * FROM pathwayInfo "
                                            "WHERE Id = ? AND pathwayName LIKE ?", t).fetchall()
                    if res:
                        pathwayName = pathwayName + [r[1] for r in res]
                        externalId = externalId + [r[2] for r in res]
                        source = source + [r[3] for r in res]
                        dblink = dblink + [r[4] for r in res]
                        counts.append(pathlist.count(pth))
                #sort
                if len(counts):
                    counts,pathwayName,externalId,source,dblink = \
                        list(zip(*sorted(zip(counts,pathwayName,externalId,source,dblink),reverse=True)))
                else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
        return pathwayName,externalId,source,dblink,counts

    def find_targets(self,tf_names):
        """
        Return Targets regulated by the tf or tf list
        """
        dbname = dict()
        if self.tfdb is not None:
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
            if len(target_names):
                target_names.sort()
        else:
            target_names = []
        return target_names,dbname

    def find_overlap_targets_tfs_genes(self,tf_names,target_names):
        """
        Return Targets which are regulated by all the given tfs and are 
        also in the given target list
        """
        if self.tfdb is not None:
            targets = []
            for tf_name in tf_names:
                t = (tf_name,)
                res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                        "WHERE TF = ? ", t).fetchall()
                if res:
                    targets = targets + [r[0] for r in res]
                else:
                    raise TFNotFoundException
            #common regulated targets by the given tf list
            com_targets = []
            for target in set(targets):
                if targets.count(target) == len(tf_names):
                    com_targets.append(target)
            #overlap with the given target list
            overlap_targets = list(set(com_targets) & set(target_names))
            if len(overlap_targets):
                overlap_targets.sort()
        else:
            overlap_targets = []
        return overlap_targets

    def find_targets_tissue(self,tf_names, tissue_name):
        """
        Return Targets regulated by the tf list in a given tissue
        """
        target_names = []
        if self.tfdb is not None:
            regstr = '%' + tissue_name + '%'
            t = (tf_names[0], regstr)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                    "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                #check if the TF is in the database
                t = (tf_names[0],)
                res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                    "WHERE TF = ? ", t).fetchall()
                if res:
                    return target_names
                else:
                    raise TFNotFoundException
                
            if len(tf_names) > 1:
                for i in range(1,len(tf_names)):
                    t = (tf_names[i],regstr)
                    res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                            "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        #check if the TF is in the database
                        t = (tf_names[i],)
                        res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                                "WHERE TF = ? ", t).fetchall()
                        if res:
                            return target_names
                        else:
                            raise TFNotFoundException	
            if len(target_names):
                target_names.sort()
        return target_names
 
    def find_targets_1(self,tf_name):
        """
        Return Targets regulated by the tf
        Accept only one tf name
        """
        if self.tfdb is not None:
            t = (tf_name,)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                    "WHERE TF = ? ORDER BY Target", t).fetchall()
            target_names = [r[0] for r in res]
        else:
            target_names = []
        return target_names

    def find_tfs_1(self,target_name):
        """
        Return TFs regulating the target in a given tissue
        Accept only one target
        """
        if self.tfdb is not None:
            t = (target_name,)
            res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
                                    "WHERE Target = ? ORDER BY TF", t).fetchall()
            tf_names = [r[0] for r in res]
        else:
            tf_names = []
        return tf_names

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

    def Is_miRNA_target(self, miRNA_name, target_name):
        """
        Return True if the miRNA regulates the target, and False if not
        query example: Does miR-20b-5p target STAT3?
        """
        if self.tfdb is not None:
             t = (miRNA_name,)
             res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo WHERE mirna LIKE ? ", t).fetchall()
             if not res:
                 raise miRNANotFoundException
             for r in res:
                 if r[0].upper() == target_name.upper():
                     return True
             #check if target_name in the database
             t = (target_name,)
             res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo WHERE target = ? ", t).fetchall()
             if not res:
                 raise TargetNotFoundException       
        return False
        
    def Is_miRNA_target2(self, miRNA_name_dict, target_name):
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
        
    def find_miRNA_target(self, target_names):
        """
        Return miRNAs regulating all the given targets
        example: What microRNAs target STAT3?
        """
        miRNAs = []
        if self.tfdb is not None:
            t = (target_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                    "WHERE target = ? ", t).fetchall()
            if res:
                miRNAs = [r[0] for r in res]
            else:
                raise TargetNotFoundException
            
            if len(target_names)>1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                    "WHERE target = ? ", t).fetchall()
                    if res:
                        miRNAs = list(set(miRNAs) & set([r[0] for r in res]))
                    else:
                        raise TargetNotFoundException
            if len(miRNAs):
                miRNAs.sort()
        return miRNAs
        
    def find_target_miRNA(self, miRNA_name_dict):
        """
        Return Targets regulated by the given miRNAs
        example: What genes does miR-20b-5p target?
        """
        miRNA_mis = dict()
        target_names = []
        miRNA_names = list(miRNA_name_dict.keys())
        if self.tfdb is not None:
            t = (miRNA_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                #raise miRNANotFoundException
                miRNA_mis[miRNA_names[0]] = miRNA_name_dict[miRNA_names[0]]
                return target_names,miRNA_mis 
            if len(miRNA_names)>1:
                for i in range(1, len(miRNA_names)):
                    t = (miRNA_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        #raise miRNANotFoundException
                        miRNA_mis[miRNA_names[i]] = miRNA_name_dict[miRNA_names[i]]
                        return target_names,miRNA_mis
                    if not len(target_names):
                        break
            if len(target_names):
                target_names.sort()       
        return target_names,miRNA_mis
        
    def find_target_miRNA_strength(self, miRNA_name_dict, evidence_strength):
        """
        Return Targets regulated by the given miRNAs
        example: What genes does miR-20b-5p target?
        """
        miRNA_mis = dict()
        target_names = []
        miRNA_names = list(miRNA_name_dict.keys())
        if self.tfdb is not None:
            if evidence_strength == 'strong':
                t = (miRNA_names[0], '%Weak%')
                res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType NOT LIKE ? ", t).fetchall()
                if res:
                    target_names = [r[0] for r in res]
                else:
                    #raise miRNANotFoundException
                    miRNA_mis[miRNA_names[0]] = miRNA_name_dict[miRNA_names[0]]
                    return target_names,miRNA_mis
                if len(miRNA_names)>1:
                    for i in range(1, len(miRNA_names)):
                        t = (miRNA_names[i], '%Weak%')
                        res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType NOT LIKE ?", t).fetchall()
                        if res:
                            target_names = list(set(target_names) & set([r[0] for r in res]))
                        else:
                            #raise miRNANotFoundException
                            miRNA_mis[miRNA_names[i]] = miRNA_name_dict[miRNA_names[i]]
                            return target_names,miRNA_mis
                        if not len(target_names):
                            break
            else:
                t = (miRNA_names[0], '%Weak%')
                res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType LIKE ? ", t).fetchall()
                if res:
                    target_names = [r[0] for r in res]
                else:
                    #raise miRNANotFoundException
                    miRNA_mis[miRNA_names[0]] = miRNA_name_dict[miRNA_names[0]]
                    return target_names,miRNA_mis
                if len(miRNA_names)>1:
                    for i in range(1, len(miRNA_names)):
                        t = (miRNA_names[i], '%Weak%')
                        res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna LIKE ? AND supportType LIKE ?", t).fetchall()
                        if res:
                            target_names = list(set(target_names) & set([r[0] for r in res]))
                        else:
                            #raise miRNANotFoundException
                            miRNA_mis[miRNA_names[i]] = miRNA_name_dict[miRNA_names[i]]
                            return target_names,miRNA_mis
                        if not len(target_names):
                            break
            if len(target_names):
                target_names.sort()       
        return target_names,miRNA_mis
        
    def find_tf_miRNA(self, miRNA_names):
        """
        Return TFs regulated by the given miRNAs
        example: what transcription factors does miR-124-3p regulate?
        miRNA_names: dict
        """
        target_names,miRNA_mis = self.find_target_miRNA(miRNA_names)
        tf_names = []
        if not miRNA_mis:
            if len(target_names):
                if not self.trans_factor:
                    if self.tfdb is not None:
                        res = self.tfdb.execute("SELECT DISTINCT tf FROM transFactor").fetchall()
                        self.trans_factor = set([r[0] for r in res])
                tf_names = list(set(target_names) & self.trans_factor)
            else:
                raise TFNotFoundException
            if len(tf_names):
                tf_names.sort()
            else:
                raise TFNotFoundException
        else:
            return tf_names,miRNA_mis
        return tf_names,miRNA_mis
        
    def find_evidence_miRNA_target(self, miRNA_name_dict, target_name):
        """
        miRNA_name: str
        Return experimental evidence that the given miRNA regulates the target
        example: What is the evidence that miR-148a-3p targets DNMT1?
        """
        miRNA_mis = dict()
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

    def find_miRNA_count_gene(self, gene_names):
        """
        For a given gene list, find the miRNAs regulate one of the genes, and the
        frequency of the miRNAs
        What mirs most frequently or commonly regulate a list of genes
        """
        mirna = []
        counts = []
        temp = []
        targets = defaultdict(list)
        if self.tfdb is not None:
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                         "WHERE target = ? ", t).fetchall()
                if res1:
                    temp = temp + [r[0] for r in res1]
                    for m in [r[0] for r in res1]:
                        targets[m].append(gene_name)
            #handle the case that all input targets are not in the database
            if not len(temp):
                raise TargetNotFoundException 
            else:               
                #miRNA count
                um = list(set(temp))
                if len(gene_names)>1:
                    for u in um:
                        ct = temp.count(u)
                        if ct > 1:
                            mirna.append(u)
                            counts.append(ct)
                    #sort
                    if len(mirna) > 1:
                        counts, mirna = list(zip(*sorted(zip(counts,mirna), reverse=True)))
                else:
                    mirna = um
                    counts = np.ones(len(um), dtype=np.int)
            if not len(mirna):
                raise miRNANotFoundException
        return mirna,counts,targets
                         
    def find_gene_count_miRNA(self, miRNA_name_dict):
        """
        For a given miRNA list, find the genes regulated by one of the miRNAs, and the
        frequency of the genes
        What genes are most frequently (or commonly) regulated by a list of mIRs?
        """
        genes = []
        counts = dict()
        temp = []
        mrna = defaultdict(list)
        miRNA_mis = dict()
        miRNA_names = list(miRNA_name_dict.keys())
        if self.tfdb is not None:
            for mir in miRNA_names:
                t = (mir,)
                res1 = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                         "WHERE mirna LIKE ? ", t).fetchall()
                if res1:
                    temp = temp + [r[0] for r in res1]
                    for m in [r[0] for r in res1]:
                        mrna[m].append(mir)
                else:
                     miRNA_mis[mir] = miRNA_name_dict[mir]
            if temp:        
                ug = list(set(temp))
                if len(miRNA_names) > 1:
                    for u in ug:
                        ct = temp.count(u)
                        if ct > 1:
                            genes.append(u)
                            counts[u] = ct
                else:
                    genes = ug
                    for u in ug:
                        counts[u] = 1
                    #counts = np.ones(len(ug), dtype=np.int)
                #sort
                #if len(genes) > 1:
                    #counts, genes = list(zip(*sorted(zip(counts,genes), reverse=True)))
        return genes,counts,mrna,miRNA_mis

    def find_pathway_db_keyword(self, db_source, pathway_names):
        """
        Return a list of pathways which come from given db and whose name containing the given keywords
        Which reactome pathways involve immune signaling?
        """
        pathwayName = dict()
        pw_link = dict()
        if self.tfdb is not None:
            pn = []
            pids = []
            plink = []
            for pathway_name in pathway_names:
                regstr = '%' + pathway_name + '%'
                t = (regstr, db_source)
                res = self.tfdb.execute("SELECT Id,pathwayName,dblink FROM pathwayInfo "
                                        "WHERE pathwayName LIKE ? AND source LIKE ?", t).fetchall()
                if res:
                    pids = pids + [r[0] for r in res]
                    pn = pn + [r[1] for r in res]
                    plink = plink + [r[2] for r in res]
            if len(pids):
                pathwayId = list(set(pids))
                for i in range(len(pids)):
                    pathwayName[pids[i]] = pn[i]
                    pw_link[pids[i]] = plink[i]
            else:
                raise PathwayNotFoundException
        return pathwayId,pathwayName,pw_link

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
        return a list of miRNAs beginning with the given miRNA which 
        is not in the database, for user clarification purpose
        miRNA_name is the given miRNA name
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
                raise miRNANotFoundException
        return clari_miRNA
        
    def get_similar_miRNAs2(self, miRNA_names):
        """
        For each miRNA, return a list of miRNAs beginning with the given miRNA which 
        is not in the database, for user clarification purpose
        miRNA_names is a list of given miRNA name
        """
        clari_miRNA = dict()
        if self.tfdb is not None:
            for miRNA_name in miRNA_names:
                regstr = miRNA_name + '-%'
                t = (regstr,)
                res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                    "WHERE mirna LIKE ? ", t).fetchall()
                if res:
                    clari_miRNA[miRNA_name] = [r[0] for r in res]
            if not clari_miRNA:
                raise miRNANotFoundException
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
            if len(kinase_names):
                kinase_names.sort()
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
                    t = (target_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT kinase FROM kinaseReg "
                                            "WHERE target = ? AND direction LIKE ? ", t).fetchall()
                    if res:
                        kinase_names = list(set(kinase_names) & set([r[0] for r in res]))
                    else:
                        raise KinaseNotFoundException
            if len(kinase_names):
                kinase_names.sort()
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
        #gene_names = list(set(gene_names) & hgnc_genes_set)
        if len(gene_names):
            gene_names.sort()
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
        #gene_names = list(set(gene_names) & hgnc_genes_set)
        if len(gene_names):
            gene_names.sort()
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
            #print('tissues: ', ','.join(tissues))
            #print('number of tissues: ', len(tissues))
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
                
    def get_onto_set(self, go_name):
        """
        Return the genes which are in the category of go_name
        """
        go_genes = []
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
                    go_genes += [r[0] for r in res]
                go_genes = set(go_genes)
        return go_genes
        
    def get_hgnc_mapping(self):
        return hgnc_symbol_id
        
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
            return []
        statements = []
        try:
            for stype in stmt_types:
                stmts = get_statements(subject=subj, object=obj, stmt_type=stype, simple_response=True)
                if len(stmts):
                    stmts_filtered = filter_evidence_source(stmts, ['reach'], policy='one')
                    if len(stmts_filtered):
                        statements += stmts_filtered
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
        mirnas = nontfs.intersection(mirna_indra_set)
        others = nontfs - mirnas
        genes = others.intersection(hgnc_genes_set)
        other = others - genes
        return tfs, genes, mirnas, other
        
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
        stmts: dict
        return the list of genes, as well as the filtered statements
        """
        regulators = list(stmts_d.keys())
        stmt_list = []
        for r in regulators:
            stmt_list += stmts_d[r]
        genes = self.find_target_indra(stmts_d[regulators[0]])
        if of_those:
            genes = genes.intersection(of_those)
        target_type_set = self.get_onto_set(target_type)
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
            stmt_list += stmts_d[t]
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
            
    
        
#test functions
#if __name__ == "__main__":
    #a = TFTA()
    
