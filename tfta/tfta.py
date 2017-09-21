#TFTA stands for TF-target agent whose task is to search for targets known to a TF and vice versa

import re
import os
import logging
import sqlite3
import numpy as np

 
logger = logging.getLogger('TFTA')
 
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
 	
class TFTA:
    def __init__(self):
        logger.debug('Using resource folder: %s' % _resource_dir)
 	#Load TF_target database
        tf_db_file = _resource_dir + 'TF_target_v5.db'
        if os.path.isfile(tf_db_file):
            self.tfdb = sqlite3.connect(tf_db_file, check_same_thread=False)
	    logger.info('Loaded TF-target database')
        else:
            logger.error('TFTA could not load TF-target database.')
            self.tfdb = None;
		
 			
    def __del__(self):
        self.tfdb.close()
 				
 		
    def Is_tf_target(self,tf_name,target_name):
        """
        Return True if the tf regulates the target, and False if not
        """
        if self.tfdb is not None:
            #check if target_name in the database
            t = (target_name,)
            res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB WHERE Target = ? ", t).fetchall()
            if not res:
                raise TargetNotFoundException
            else:
                t = (tf_name,)
                res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB WHERE TF = ? ", t).fetchall()
                if not res:
                    raise TFNotFoundException
                for r in res:
                    if r[0].upper() == target_name.upper():
                        return True
        return False
 		
    def Is_tf_target_tissue(self,tf_name,target_name,tissue_name):
        """
        Return True if the tf regulates the target in a given tissue, and False if not
        """
        if self.tfdb is not None:
            t = (tf_name, tissue_name)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                    "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
            if not res:
                raise TFNotFoundException
            for r in res:
                if r[0].upper() == target_name.upper():
                    return True   
            #check if target_name in database
            t = (target_name, tissue_name)
            res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                    "WHERE Target = ? AND Tissue LIKE ? ", t).fetchall()
            if not res:
                raise TargetNotFoundException
        return False
 		

    def find_tfs(self,target_names):
        """
        Return TFs regulating all the given targets
        """		
 	#query
        if self.tfdb is not None:
            t = (target_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
                                    "WHERE Target = ? ", t).fetchall()
            if res:
                tf_names = [r[0] for r in res]
            else:
                raise TargetNotFoundException 
            if len(target_names)>1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
                                            "WHERE Target = ? ", t).fetchall()
                    if res:
                        tf_names = list(set(tf_names) & set([r[0] for r in res]))
                    else:
                        raise TargetNotFoundException
            if len(tf_names):
                tf_names.sort()
        else:
            tf_names = []
	#print tf_names
        return tf_names
 		
    def find_tfs_count(self,target_names):
        """
        Return TFs regulating the given target list as well as the frequency (count)
	of each TF
	"""		
 	#query
	tf_counts = []
	tf_names = []
	unique_tfs = []
        if self.tfdb is not None:
            for target_name in target_names:
                t = (target_name,)
                res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
                                        "WHERE Target = ? ", t).fetchall()
                if res:
                    tf_names = tf_names + [r[0] for r in res]
	    #count for each TF
	    if len(tf_names):
                unique_tfs = list(set(tf_names))
                for i in range(len(unique_tfs)):
                    tf_counts.append(tf_names.count(unique_tfs[i]))		
 	        #sort
                tf_counts,unique_tfs = zip(*sorted(zip(tf_counts, unique_tfs),reverse=True))
	    else:
		raise TFNotFoundException	
        return unique_tfs, tf_counts
 		
    def find_tfs_tissue(self,target_names,tissue_name):
        """
        Return TFs regulating the targets in a given tissue
        """		
 	#query
        if self.tfdb is not None:
            t = (target_names[0],tissue_name)
            res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                    "WHERE Target = ? AND Tissue LIKE ? ", t).fetchall()
            if res:
                tf_names = [r[0] for r in res]
            else:
                raise TargetNotFoundException
 			
            if len(target_names) > 1:
                for i in range(1,len(target_names)):
                    t = (target_names[i],tissue_name)
                    res = self.tfdb.execute("SELECT DISTINCT TF FROM Target2TF2Tissue "
                                            "WHERE Target = ? AND Tissue LIKE ? ", t).fetchall()
                    if res:
                        tf_names = list(set(tf_names) & set([r[0] for r in res]))
                    else:
                        raise TargetNotFoundException
            if len(tf_names):			
                tf_names.sort()		
        else:
            tf_names = []			
        return tf_names
 		
    def find_pathways_from_name(self, pathway_name):
        """
        return pathway information related to pathway_name
        """
 	#query
        if self.tfdb is not None:
            regstr = '%' + pathway_name + '%'
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
	#sort
	#pathwayName,pathwayId,externalId,source,dblink = \
	    #zip(*sorted(zip(pathwayName,pathwayId,externalId,source,dblink)))
        return pathwayId,pathwayName,externalId,source,dblink

    def find_pathways_from_dbsource_geneName1(self, dbsource,gene_name):
        """
        return pathway information for given dbsource and gene_name
        """
 	#query
        pathwayId = []
        pathwayName = []
        externalId = []
        psource = []
        dblink = []
        if self.tfdb is not None:
 	    #regstr='%'+pathway_name+'%'
            t = (gene_name[0],)
 	    #get pathwayId
            #res = self.tfdb.execute("SELECT * FROM pathwayInfo "
            #                       "WHERE Id in (SELECT DISTINCT pathwayID FROM pathway2Genes "
              #                      "WHERE genesymbol = ?) AND source LIKE ? ", t).fetchall()
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
	            zip(*sorted(zip(pathwayName,pathwayId,externalId,psource,dblink)))
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
	            zip(*sorted(zip(pathwayName,pathwayId,externalId,psource,dblink)))				
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
 	
    def find_tfs_from_pathwayName(self,pathway_names):
        """
        Return TFs within the given pathway
        """		
 	pathwayId = []
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
			pathwayId.remove(pthID)
            else:
                raise PathwayNotFoundException				
        return pathwayId,pathwayName,tflist,dblink
 		
    def find_pathways_from_genelist(self,gene_names):
        """
        return pathways having given genes
        """
 	#query
        pathwayId=[]
        pathwayName=[]
        #externalId=[]
        #source=[]
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
	            zip(*sorted(zip(pathwayName,dblink)))
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
            for pthID in upid:
                t = (pthID,)
                res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                        "WHERE pathwayID = ? ORDER BY genesymbol", t).fetchall()
                genes = [r[0] for r in res1]
                overlap_genes = list(set(gene_names) & set(genes))
                if len(overlap_genes) < len(gene_names):
                    del pname[pthID]
                    del plink[pthID]
                    upid.remove(pthID)    
            return upid, pname, plink

    def find_pathways_from_genelist_keyword(self,gene_names, keyword):
        """
        return pathways having given genes and some information in pathway name
        """
 	#query
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
	                zip(*sorted(zip(pathwayName1,pathwayId,externalId,source,dblink)))
		else:
		    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException    
        return pathwayId,pathwayName1,externalId,source,dblink
 		
    def find_pathways_from_chemical(self, chemical_name):
        """
        return pathways containing the given chemical
        """
 	#query
        if self.tfdb is not None:
            regstr = '%' + chemical_name + '%'
            t = (regstr,)
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
 		
    def find_tfs_from_pathwaysWithChemical(self, chemical_name):
        """
        return pathways containing the given chemical and related tfs 
        """
 	#query
        if self.tfdb is not None:
            regstr = '%' + chemical_name + '%'
            t = (regstr,)
            res = self.tfdb.execute("SELECT Id,pathwayName FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
            if res:
                pathwayId = [r[0] for r in res]
                pathwayName = [r[1] for r in res]
 				
 		#search tfs
                tflist = dict()
                newpathwayName = []
                newpathwayId = []
                for pn,pthID in zip(pathwayName,pathwayId):
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                             "WHERE pathwayID = ? AND isTF = 1 ORDER BY genesymbol", t).fetchall()
                    if res1:
                        tfs = [r[0] for r in res1]
                        tflist[pthID] = tfs
                        newpathwayId.append(pthID)
                        newpathwayName.append(pn)     		
            else:
                raise PathwayNotFoundException		
        else:
            newpathwayId = []
            newpathwayName = []
            tflist = dict()		
        return newpathwayId,newpathwayName,tflist
 		
    def find_pathway_count_genes(self, gene_names):
        """
        For a given gene list, find the pathways containing one of the genes,
 	and the frequency of the pathways 
 	"""
 	#query
 	#pathwayId=[]
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
 		    #pathwayId=pathwayId+[r[0] for r in res]
                    pathwayName = pathwayName + [r[1] for r in res]
                    externalId = externalId + [r[2] for r in res]
                    source = source + [r[3] for r in res]
                    dblink = dblink + [r[4] for r in res]			
 		#sort
                counts,pathwayName,externalId,source,dblink = \
		    zip(*sorted(zip(counts,pathwayName,externalId,source,dblink),reverse=True))
            else:
                raise PathwayNotFoundException		
        return pathwayName,externalId,source,dblink,counts

    def find_pathway_count_genes_keyword(self, gene_names, keyword):
        """
        For a given gene list and keyword, find the pathways containing one of the genes,
 	and the frequency of the pathways 
 	"""
 	#query
 	#pathwayId=[]
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
		regstr = '%' + keyword + '%'
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
		        zip(*sorted(zip(counts,pathwayName,externalId,source,dblink),reverse=True))
		else:
                    raise PathwayNotFoundException
            else:
                raise PathwayNotFoundException
        return pathwayName,externalId,source,dblink,counts
 			
    def find_targets(self,tf_names):
        """
        Return Targets regulated by the tf or tf list
        """
 	#query
        if self.tfdb is not None:
            t = (tf_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                    "WHERE TF = ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                raise TFNotFoundException
                
            if len(tf_names) > 1:
                for i in range(1,len(tf_names)):
                    t = (tf_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
                                            "WHERE TF = ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        raise TFNotFoundException
 	    if len(target_names):	
                target_names.sort()	
        else:
            target_names = []	
        #print "target_names=",target_names
        return target_names
 		
    def find_overlap_targets_tfs_genes(self,tf_names,target_names):
        """
        Return Targets which are regulated by all the given tfs and are 
	also in the given target list
        """
 	#query
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
 	#query
        if self.tfdb is not None:
            t = (tf_names[0],tissue_name)
            res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                    "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                raise TFNotFoundException
                
            if len(tf_names) > 1:
                for i in range(1,len(tf_names)):
                    t = (tf_names[i],tissue_name)
                    res = self.tfdb.execute("SELECT DISTINCT Target FROM Target2TF2Tissue "
                                            "WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        raise TFNotFoundException	
            target_names.sort()		
        else:
            target_names = []		
        return target_names
 			
    def find_targets_1(self,tf_name):
        """
        Return Targets regulated by the tf
        Accept only one tf name
        """
 	#query
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
 	#query
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
            #print 'res=',res
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
        query example: Does miR-20b-5p target STAT3? miRNANotFoundException
        """
        if self.tfdb is not None:
             t = (miRNA_name,)
             res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo WHERE mirna = ? ", t).fetchall()
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
        
    def find_target_miRNA(self, miRNA_names):
        """
        Return Targets regulated by the given miRNAs
        example: What genes does miR-20b-5p target?
        """
        target_names = []
        if self.tfdb is not None:
            t = (miRNA_names[0],)
            res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna = ? ", t).fetchall()
            if res:
                target_names = [r[0] for r in res]
            else:
                raise miRNANotFoundException
            if len(miRNA_names)>1:
                for i in range(1, len(miRNA_names)):
                    t = (miRNA_names[i],)
                    res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                    "WHERE mirna = ? ", t).fetchall()
                    if res:
                        target_names = list(set(target_names) & set([r[0] for r in res]))
                    else:
                        raise miRNANotFoundException
            if len(target_names):
                target_names.sort()       
        return target_names
        
    def find_evidence_miRNA_target(self, miRNA_name, target_name):
        """
        Return experimental evidence that the given miRNA regulates the target
        example: What is the evidence that miR-148a-3p targets DNMT1?
        """
        experiments = []
        support_types = []
        pmid_link = []
        if self.tfdb is not None:
            t = (miRNA_name, target_name)
            res = self.tfdb.execute("SELECT * FROM mirnaInfo "
                                    "WHERE mirna = ? AND target = ? ", t).fetchall()
            if res:
                experiments = [r[3] for r in res]
                support_types = [r[4] for r in res]
                pmid_link = [pmid_sublink+str(r[5]) for r in res]
            else:
                #check if miRNA_name in the database
                t = (miRNA_name, )
                res = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                        "WHERE mirna = ? ", t).fetchall()
                if not res:
                    raise miRNANotFoundException
                else:
                    #check if target_name in the database
                    t = (target_name,)
                    res = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                            "WHERE target = ? ", t).fetchall()
                    if not res:
                        raise TargetNotFoundException            
        return experiments, support_types, pmid_link

    def find_miRNA_count_gene(self, gene_names):
        """
        For a given gene list, find the miRNAs regulate one of the genes, and the
        frequency of the miRNAs
        What mirs most frequently or commonly regulate a list of genes
        """
        mirna = []
        counts = []
        temp = []
        if self.tfdb is not None:
            for gene_name in gene_names:
                t = (gene_name,)
                res1 = self.tfdb.execute("SELECT DISTINCT mirna FROM mirnaInfo "
                                         "WHERE target = ? ", t).fetchall()
                if res1:
                    temp = temp + [r[0] for r in res1]
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
                        counts,mirna = zip(*sorted(zip(counts,mirna), reverse=True))
                else:
                    mirna = um
                    counts = np.ones(len(um), dtype=np.int)
            if not len(mirna):
                raise miRNANotFoundException
        return mirna,counts
                         
    def find_gene_count_miRNA(self, miRNA_names):
        """
        For a given miRNA list, find the genes regulated by one of the miRNAs, and the
        frequency of the genes
        What genes are most frequently (or commonly) regulated by a list of mIRs?
        """
        genes = []
        counts = []
        temp = []
        if self.tfdb is not None:
            for mir in miRNA_names:
                t = (mir,)
                res1 = self.tfdb.execute("SELECT DISTINCT target FROM mirnaInfo "
                                         "WHERE mirna = ? ", t).fetchall()
                if res1:
                    temp = temp + [r[0] for r in res1]
            #handle the case that all input miRNAs are not in the database
            if not len(temp):
                raise miRNANotFoundException
            else:        
                #gene count
                ug = list(set(temp))
                if len(miRNA_names)>1:
                    for u in ug:
                        ct = temp.count(u)
                        if ct > 1:
                            genes.append(u)
                            counts.append(ct)
                    #sort
                    if len(genes) > 1:
                        counts,genes = zip(*sorted(zip(counts,genes), reverse=True))  
                else:
                    genes = ug
                    counts = np.ones(len(ug), dtype=np.int)
                    
            if not len(genes):
                raise TargetNotFoundException     
        return genes,counts

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
		 							 

#test functions
#if __name__ == "__main__":
    #a=TFTA()	
