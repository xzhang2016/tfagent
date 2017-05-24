#TFTA stands for TF-target agent whose task is to search for targets known to a TF and vice versa

import re
import os
import logging
import sqlite3

 
logger = logging.getLogger('TFTA')
 
_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/../resources/'
 
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
 	
class TFTA:
    def __init__(self):
        logger.debug('Using resource folder: %s' % _resource_dir)
 	#Load TF_target database
        tf_db_file = _resource_dir + 'TF_target_v3.db'
        if os.path.isfile(tf_db_file):
            self.tfdb = sqlite3.connect(tf_db_file, check_same_thread=False)
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
 	
	print tf_names
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
            unique_tfs = list(set(tf_names))
	    if len(unique_tfs):
                for i in range(len(unique_tfs)):
                    tf_counts.append(tf_names.count(unique_tfs[i]))
 				
 	        #sort
                tf_counts,unique_tfs = zip(*sorted(zip(tf_counts, unique_tfs),reverse=True))
 			
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
 			
 		
    def find_genes_from_pathwayName(self, pathway_name):
        """
        return genes related to pathway_name
        """
 	#query
        if self.tfdb is not None:
            regstr = '%' + pathway_name + '%'
            t = (regstr,)
 	    #get pathwayId
            res = self.tfdb.execute("SELECT Id,pathwayName FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
            if res:
                pathwayId = [r[0] for r in res]
                pathwayName = [r[1] for r in res]
 				
                genelist = dict()
                for pthID in pathwayId:
                    t = (pthID,)
                    res1 = self.tfdb.execute("SELECT DISTINCT genesymbol FROM pathway2Genes "
                                           "WHERE pathwayID = ? ORDER BY genesymbol", t).fetchall()
                    genes = [r[0] for r in res1]
                    genelist[pthID] = genes
 					
            else:
                raise PathwayNotFoundException
 				
        else:
            pathwayId = []
            pathwayName = []
            genelist = dict()
 			
        return pathwayId,pathwayName,genelist
 	
    def find_tfs_from_pathwayName(self,pathway_name):
        """
        Return TFs within the given pathway
        """
 			
 	#query
        if self.tfdb is not None:
            regstr = '%' + pathway_name + '%'
            t = (regstr,)
 	    #get pathwayId and TF symbols
            res = self.tfdb.execute("SELECT Id,pathwayName FROM pathwayInfo "
                                    "WHERE pathwayName LIKE ? ", t).fetchall()
 			
            if res:
                pathwayId = [r[0] for r in res]
                pathwayName = [r[1] for r in res]
 			
                tflist = dict()
                newpathwayId = []
                newpathwayName = []
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
 		
    def find_pathways_from_genelist(self,gene_names):
        """
        return pathways having given genes
        """
 	#query
        pathwayId=[]
        pathwayName=[]
        externalId=[]
        source=[]
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
            pathIDs=[]
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
                    source = source + [r[3] for r in res]
                    dblink = dblink + [r[4] for r in res]
		
		#sort
	        pathwayName,pathwayId,externalId,source,dblink = \
	            zip(*sorted(zip(pathwayName,pathwayId,externalId,source,dblink)))
            else:
                raise PathwayNotFoundException
 			
        return pathwayId,pathwayName,externalId,source,dblink

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
	    if len(uniq_path):
	        uniq_path2 = []
                for pth in uniq_path:
                    ct = pathlist.count(pth)
		    if ct > 1:
                        counts.append(ct)
		        uniq_path2.append(pth)
 			
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
 			
        print "target_names=",target_names
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
                go_ids = [r[0] for r in res]
                go_types = [r[1] for r in res]
                go_names = [r[2] for r in res]
                
                #search genes
                for i in range(len(go_ids)):
                    t = (go_ids[i],)
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
		 							 

#test functions
if __name__ == "__main__":
	a=TFTA()
	"""
	ss2 = a.Is_tf_target('FOS','ELF3');
	if ss2 is not None:
		print 'ss2='+str(ss2)
	else:
		print 'ss2 is none'
	
	ss3 = a.find_tfs(['KRAS','USP1']) 
	if ss3 is not None:
		print 'lenth(ss3)='+str(len(ss3))
		print 'ss3='+','.join(ss3)
		
	ss4 = a.find_targets(['FOS','STAT3'])
        print 'length(ss4)='+str(len(ss4))
        print 'ss4='+','.join(ss4)
		
	ss6=a.Is_tf_target_tissue('BSAP','USP1','bladder')
        print 'ss6='+str(ss6)
		
	ss7=a.find_tfs_tissue(['MAPK14','SYPL1'],'bladder')
        print 'lenth(ss7)='+str(len(ss7))
        print 'ss7='+','.join(ss7)
	
	ss8=a.find_targets_tissue(['SREBP-1','USF'],'bladder')
        print 'length(ss8)=', len(ss8)
        print 'ss8='+','.join(ss8)
 	
 	pathId,pathName,extId,pathsource,pathlink=a.find_pathways_from_name('MAPK signaling pathway')
 	print 'pathId='+','.join(str(i) for i in pathId)
 	#imap
 	print 'pathId='+','.join(map(str,pathId))
 	print 'pathName='+','.join(pathName)
 	print 'extId='+','.join(extId)
 	print 'pathsource='+','.join(pathsource)
 	print 'pathlink='+','.join(pathlink)
 	
 	pathId,pathName,genes=a.find_tfs_from_pathwayName('MAPK signaling pathway')
 	print 'pathId='+','.join(map(str,pathId))
 	print 'pathName='+','.join(pathName)
 	print genes
 	print genes.keys()
 	
 	pathId,pathName,extId,pathsource,pathlink=a.find_pathways_from_dbsource_geneName('reactome','SRF')
 	print 'pathId='+','.join(map(str,pathId))
 	print 'pathName='+','.join(pathName)
 	print 'extId='+','.join(extId)
 	print 'pathsource='+','.join(pathsource)
 	print 'pathlink='+','.join(pathlink)
 	
 	pathId,pathName,genelist=a.find_genes_from_pathwayName('MAPK signaling pathway')
 	print 'pathId='+','.join(map(str,pathId))
 	print genelist
 	print genelist.keys()
 	
 	pathId,pathName,extId,pathsource,pathlink=a.find_pathways_from_genelist(['KRAS','ELK1'])
 	print 'pathId='+','.join(map(str,pathId))
 	print 'pathName='+','.join(pathName)
 	print 'extId='+','.join(extId)
 	print 'pathsource='+','.join(pathsource)
 	print 'pathlink='+','.join(pathlink)
 	
 	pathId,pathName,extId,pathsource,pathlink=a.find_pathways_from_chemical('calcium')
 	print 'pathId='+','.join(map(str,pathId))
 	print 'pathName='+','.join(pathName)
 	print 'extId='+','.join(extId)
 	print 'pathsource='+','.join(pathsource)
 	print 'pathlink='+','.join(pathlink)
 	
 	pathId,pathName,genes=a.find_tfs_from_pathwaysWithChemical('calcium')
 	print 'pathId='+','.join(map(str,pathId))
 	print 'pathName='+','.join(pathName)
 	print genes
 	print genes.keys()
 	
 	
 	targets=['AK5','ALDH1L2','ALPK2','ATOH8','BEGAIN','C14orf37','CCDC858','CCL28','CD2A',
 			'CEMP','COL1A1','COL1A2','COL5A1','COL5A2','COL8A1','CPA4','CRYAB','CTGF','CYP181',
 			'ELFN2','EMILIN2','EPHX1','F2R','F2RL2','FBN1','FBX032','FIBCD1','GALNT5','GREM1',
 			'HHIPL2','HSP86','IGFBP3','IGFBP4','IGFBP7','KIF268','LOX','LOXLl','LOXL3','LRRC15',
 			'LYPD1','MAGED1','MDGA1','MFAP4','MOK','MYLK','NI02','OLFML3','PADI2','PAMR1',
 			'PCOLCE2','PLXDC1','PLXNA4','PTK7','PTPRS','RASD2','SAMD11','SIPA1L2','SLC30A4',
 			'SLC7A8','SMADG','SPARC','SSC50','STAC2','STC2','SYNC','TENM2','TGFBI','THBS1',
 			'TNFRSF21','TNFRSF9','TUBA1A']
 	
 	tfs,counts=a.find_tfs_count(targets)
 	for i in range(len(tfs)):
 		print tfs[i]+'('+str(counts[i])+')'
 	print 'lenth(count)='+str(len(counts))
 	print 'length(tfs)='+str(len(tfs))
 	
 	targetlist=a.find_overlap_targets_tfs_genes(['IRF3','MAX'],targets)
 	print 'targetlist='+','.join(targetlist)
 	print 'length(targetlist)='+str(len(targetlist))
 	
 	pathwayName,externalId,source,dblink,counts=a.find_pathway_count_genes(targets)
 	for pn, ct in zip(pathwayName, counts):
 	    print pn, '(',ct,')'
	"""
 	
 	
