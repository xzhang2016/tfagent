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
 	
class TFTA:
	def __init__(self):
		logger.debug('Using resource folder: %s' % _resource_dir)
 		#Load TF_target database
 		tf_db_file = _resource_dir + 'TF_target.db'
 		if os.path.isfile(tf_db_file):
 			self.tfdb = sqlite3.connect(tf_db_file, check_same_thread=False)
 		else:
 			logger.error('TFTA could not load TF-target database.')
 			self.tfdb = None;
		
 			
	def __del__(self):
		self.tfdb.close()
 				
 		
	def Is_tf_target(self,tf_name,target_name):
 		'''
 		Return True if the tf regulates the target, and False if not
 		'''
 		if self.tfdb is not None:
 			t = (tf_name,)
 			res = self.tfdb.execute("SELECT Target FROM CombinedDB WHERE TF = ? ", t).fetchall()
 			if not res:
 				raise TFNotFoundException
 			for r in res:
 				if r[0].upper() == target_name.upper():
 					return True
 		
 		return False
 		
 	def Is_tf_target_tissue(self,tf_name,target_name,tissue_name):
 		'''
 		Return True if the tf regulates the target in a given tissue, and False if not
 		'''
 		if self.tfdb is not None:
 			t = (tf_name, tissue_name)
 			res = self.tfdb.execute("SELECT Target FROM Target2TF2Tissue WHERE TF = ? AND Tissue LIKE ? ", t).fetchall()
 			if not res:
 				raise TFNotFoundException
 			for r in res:
 				if r[0].upper() == target_name.upper():
 					return True
 		
 		return False
 		

	def find_tfs(self,target_names):
		'''
		Return TFs regulating all the given targets
		'''
 			
 		#query
 		if self.tfdb is not None:
 			t = (target_names[0],)
 			res = self.tfdb.execute("SELECT TF FROM CombinedDB "
 									"WHERE Target = ?  ORDER BY TF", t).fetchall()
 			tf_names = [r[0] for r in res]
 			
 			if len(target_names)>1:
 				for i in range(1,len(target_names)):
 					t = (target_names[i],)
 					res = self.tfdb.execute("SELECT TF FROM CombinedDB "
 									"WHERE Target = ?  ORDER BY TF", t).fetchall()
 					tf_names = list(set(tf_names) & set([r[0] for r in res]))
 		else:
 			tf_names = []
 				
 		return tf_names
 		
 	def find_tfs_count(self,target_names):
		'''
		Return TFs regulating the given target list as well as the frequency (count)
		of each TF
		'''
 			
 		#query
 		if self.tfdb is not None:
 			tf_names = []
 			for target_name in target_names:
 				t = (target_name,)
 				res = self.tfdb.execute("SELECT DISTINCT TF FROM CombinedDB "
 									"WHERE Target = ? ", t).fetchall()
 				tf_names = tf_names + [r[0] for r in res]
 			#count for each TF
 			unique_tfs=list(set(tf_names))
 			tf_counts=[]
 			for i in range(len(unique_tfs)):
 				tf_counts.append(tf_names.count(unique_tfs[i]))
 				
 			#sort
 			tf_counts,unique_tfs = zip(*sorted(zip(tf_counts, unique_tfs),reverse=True))
 			
 		return unique_tfs, tf_counts
 		
 	def find_tfs_tissue(self,target_name,tissue_name):
		'''
		Return TFs regulating the target in a given tissue
		'''
 			
 		#query
 		if self.tfdb is not None:
 			t = (target_name,tissue_name)
 			res = self.tfdb.execute("SELECT TF FROM Target2TF2Tissue "
 									"WHERE Target = ? and Tissue LIKE ? ORDER BY TF", t).fetchall()
 			tf_names = [r[0] for r in res]
 		else:
 			tf_names = []
 				
 		return tf_names
 		
 	def find_pathways_from_name(self, pathway_name):
 		'''
 		return pathway information related to pathway_name
 		'''
 		#query
 		if self.tfdb is not None:
 			regstr='%'+pathway_name+'%'
 			t = (regstr,)
 			#get pathwayId
 			res = self.tfdb.execute("SELECT * FROM pathwayInfo "
 									"WHERE pathwayName LIKE ? ", t).fetchall()
 			if res is not None:
 				pathwayId=[r[0] for r in res]
 				pathwayName=[r[1] for r in res]
 				externalId=[r[2] for r in res]
 				source=[r[3] for r in res]
 				dblink=[r[4] for r in res]
 			else:
 				raise PathwayNotFoundException
 				
 		else:
 			pathwayId = []
 			pathwayName = []
 			externalId = []
 			source = []
 			dblink = []
 			
 		return pathwayId,pathwayName,externalId,source,dblink
 		
 	def find_pathways_from_dbsource_geneName(self, dbsource,gene_name):
 		'''
 		return pathway information for given dbsource and gene_name
 		'''
 		#query
 		if self.tfdb is not None:
 			#regstr='%'+pathway_name+'%'
 			t = (gene_name,dbsource)
 			#get pathwayId
 			res = self.tfdb.execute("SELECT * FROM pathwayInfo "
 									"WHERE Id in (SELECT DISTINCT pathwayID FROM pathway2Genes "
 									"WHERE genesymbol = ?) AND source LIKE ? ", t).fetchall()
 			if res is not None:
 				pathwayId=[r[0] for r in res]
 				pathwayName=[r[1] for r in res]
 				externalId=[r[2] for r in res]
 				source=[r[3] for r in res]
 				dblink=[r[4] for r in res]
 			else:
 				raise PathwayNotFoundException
 				
 		else:
 			pathwayId = []
 			pathwayName = []
 			externalId = []
 			source = []
 			dblink = []
 			
 		return pathwayId,pathwayName,externalId,source,dblink
 		
 	def find_genes_from_pathwayName(self, pathway_name):
 		'''
 		return genes related to pathway_name
 		'''
 		#query
 		if self.tfdb is not None:
 			regstr='%'+pathway_name+'%'
 			t = (regstr,)
 			#get pathwayId
 			res = self.tfdb.execute("SELECT Id,pathwayName FROM pathwayInfo "
 									"WHERE pathwayName LIKE ? ", t).fetchall()
 			if res is not None:
 				pathwayId=[r[0] for r in res]
 				pathwayName=[r[1] for r in res]
 				'''
 				externalId=[r[2] for r in res]
 				source=[r[3] for r in res]
 				dblink=[r[4] for r in res]
 				'''
 				genelist=dict()
 				for pthID in pathwayId:
 					t=(pthID,)
 					res1=self.tfdb.execute("SELECT genesymbol FROM pathway2Genes "
 									"WHERE pathwayID = ? ", t).fetchall()
 					genes=[r[0] for r in res1]
 					genelist[pthID]=genes
 					
 			else:
 				raise PathwayNotFoundException
 				
 		else:
 			pathwayId = []
 			pathwayName = []
 			#externalId = []
 			#source = []
 			#dblink = []
 			genelist=dict()
 			
 		return pathwayId,pathwayName,genelist
 	
 	def find_tfs_from_pathwayName(self,pathway_name):
		'''
		Return TFs within the given pathway
		'''
 			
 		#query
 		if self.tfdb is not None:
 			regstr='%'+pathway_name+'%'
 			t = (regstr,)
 			#get pathwayId and TF symbols
 			res = self.tfdb.execute("SELECT Id,pathwayName FROM pathwayInfo "
 									"WHERE pathwayName LIKE ? ", t).fetchall()
 			
 			if res is not None:
 				pathwayId=[r[0] for r in res]
 				pathwayName=[r[1] for r in res]
 			
 				tflist=dict()
 				for pthID in pathwayId:
 					t=(pthID,)
 					res1 = self.tfdb.execute("SELECT genesymbol FROM pathway2Genes "
 									"WHERE pathwayID = ? AND isTF=1 ", t).fetchall()
 					tfs=[r[0] for r in res1]
 					tflist[pthID]=tfs
 					
			else:
				raise PathwayNotFoundException
				
 		else:
 			pathwayId = []
 			pathwayName = []
 			tflist = dict()
 				
 		return pathwayId,pathwayName,tflist
 		
 	def find_pathways_from_genelist(self,gene_names):
 		'''
 		return pathways having given genes
 		'''
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
 				pathlist=pathlist+[r[0] for r in res1]
 				#pathIDs=set(pathIDs).intersection(set([r[0] for r in res1]))
 			#intersection
 			pathIDs=[]
 			for pth in set(pathlist):
 				if pathlist.count(pth)==len(gene_names):
 					pathIDs.append(pth)
 			
 			if len(pathIDs)>0:
 				
 				for pth in pathIDs:
 					t=(pth,)	
 					res = self.tfdb.execute("SELECT * FROM pathwayInfo "
 									"WHERE Id = ? ", t).fetchall()
 			
 					pathwayId=pathwayId+[r[0] for r in res]
 					pathwayName=pathwayName+[r[1] for r in res]
 					externalId=externalId+[r[2] for r in res]
 					source=source+[r[3] for r in res]
 					dblink=dblink+[r[4] for r in res]
 			else:
 				raise PathwayNotFoundException
 			
 		return pathwayId,pathwayName,externalId,source,dblink
 		
 	def find_pathways_from_chemical(self, chemical_name):
 		'''
 		return pathways containing the given chemical
 		'''
 		#query
 		if self.tfdb is not None:
 			regstr='%'+chemical_name+'%'
 			t = (regstr,)
 			res = self.tfdb.execute("SELECT * FROM pathwayInfo "
 									"WHERE pathwayName LIKE ? ", t).fetchall()
 			if res is not None:
 				pathwayId=[r[0] for r in res]
 				pathwayName=[r[1] for r in res]
 				externalId=[r[2] for r in res]
 				source=[r[3] for r in res]
 				dblink=[r[4] for r in res]
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
 		'''
 		return pathways containing the given chemical and related tfs 
 		'''
 		#query
 		if self.tfdb is not None:
 			regstr='%'+chemical_name+'%'
 			t = (regstr,)
 			res = self.tfdb.execute("SELECT Id,pathwayName FROM pathwayInfo "
 									"WHERE pathwayName LIKE ? ", t).fetchall()
 			if res is not None:
 				pathwayId=[r[0] for r in res]
 				pathwayName=[r[1] for r in res]
 				
 				#search tfs
 				tflist=dict()
 				for pthID in pathwayId:
 					t=(pthID,)
 					res1 = self.tfdb.execute("SELECT genesymbol FROM pathway2Genes "
 									"WHERE pathwayID = ? AND isTF=1 ", t).fetchall()
 					tfs=[r[0] for r in res1]
 					tflist[pthID]=tfs
 				
 			else:
 				raise PathwayNotFoundException
 				
 		else:
 			pathwayId = []
 			pathwayName = []
 			tflist=dict()
 			
 		return pathwayId,pathwayName,tflist
 			
 	def find_targets(self,tf_names):
 		'''
		Return Targets regulated by the tf or tf list
		'''
		
 		#query
 		if self.tfdb is not None:
 			t = (tf_names[0],)
 			res = self.tfdb.execute("SELECT Target,TargetEntrezID FROM CombinedDB "
 									"WHERE TF = ?  ORDER BY Target", t).fetchall()
 			target_names = [r[0] for r in res]
 			if len(tf_names) > 1:
 				for i in range(1,len(tf_names)):
 					t = (tf_names[i],)
 					res = self.tfdb.execute("SELECT Target,TargetEntrezID FROM CombinedDB "
 									"WHERE TF = ?  ORDER BY Target", t).fetchall()
 					target_names = list(set(target_names) & set([r[0] for r in res]))
 			#targetEntrez = [r[1] for r in res]
 		else:
 			target_names = []
 			#targetEntrez = []
 				
 		return target_names
 		
 	def find_overlap_targets_tfs_genes(self,tf_names,target_names):
 		'''
		Return Targets which are regulated by all the given tfs and are 
		also in the given target list
		'''
		
 		#query
 		if self.tfdb is not None:
 			targets = []
 			for tf_name in tf_names:
 				t = (tf_name,)
 				res = self.tfdb.execute("SELECT DISTINCT Target FROM CombinedDB "
 									"WHERE TF = ? ", t).fetchall()
 				targets = targets + [r[0] for r in res]
 				
 			#common regulated targets by the given tf list
 			com_targets = []
 			for target in set(targets):
 				if targets.count(target) == len(tf_names):
 					com_targets.append(target)
 					
 			#overlap with the given target list
 			overlap_targets = list(set(com_targets) & set(target_names))
 			overlap_targets.sort()
 			
 		else:
 			overlap_targets = []
 				
 		return overlap_targets
 		
 	def find_targets_tissue(self,tf_name, tissue_name):
 		'''
		Return Targets regulated by the tf in a given tissue
		'''
		
 		#query
 		if self.tfdb is not None:
 			t = (tf_name, tissue_name)
 			res = self.tfdb.execute("SELECT Target FROM Target2TF2Tissue "
 									"WHERE TF = ? and Tissue = ? ORDER BY Target", t).fetchall()
 			target_names = [r[0] for r in res]	
 		else:
 			target_names = []
 				
 		return target_names
 		
 	
		 							 

#test functions
if __name__ == "__main__":
	a=TFTA()
	'''
	ss2 = a.Is_tf_target('FOS','ELF3');
	if ss2 is not None:
		print 'ss2='+str(ss2)
	else:
		print 'ss2 is none'
	'''
	ss3 = a.find_tfs(['ELF3','KRAS']) 
	if ss3 is not None:
		print 'lenth(ss3)='+str(len(ss3))
		print 'ss3='+','.join(ss3)
		
	ss4 = a.find_targets(['FOS','STAT3'])
	if ss4 is not None:
		print 'ss4='+','.join(ss4)
		
	else:
		print 'ss4 is none'
	'''
	ss6=a.Is_tf_target_tissue('BSAP','USP1','bladder')
	if ss6 is not None:
		print 'ss6='+str(ss6)
	else:
		print 'ss6 is none'
		
	ss7=a.find_tfs_tissue('USP1','bladder')
	if ss7 is not None:
		print 'lenth(ss7)='+str(len(ss7))
		print 'ss7='+','.join(ss7)
		
	ss8=a.find_targets_tissue('BSAP','bladder')
	if ss8 is not None:
		print 'ss8='+','.join(ss8)
		
	else:
		print 'ss8 is none'
 	
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
 	'''
 	'''
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
 	
 	targetlist=a.find_targets_tfs_genes(['IRF3','MAX'],targets)
 	print 'targetlist='+','.join(targetlist)
 	print 'length(targetlist)='+str(len(targetlist))
 	'''
 	
 	