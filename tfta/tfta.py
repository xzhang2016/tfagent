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
 			res = self.tfdb.execute("SELECT Target FROM Target2TF2Tissue WHERE TF = ? AND Tissue = ?", t).fetchall()
 			if not res:
 				raise TFNotFoundException
 			for r in res:
 				if r[0].upper() == target_name.upper():
 					return True
 		
 		return False
 		

	def find_tfs(self,target_name):
		'''
		Return TFs regulating the target
		'''
 			
 		#query
 		if self.tfdb is not None:
 			t = (target_name,)
 			res = self.tfdb.execute("SELECT TF FROM CombinedDB "
 									"WHERE Target = ?  ORDER BY TF", t).fetchall()
 			tf_names = [r[0] for r in res]
 		else:
 			tf_names = []
 				
 		return tf_names
 		
 	def find_tfs_tissue(self,target_name,tissue_name):
		'''
		Return TFs regulating the target in a given tissue
		'''
 			
 		#query
 		if self.tfdb is not None:
 			t = (target_name,tissue_name)
 			res = self.tfdb.execute("SELECT TF FROM Target2TF2Tissue "
 									"WHERE Target = ? and Tissue = ? ORDER BY TF", t).fetchall()
 			tf_names = [r[0] for r in res]
 		else:
 			tf_names = []
 				
 		return tf_names
 		
 	def find_targets(self,tf_name):
 		'''
		Return Targets regulated by the tf
		'''
		
 		#query
 		if self.tfdb is not None:
 			t = (tf_name,)
 			res = self.tfdb.execute("SELECT Target,TargetEntrezID FROM CombinedDB "
 									"WHERE TF = ?  ORDER BY Target", t).fetchall()
 			target_names = [r[0] for r in res]
 			targetEntrez = [r[1] for r in res]
 		else:
 			target_names = []
 			targetEntrez = []
 				
 		return target_names, targetEntrez
 		
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
	
	ss2 = a.Is_tf_target('FOS','ELF3');
	if ss2 is not None:
		print 'ss2='+str(ss2)
	else:
		print 'ss2 is none'
	
	ss3 = a.find_tfs('ELF3') 
	if ss3 is not None:
		print 'lenth(ss3)='+str(len(ss3))
		print 'ss3='+','.join(ss3)
		
		ss33=set(ss3)
		print 'length(ss33)='+str(len(ss33))
		print 'ss33='+';'.join(ss33)
	else:
		print 'ss3 is none'
		
	ss4,ss5 = a.find_targets('FOS')
	if ss4 is not None:
		print 'ss4='+','.join(ss4)
		
	else:
		print 'ss4 is none'
	
	ss6=a.Is_tf_target_tissue('BSAP','USP1','bladder')
	if ss6 is not None:
		print 'ss6='+str(ss6)
	else:
		print 'ss6 is none'
		
	ss7=a.find_tfs_tissue('USP1','bladder')
	if ss7 is not None:
		print 'lenth(ss7)='+str(len(ss7))
		print 'ss7='+','.join(ss7)
		
		ss77=set(ss7)
		print 'length(ss77)='+str(len(ss77))
		print 'ss77='+';'.join(ss77)
	else:
		print 'ss7 is none'
		
	ss8=a.find_targets_tissue('BSAP','bladder')
	if ss8 is not None:
		print 'ss8='+','.join(ss8)
		
	else:
		print 'ss8 is none'
 			