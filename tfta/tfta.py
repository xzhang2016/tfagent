#TFTA stands for TF-target agent whose task is to search for targets known to a TF and vice versa

import MySQLdb as mdb
import sys
import logging
 
logger = logging.getLogger('TFTA')
 
socket_path='/Applications/XAMPP/xamppfiles/var/mysql/mysql.sock';
 
class TFNotFoundException(Exception):
	def __init__(self, *args, **kwargs):
		Exception.__init__(self, *args, **kwargs)
 
class TargetNotFoundException(Exception):
	def __init__(self, *args, **kwargs):
		Exception.__init__(self, *args, **kwargs)
 	
class TFTA:
	def __init__(self):
 		#connect the database of tf-target:tfdb
 		
		self.tfdb=None;
		try:
			self.tfdb = mdb.connect(unix_socket=socket_path, host='localhost', user='user01', passwd='asd2016*', db='tfdb');
 			
		except mdb.Error,e:
			print "Error %d: %s" % (e.args[0],e.args[1])
			logger.error('TFTA could not load tf-target database.')
			sys.exit(1)
 			
	def __del__(self):
		self.tfdb.close()
 		
	def aliasToOfficial(self,gene):
 		'''
 		Search table symbol2alias, convert gene to its official gene symbol
 		'''
 		#query = 'SELECT HGNCSymbol FROM symbol2alias where HGNCSymbol= %s;' % gene
 		
		with self.tfdb:
			cur = self.tfdb.cursor()
			cur.execute("SELECT HGNCSymbol FROM symbol2alias where HGNCSymbol = %s OR alias = %s;", (gene, gene))
 			
 			numrows = int(cur.rowcount)
 			if numrows > 0:
				row = cur.fetchone()
				return row[0]
			else:
				return None
 				
 		
	def Is_tf_target(self,tf,target):
 		'''
 		Return True if the tf regulates the target, and False if not
 		'''
 		tf_official=tf
 		target_official=target
 		
 		#query the database
 		with self.tfdb:
 			cur = self.tfdb.cursor()
 			#convert tf and target to official symbols
 			cur.execute("SELECT HGNCSymbol FROM symbol2alias where HGNCSymbol = %s OR alias = %s;", (target, target))
 			numrows = int(cur.rowcount)
 			if numrows > 0:
				target_official = cur.fetchone()
				
			cur.execute("SELECT HGNCSymbol FROM symbol2alias where HGNCSymbol = %s OR alias = %s;", (tf, tf))
 			numrows = int(cur.rowcount)
 			if numrows > 0:
				tf_official = cur.fetchone()
 			
 			cur.execute("SELECT * FROM CombinedDB where TF = %s AND Target = %s;", (tf_official, target_official))
 			numrows = int(cur.rowcount)
 			if numrows>0 :
 				return True
 			else:
 				return False
 		

	def find_tfs(self,target):
		'''
		Return TFs regulating the target
		'''
 			
 		#query
 		tfs = []
 		target_official=target
 		
 		with self.tfdb:
 			cur = self.tfdb.cursor()
 			#convert to official symbol
 			cur.execute("SELECT HGNCSymbol FROM symbol2alias where HGNCSymbol = %s OR alias = %s;", (target, target))
 			numrows = int(cur.rowcount)
 			if numrows > 0:
				target_official = cur.fetchone()
 			
 			cur.execute("SELECT TF FROM CombinedDB where Target = %s;", (target_official))
 			numrows = int(cur.rowcount)
 			if numrows > 0:
 				rows = cur.fetchall()
 				tfs = [row[0] for row in rows]
 				
 		return tfs
 		
 	def find_targets(self,tf):
 		'''
		Return Targets regulated by the tf
		'''
		
 		#query
 		targets = []
 		tf_official = tf
 		
 		with self.tfdb:
 			cur = self.tfdb.cursor()
 			#convert the tf to its official symbol
 			cur.execute("SELECT HGNCSymbol FROM symbol2alias where HGNCSymbol = %s OR alias = %s;", (tf, tf))
 			numrows = int(cur.rowcount)
 			if numrows > 0:
				tf_official = cur.fetchone()
 			
 			cur.execute("SELECT Target FROM CombinedDB where TF = %s;", (tf_official))
 			numrows = int(cur.rowcount)
 			if numrows > 0:
 				rows = cur.fetchall()
 				targets = [row[0] for row in rows]
 				
 		return targets
		 					

 		 

#test functions
if __name__ == "__main__":
	a=TFTA()
	ss = a.aliasToOfficial('JTK7')
	if ss is not None:
		print 'ss='+ss
	else:
		print 'ss is none'
	
	ss1 = a.aliasToOfficial('ecs')
	if ss1 is None:
		print 'No official symbol exists according to our database for your input '+'ecs'
	else:
		print 'ss1='+ss1
		
	ss2 = a.Is_tf_target('abl','bax');
	if ss2 is not None:
		print 'ss2='+str(ss2)
	else:
		print 'ss2 is none'
	
	ss3 = a.find_tfs('bcl2l4') #bax
	if ss3 is not None:
		print 'lenth(ss3)='+str(len(ss3))
		print 'ss3='+','.join(ss3)
		
		ss33=set(ss3)
		print 'length(ss33)='+str(len(ss33))
		print 'ss33='+';'.join(ss33)
	else:
		print 'ss3 is none'
		
	ss4 = a.find_targets('abl1')
	if ss4 is not None:
		print 'ss4='+','.join(ss4)
	else:
		print 'ss4 is none'
 			