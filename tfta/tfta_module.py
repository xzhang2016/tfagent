'''
TFTA module is to receive and decode messages and send responses from and to other agents in the system
'''

import sys
import logging
from kqml import KQMLModule, KQMLPerformative, KQMLList
from tfta import TFTA, TFNotFoundException, TargetNotFoundException
from indra.trips.processor import TripsProcessor

logger = logging.getLogger('TFTA')

class TFTA_Module(KQMLModule):
	'''
	TFTA module is used to receive and decode messages and send responses from and to other agents in the system
    '''
	def __init__(self, argv):
		super(TFTA_Module, self).__init__(argv)
		self.tasks = ['IS-TF-TARGET', 'FIND-TF-TARGET', 'FIND-TARGET-TF',
		              'FIND_PATHWAY_GENE', 'FIND_PATHWAY_DB_GENE', 'FIND_TF_PATHWAY',
		              'FIND_GENE_PATHWAY', 'FIND_PATHWAY_CHEMICAL', 'FIND_TF_CHEMICAL',
		              'FIND_COMMON_TF_GENES', 'FIND_OVERLAP_TARGETS_TF_GENES',
		              'IS-TF-TARGET-TISSUE', 'FIND-TF-TARGET-TISSUE', 'FIND-TARGET-TF-TISSUE']
		
		#Send subscribe messages
		for task in self.tasks:
			msg_txt = '(subscribe :content (request &key :content (%s . *)))' % task
			self.send(KQMLPerformative.from_string(msg_txt))
		
		#Instantiate a singleton TFTA agent
		self.tfta = TFTA()
		#send ready message
		self.ready()
		super(TFTA_Module, self).start()
		
	def receive_request(self, msg, content):
		'''
		If a "request" message is received, decode the task and the content
        and call the appropriate function to prepare the response. A reply
        message is then sent back.
        '''
		content_list = content
		task_str = content_list[0].to_string().upper()
		if task_str == 'IS-TF-TARGET':
			reply_content = self.respond_is_tf_target(content_list)
		elif task_str == 'FIND-TF-TARGET':
			reply_content = self.respond_find_tf_targets(content_list)
		elif task_str == 'FIND-TARGET-TF':
			reply_content = self.respond_find_target_tfs(content_list)
		elif task_str == 'FIND_PATHWAY_GENE':
			reply_content = self.respond_find_pathway_gene(content_list)
		elif task_str == 'FIND_PATHWAY_DB_GENE':
			reply_content = self.respond_find_pathway_db_gene(content_list)
		elif task_str == 'FIND_TF_PATHWAY':
			reply_content = self.respond_find_tf_pathway(content_list)
		elif task_str == 'FIND_GENE_PATHWAY':
			reply_content = self.respond_find_gene_pathway(content_list)
		elif task_str == 'FIND_PATHWAY_CHEMICAL':
			reply_content = self.respond_find_pathway_chemical(content_list)
		elif task_str == 'FIND_TF_CHEMICAL':
			reply_content = self.respond_find_tf_chemical(content_list)
		elif task_str == 'FIND_COMMON_TF_GENES':
			reply_content = self.respond_find_common_tfs_genes(content_list)
		elif task_str == 'FIND_OVERLAP_TARGETS_TF_GENES':
			reply_content = self.respond_find_overlap_targets_tf_genes(content_list)
		elif task_str == 'IS-TF-TARGET-TISSUE':
			reply_content = self.respond_is_tf_target_tissue(content_list)
		elif task_str == 'FIND-TF-TARGET_TISSUE':
			reply_content = self.respond_find_tf_targets_tissue(content_list)
		elif task_str == 'FIND-TARGET-TF_TISSUE':
			reply_content = self.respond_find_target_tfs_tissue(content_list)
		else:
			self.error_reply(msg, 'unknown request task ' + task_str)
			return
        	
		reply_msg = KQMLPerformative('reply')
		reply_msg.set_parameter(':content', reply_content)
		self.reply(msg, reply_msg)
	
	'''	
	def respond_dont_know(self, msg, content_string):
		#TODO: need to rewrite this function
		resp = '(ONT::TELL :content (ONT::DONT-KNOW :content %s))' %\
            content_string
        resp_list = KQMLList.from_string(resp)
        reply_msg = KQMLPerformative('reply')
        reply_msg.set_parameter(':content', resp_list)
		self.reply(msg, reply_msg)
	'''
		
	def respond_is_tf_target(self, content_list):
		'''
		Response content to is-tf-target request
		'''
		tf_arg = content_list.get_keyword_arg(':tf')
		tf = self._get_target(tf_arg)
		tf_name = tf.name
		#print 'tf='+tf.name
		
		target_arg = content_list.get_keyword_arg(':target')
		target = self._get_target(target_arg)
		target_name = target.name
		
		reply_content = KQMLList()
		try:
			is_target = self.tfta.Is_tf_target(tf_name, target_name)
		except TFNotFoundException:
			reply_content.add('FAILURE :reason TF_NOT_FOUND')
			return reply_content
			
		status = 'SUCCESS'
		if is_target:
			is_target_str = 'TRUE'
		else:
			is_target_str = 'FALSE'
		msg_str = '%s :is-tf-target %s' %\
                  (status, is_target_str)
		reply_content.add(msg_str)
		return reply_content
	
	
	tissue_list = ['bladder','blood','bone','bone_marrow','brain','cervix','colon','eye',
	               'heart','kidney','larynx','liver','lung','lymph_node','mammary_gland',
	               'muscle','ovary','pancreas','peripheral_nervous_system','placenta',
	               'prostate','skin','small_intestine','soft_tissue','spleen','stomach',
	               'testis','thymus','tongue','uterus']
	               
	
	def respond_is_tf_target_tissue(self,content_list):
		'''
		Response content to is-tf-target-tissue request
		'''
		tf_arg = content_list.get_keyword_arg(':tf')
		tf = self._get_target(tf_arg)
		tf_name = tf.name
		print 'tf='+tf.name
		
		target_arg = content_list.get_keyword_arg(':target')
		target = self._get_target(target_arg)
		target_name = target.name
		
		tissue_arg = content_list.get_keyword_arg(':tissue')
		tissue_name = tissue_arg[0].to_string()
		
			
		reply_content = KQMLList()
		
		if tissue_name not in tissue_list:
			reply_content.add('FAILURE :reason INVALID_TISSUE')
			return reply_content
		
		try:
			is_target = self.tfta.Is_tf_target_tissue(tf_name, target_name, tissue_name)
		except TFNotFoundException:
			reply_content.add('FAILURE :reason TF_NOT_FOUND')
			return reply_content
			
		status = 'SUCCESS'
		if is_target:
			is_target_str = 'TRUE'
		else:
			is_target_str = 'FALSE'
		msg_str = '%s :is-tf-target %s' %\
                  (status, is_target_str)
		reply_content.add(msg_str)
		return reply_content
		
	def respond_find_tf_targets(self, content_list):
		'''
		Response content to find-tf-target request
		For a tf, reply with the targets found
		'''
		tf_arg = content_list.get_keyword_arg(':tf')
		tf = self._get_target(tf_arg)
		tf_name = tf.name
		target_names, targetEntrezIDs = self.tfta.find_targets(tf_name)
		
		target_list_str = ''
		for tg, ez in zip(target_names,targetEntrezIDs):
			target_list_str += '(:name %s :EntrezID %s) ' % (tg, ez)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :targets (' + target_list_str + '))')
		return reply_content
		
	def respond_find_tf_targets_tissue(self, content_list):
		'''
		Response content to find-tf-target request
		For a tf, reply with the targets found within given tissue
		'''
		tf_arg = content_list.get_keyword_arg(':tf')
		tf = self._get_target(tf_arg)
		tf_name = tf.name
		
		tissue_arg = content_list.get_keyword_arg(':tissue')
		tissue_name = tissue_arg[0].to_string()
		
		if tissue_name not in tissue_list:
			reply_content = KQMLList()
			reply_content.add('FAILURE :reason INVALID_TISSUE')
			return reply_content
		
		target_names = self.tfta.find_targets_tissue(tf_name, tissue_name)
		
		target_list_str = ''
		for tg in target_names:
			target_list_str += '(:name %s) ' % tg.encode('ascii', 'ignore')
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :targets (' + target_list_str + '))')
		return reply_content
		
	def respond_find_target_tfs(self, content_list):
		'''
		Response content to find-target-tf request
		For a target, reply the tfs found
		'''
		target_arg = content_list.get_keyword_arg(':target')
		target = self._get_target(target_arg)
		target_name = target.name
		tf_names = self.tfta.find_tfs(target_name)
		tf_list_str = ''
		for tf in tf_names:
			tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')
		reply_content = KQMLList.from_string(
            '(SUCCESS :tfs (' + tf_list_str + '))')
		return reply_content
		
	def respond_find_target_tfs_tissue(self, content_list):
		'''
		Response content to find-target-tf request
		For a target, reply the tfs found within a given tissue
		'''
		target_arg = content_list.get_keyword_arg(':target')
		target = self._get_target(target_arg)
		target_name = target.name
		
		tissue_arg = content_list.get_keyword_arg(':tissue')
		tissue_name = tissue_arg[0].to_string()
		
		if tissue_name not in tissue_list:
			reply_content = KQMLList()
			reply_content.add('FAILURE :reason INVALID_TISSUE')
			return reply_content
		
		tf_names = self.tfta.find_tfs_tissue(target_name, tissue_name)
		tf_list_str = ''
		for tf in tf_names:
			tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :tfs (' + tf_list_str + '))')
		return reply_content
		
	def respond_find_pathway_gene(self,content_list):
		'''
		Response content to find_pathway_gene request
		For a given gene list, reply the related pathways information
		'''
		gene_arg = content_list.get_keyword_arg(':gene')
		genes = self._get_targets(gene_arg)
		gene_names = []
		for gene in genes:
			gene_names.append(gene.name)
		
		reply_content = KQMLList()
		try:
			pathwayId,pathwayName,externalId,source,dblink = self.tfta.find_pathways_from_genelist(gene_names)
		except PathwayNotFoundException:
			reply_content.add('FAILURE :reason PathwayNotFoundException')
			return reply_content
		
		pathway_list_str = ''
		for pn, eid, src, dbl in zip(pathwayName,externalId,source,dblink):
			pathway_list_str += '(:name %s :externalId %s :source %s :dblink %s) ' % (pn, eid ,src, dbl)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
		return reply_content
		
	def respond_find_pathway_db_gene(self, content_list):
		'''
		Response content to FIND_PATHWAY_DB_GENE request
		For a given gene list and certain db source, reply the related pathways information
		'''
		db_arg = content_list.get_keyword_arg(':database')
		db_name = db_arg[0].to_string()
		
		gene_arg = content_list.get_keyword_arg(':gene')
		genes = self._get_targets(gene_arg)
		gene_names=[]
		for gene in genes:
			gene_names.append(gene.name)
		
		reply_content = KQMLList()
		try:
			pathwayId,pathwayName,externalId,source,dblink = self.tfta.find_pathways_from_dbsource_geneName(dbsource,gene_name)
		except PathwayNotFoundException:
			reply_content.add('FAILURE :reason PathwayNotFoundException')
			return reply_content
		
		pathway_list_str = ''
		for pn, eid, src, dbl in zip(pathwayName,externalId,source,dblink):
			pathway_list_str += '(:name %s :externalId %s :source %s :dblink %s) ' % (pn, eid ,src, dbl)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
		return reply_content
		
	def respond_find_tf_pathway(self, content_list):
		'''
		Response content to FIND_TF_PATHWAY request
		For a given pathway name, reply the tfs within the pathway
		'''
		pathway_arg = content_list.get_keyword_arg(':pathway')
		pathway_name = pathway_arg[0].to_string()
		
		reply_content = KQMLList()
		try:
			pathwayId,pathwayName,tflist = self.tfta.find_tfs_from_pathwayName(pathway_name)
		except PathwayNotFoundException:
			reply_content.add('FAILURE :reason PathwayNotFoundException')
			return reply_content
		
		pathway_list_str = ''
		for pid, pn in zip(pathwayId, pathwayName):
			tf_list_str = ''
			for tf in tflist[pid]:
				tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')
				
			tf_list_str = '(' + tf_list_str + ')'
			pathway_list_str += '(:name %s :tfs %s) ' % (pn, tf_list_str)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
		return reply_content
		
	def respond_find_gene_pathway(self, content_list):
		'''
		Response content to FIND_GENE_PATHWAY request
		For a given pathway name, reply the genes within the pathway
		'''
		pathway_arg = content_list.get_keyword_arg(':pathway')
		pathway_name = pathway_arg[0].to_string()
		
		reply_content = KQMLList()
		try:
			pathwayId,pathwayName,genelist = self.tfta.find_genes_from_pathwayName(pathway_name)
		except PathwayNotFoundException:
			reply_content.add('FAILURE :reason PathwayNotFoundException')
			return reply_content
		
		pathway_list_str = ''
		for pid, pn in zip(pathwayId, pathwayName):
			gene_list_str = ''
			for gene in genelist[pid]:
				gene_list_str += '(:name %s) ' % gene.encode('ascii', 'ignore')
				
			gene_list_str = '(' + gene_list_str + ')'
			pathway_list_str += '(:name %s :genes %s) ' % (pn, gene_list_str)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
		return reply_content
		
	def respond_find_pathway_chemical(self, content_list):
		'''
		Response content to FIND_PATHWAY_CHEMICAL request
		For a given chemical name, reply the pathways involving the chemical
		'''
		chemical_arg = content_list.get_keyword_arg(':chemical')
		chemical_name = chemical_arg[0].to_string()
		
		reply_content = KQMLList()
		try:
			pathwayId,pathwayName,externalId,source,dblink = self.tfta.find_pathways_from_chemical(chemical_name)
		except PathwayNotFoundException:
			reply_content.add('FAILURE :reason PathwayNotFoundException')
			return reply_content
		
		pathway_list_str = ''
		for pn, eid, src, dbl in zip(pathwayName,externalId,source,dblink):
			pathway_list_str += '(:name %s :externalId %s :source %s :dblink %s) ' % (pn, eid ,src, dbl)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
		return reply_content
		
	def respond_find_tf_chemical(self, content_list):
		'''
		Response content to FIND_TF_CHEMICAL request
		For a given chemical name, reply the tfs within the pathways involving the chemical
		'''
		chemical_arg = content_list.get_keyword_arg(':chemical')
		chemical_name = chemical_arg[0].to_string()
		
		reply_content = KQMLList()
		try:
			pathwayId,pathwayName,tflist = self.tfta.find_tfs_from_pathwaysWithChemical(chemical_name)
		except PathwayNotFoundException:
			reply_content.add('FAILURE :reason PathwayNotFoundException')
			return reply_content
			
		pathway_list_str = ''
		for pid, pn in zip(pathwayId, pathwayName):
			tf_list_str = ''
			for tf in tflist[pid]:
				tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')
			
			tf_list_str = '(' + tf_list_str + ')'
			pathway_list_str += '(:name %s :tfs %s) ' % (pn, tf_list_str)
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
		return reply_content
		
	def respond_find_common_tfs_genes(self, content_list):
		'''
		Response content to FIND_COMMON_TF_GENES request
		For a given target list, reply the tfs regulating these genes
		and the frequency of each TF
		'''
		target_arg = content_list.get_keyword_arg(':target')
		targets = self._get_targets(target_arg)
		target_names = []
		for target in targets:
			target_names.append(target.name)
			
		tf_names, tf_counts = self.tfta.find_tfs_count(target_names)
		tf_list_str = ''
		for tf,ct in zip(tf_names, tf_counts):
			tf_list_str += '(:name %s :count %s) ' % (tf, ct)
		reply_content = KQMLList.from_string(
            '(SUCCESS :tfs (' + tf_list_str + '))')
		return reply_content
		
	def respond_find_overlap_targets_tf_genes(self, content_list):
		'''
		For given tf list, find the targets which are regulated by all of them;
		then return the overlap between the targets and the given gene list
		'''
		tf_arg = content_list.get_keyword_arg(':tf')
		tfs = self._get_targets(tf_arg)
		tf_names = []
		for tf in tfs:
			tf_names.append(tf.name)
			
		target_arg = content_list.get_keyword_arg(':target')
		targets = self._get_targets(target_arg)
		target_names = []
		for tg in targets:
			target_names.append(tg.name)
			
		overlap_targets = self.tfta.find_overlap_targets_tfs_genes(tf_names,target_names)
		
		target_list_str = ''
		for tg in overlap_targets:
			target_list_str += '(:name %s ) ' % tg.encode('ascii', 'ignore')
			
		reply_content = KQMLList.from_string(
            '(SUCCESS :targets (' + target_list_str + '))')
		return reply_content
		
	def _get_target(self, target_arg):
		target_str = str(target_arg)
		target_str = self.decode_description('<ekb>' + target_str + '</ekb>')
		tp = TripsProcessor(target_str)
		terms = tp.tree.findall('TERM')
		term_id = terms[0].attrib['id']
		agent = tp._get_agent_by_id(term_id, None)
		return agent
		
	def _get_targets(self, target_arg):
		target_str = str(target_arg)
		target_str = self.decode_description('<ekb>' + target_str + '</ekb>')
		tp = TripsProcessor(target_str)
		terms = tp.tree.findall('TERM')
		agent = []
		for term in terms:
			term_id = term[0].attrib['id']
			agent.append(tp._get_agent_by_id(term_id, None))
			
		return agent
		
	@staticmethod
	def decode_description(descr):
		if descr[0] == '"':
			descr = descr[1:]
		if descr[-1] == '"':
			descr = descr[:-1]
		descr = descr.replace('\\"', '"')
		return descr
		
if __name__ == "__main__":
	TFTA_Module(['-name', 'TFTA'] + sys.argv[1:])
	
	
		
		



