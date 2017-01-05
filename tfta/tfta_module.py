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
		self.tasks = ['IS-TF-TARGET', 'FIND-TF-TARGET', 'FIND-TARGET-TF']
		
		#Send subscribe messages
		for task in self.tasks:
			msg_txt = '(subscrbe :content (request &key :content (%s . *)))' % task
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
		print 'tf='+tf.name
		
		target_arg = content_list.get_keyword_arg(':target')
		target = self._get_target(target_arg)
		target_name = target.name
		
		reply_content = KQMLList()
		try:
			is_target = self.tfta.is_tf_target(tf_name, target_name)
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
		
	def _get_target(self, target_arg):
		target_str = str(target_arg)
		target_str = self.decode_description('<ekb>' + target_str + '</ekb>')
		tp = TripsProcessor(target_str)
		terms = tp.tree.findall('TERM')
		term_id = terms[0].attrib['id']
		agent = tp._get_agent_by_id(term_id, None)
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
	
	
		
		



