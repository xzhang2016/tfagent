"""TFTA module is to receive and decode messages and send responses from
and to other agents in the system"""

import sys
import logging
from kqml import KQMLModule, KQMLPerformative, KQMLList
from tfta import TFTA, TFNotFoundException, TargetNotFoundException, PathwayNotFoundException
from indra.trips.processor import TripsProcessor

logger = logging.getLogger('TFTA')

class TFTA_Module(KQMLModule):
    """TFTA module is used to receive and decode messages and send
    responses from and to other agents in the system."""
    def __init__(self, argv):
        super(TFTA_Module, self).__init__(argv)
        self.tasks = ['IS-TF-TARGET', 'FIND-TF-TARGET', 'FIND-TARGET-TF',
                      'FIND_PATHWAY_GENE', 'FIND_PATHWAY_DB_GENE',
                      'FIND_TF_PATHWAY', 'FIND_GENE_PATHWAY',
                      'FIND_PATHWAY_KEYWORD', 'FIND_TF_KEYWORD',
                      'FIND_COMMON_TF_GENES', 'FIND_OVERLAP_TARGETS_TF_GENES',
                      'FIND_COMMON_PATHWAY_GENES',
                      'IS-TF-TARGET-TISSUE', 'FIND-TF-TARGET-TISSUE',
                      'FIND-TARGET-TF-TISSUE']

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
        """If a "request" message is received, decode the task and
        the content and call the appropriate function to prepare the
        response. A reply message is then sent back.
        """
        task_str = content.head()
        if task_str == 'IS-TF-TARGET':
            reply_content = self.respond_is_tf_target(content)
        elif task_str == 'FIND-TF-TARGET':
            reply_content = self.respond_find_tf_targets(content)
        elif task_str == 'FIND-TARGET-TF':
            reply_content = self.respond_find_target_tfs(content)
        elif task_str == 'FIND_PATHWAY_GENE':
            reply_content = self.respond_find_pathway_gene(content)
        elif task_str == 'FIND_PATHWAY_DB_GENE':
            reply_content = self.respond_find_pathway_db_gene(content)
        elif task_str == 'FIND_TF_PATHWAY':
            reply_content = self.respond_find_tf_pathway(content)
        elif task_str == 'FIND_GENE_PATHWAY':
            reply_content = self.respond_find_gene_pathway(content)
        elif task_str == 'FIND_PATHWAY_KEYWORD':
            reply_content = self.respond_find_pathway_chemical(content)
        elif task_str == 'FIND_TF_KEYWORD':
            reply_content = self.respond_find_tf_chemical(content)
        elif task_str == 'FIND_COMMON_TF_GENES':
            reply_content = self.respond_find_common_tfs_genes(content)
        elif task_str == 'FIND_OVERLAP_TARGETS_TF_GENES':
            reply_content = self.respond_find_overlap_targets_tf_genes(content)
        elif task_str == 'FIND_COMMON_PATHWAY_GENES':
			reply_content = self.respond_find_common_pathway_genes(content_list)
        elif task_str == 'IS-TF-TARGET-TISSUE':
            reply_content = self.respond_is_tf_target_tissue(content)
        elif task_str == 'FIND-TF-TARGET_TISSUE':
            reply_content = self.respond_find_tf_targets_tissue(content)
        elif task_str == 'FIND-TARGET-TF_TISSUE':
            reply_content = self.respond_find_target_tfs_tissue(content)
        else:
            self.error_reply(msg, 'unknown request task ' + task_str)
            return

        reply_msg = KQMLPerformative('reply')
        reply_msg.set('content', reply_content)
        self.reply(msg, reply_msg)

    def respond_is_tf_target(self, content):
        """Response content to is-tf-target request."""
        tf_arg = content.gets('tf')
        tf = _get_target(tf_arg)
        tf_name = tf.name

        target_arg = content.gets('target')
        target = _get_target(target_arg)
        target_name = target.name

        try:
            is_target = self.tfta.Is_tf_target(tf_name, target_name)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply

        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if is_target else 'FALSE'
        reply.set('is-tf-target', is_target_str)
        return reply

    tissue_list = ['bladder','blood','bone','bone_marrow','brain','cervix',
                   'colon','eye','heart','kidney','larynx','liver','lung',
                   'lymph_node','mammary_gland','muscle','ovary','pancreas',
                   'peripheral_nervous_system','placenta','prostate','skin',
                   'small_intestine','soft_tissue','spleen','stomach',
                   'testis','thymus','tongue','uterus']


    def respond_is_tf_target_tissue(self, content):
        """Response content to is-tf-target-tissue request."""
        tf_arg = content.gets('tf')
        tf = _get_target(tf_arg)
        tf_name = tf.name
        print 'tf=' + tf.name

        target_arg = content.gets('target')
        target = _get_target(target_arg)
        target_name = target.name

        tissue_arg = content.get('tissue')
        tissue_name = tissue_arg.head()

        if tissue_name not in tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply

        try:
            is_target = self.tfta.Is_tf_target_tissue(tf_name, target_name, tissue_name)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply

        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if is_target else 'FALSE'
        reply.set('is-tf-target', is_target_str)
        return reply

    def respond_find_tf_targets(self, content):
        """Response content to find-tf-target request
        For a tf list, reply with the targets found
        """
        tf_arg = content.gets('tf')
        tfs = _get_targets(tf_arg)
        tf_names = []
        for tf in tfs:
            tf_names.append(tf.name)
            
        try:
            target_names = self.tfta.find_targets(tf_names)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        
        if len(target_names):
            target_list_str = ''
            for tg in target_names:
                target_list_str += '(:name %s ) ' % tg.encode('ascii', 'ignore')

            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        return reply

    def respond_find_tf_targets_tissue(self, content):
        """Response content to find-tf-target request
        For tf list, reply with the targets found within given tissue
        """
        tf_arg = content.gets('tf')
        tfs = _get_targets(tf_arg)
        tf_names = []
        for tf in tfs:
            tf_names.append(tf.name)

        tissue_arg = content.get('tissue')
        tissue_name = tissue_arg.head()

        if tissue_name not in tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        
        try:
            target_names = self.tfta.find_targets_tissue(tf_names, tissue_name)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        
        if len(target_names):
            target_list_str = ''
            for tg in target_names:
                target_list_str += '(:name %s) ' % tg.encode('ascii', 'ignore')

            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        return reply

    def respond_find_target_tfs(self, content):
        """Response content to find-target-tf request
        For a target list, reply the tfs found"""
        target_arg = content.gets('target')
        targets = _get_targets(target_arg)
        target_names = []
        for target in targets:
            target_names.append(target.name)
        try:
            tf_names = self.tfta.find_tfs(target_names)
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply
            
        if len(tf_names):    
            tf_list_str = ''
            for tf in tf_names:
                tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')
            reply = KQMLList.from_string(
                '(SUCCESS :tfs (' + tf_list_str + '))')
        else:
            reply = make_failure('NO_TF_FOUND')
        return reply

    def respond_find_target_tfs_tissue(self, content):
        """Response content to find-target-tf request
        For a target, reply the tfs found within a given tissue"""
        target_arg = content.gets('target')
        target = _get_target(target_arg)
        target_name = target.name

        tissue_arg = content.get('tissue')
        tissue_name = tissue_arg.head()

        if tissue_name not in tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        
        try:
            tf_names = self.tfta.find_tfs_tissue(target_name, tissue_name)
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply
            
        if len(tf_names):
            tf_list_str = ''
            for tf in tf_names:
                tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')

            reply = KQMLList.from_string(
                '(SUCCESS :tfs (' + tf_list_str + '))')
        else:
            reply = make_failure('NO_TF_FOUND')
        return reply

    def respond_find_pathway_gene(self,content):
        """Response content to find_pathway_gene request
        For a given gene list, reply the related pathways information"""
        gene_arg = content.gets('gene')
        genes = _get_targets(gene_arg)
        gene_names = []
        for gene in genes:
            gene_names.append(gene.name)

        try:
            pathwayId, pathwayName, externalId, source,dblink = \
                self.tfta.find_pathways_from_genelist(gene_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName, externalId, source, dblink):
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pn, eid ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_db_gene(self, content):
        """Response content to FIND_PATHWAY_DB_GENE request
        For a given gene list and certain db source, reply the related
        pathways information"""
        db_arg = content.get('database')
        db_name = db_arg.head()

        gene_arg = content.gets('gene')
        genes = _get_targets(gene_arg)
        gene_names = []
        for gene in genes:
            gene_names.append(gene.name)

        try:
            pathwayId, pathwayName, externalId, source, dblink = \
                self.tfta.find_pathways_from_dbsource_geneName(dbsource,gene_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName, externalId, source, dblink):
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pn, eid ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tf_pathway(self, content):
        """Response content to FIND_TF_PATHWAY request
        For a given pathway name, reply the tfs within the pathway"""
        pathway_arg = content.get('pathway')
        pathway_name = pathway_arg.head()

        try:
            pathwayId,pathwayName,tflist = \
                self.tfta.find_tfs_from_pathwayName(pathway_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pid, pn in zip(pathwayId, pathwayName):
            tf_list_str = ''
            for tf in tflist[pid]:
                tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')

            tf_list_str = '(' + tf_list_str + ')'
            pathway_list_str += '(:name %s :tfs %s) ' % (pn, tf_list_str)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_gene_pathway(self, content):
        """Response content to FIND_GENE_PATHWAY request
        For a given pathway name, reply the genes within the pathway"""
        pathway_arg = content.get('pathway')
        pathway_name = pathway_arg.head()

        try:
            pathwayId,pathwayName,genelist = \
                self.tfta.find_genes_from_pathwayName(pathway_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pid, pn in zip(pathwayId, pathwayName):
            gene_list_str = ''
            for gene in genelist[pid]:
                gene_list_str += '(:name %s) ' % gene.encode('ascii', 'ignore')

            gene_list_str = '(' + gene_list_str + ')'
            pathway_list_str += '(:name %s :genes %s) ' % (pn, gene_list_str)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_chemical(self, content):
        """Response content to FIND_PATHWAY_KEYWORD request
        For a given chemical name, reply the pathways involving the chemical"""
        chemical_arg = content.get('keyword')
        chemical_name = chemical_arg.head()

        try:
            pathwayId, pathwayName, externalId, source, dblink = \
                self.tfta.find_pathways_from_chemical(chemical_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName,externalId,source,dblink):
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pn, eid ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tf_chemical(self, content):
        """Response content to FIND_TF_KEYWORD request
        For a given chemical name, reply the tfs within the pathways
        involving the chemical"""
        chemical_arg = content.get('keyword')
        chemical_name = chemical_arg.head()

        try:
            pathwayId, pathwayName, tflist = \
                self.tfta.find_tfs_from_pathwaysWithChemical(chemical_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pid, pn in zip(pathwayId, pathwayName):
            tf_list_str = ''
            for tf in tflist[pid]:
                tf_list_str += '(:name %s) ' % tf.encode('ascii', 'ignore')

            tf_list_str = '(' + tf_list_str + ')'
            pathway_list_str += '(:name %s :tfs %s) ' % (pn, tf_list_str)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_common_tfs_genes(self, content):
        """Response content to FIND_COMMON_TF_GENES request
        For a given target list, reply the tfs regulating these genes
        and the frequency of each TF"""
        target_arg = content.gets('target')
        targets = _get_targets(target_arg)
        target_names = []
        for target in targets:
            target_names.append(target.name)

        tf_names, tf_counts = self.tfta.find_tfs_count(target_names)
        tf_list_str = ''
        for tf,ct in zip(tf_names, tf_counts):
            tf_list_str += '(:name %s :count %s) ' % (tf, ct)
        reply = KQMLList.from_string(
            '(SUCCESS :tfs (' + tf_list_str + '))')
        return reply

    def respond_find_overlap_targets_tf_genes(self, content):
        """For given tf list, find the targets which are regulated
        by all of them; then return the overlap between the targets
        and the given gene list"""
        tf_arg = content.gets('tf')
        tfs = _get_targets(tf_arg)
        tf_names = []
        for tf in tfs:
            tf_names.append(tf.name)

        target_arg = content.gets('target')
        targets = _get_targets(target_arg)
        target_names = []
        for tg in targets:
            target_names.append(tg.name)
        
        try:
            overlap_targets = \
                self.tfta.find_overlap_targets_tfs_genes(tf_names, target_names)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
	    return reply
			
        if len(overlap_targets):
            target_list_str = ''
            for tg in overlap_targets:
                target_list_str += '(:name %s ) ' % tg.encode('ascii', 'ignore')

            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        return reply
        
    def respond_find_common_pathway_genes(self, content):
	'''response content to FIND_COMMON_PATHWAY_GENES request'''
	target_arg = content.gets('target')
	targets = self._get_targets(target_arg)
	target_names = []
	for tg in targets:
            target_names.append(tg.name)
		
        try:
            pathwayName,externalId,source,dblink,counts = self.tfta.find_pathway_count_genes(target_names)
	except PathwayNotFoundException:
		reply = make_failure('PathwayNotFoundException')
		return reply
			
        path_list_str = ''
	for pn,eid,src,dbl,ct in zip(pathwayName,externalId,source,dblink,counts):
		path_list_str += '(:name %s :externalId %s :source %s :dblink %s :count %s) ' % (pn, eid ,src, dbl, ct)

	reply = KQMLList.from_string(
               '(SUCCESS :pathway (' + path_list_str + '))')
        return reply

def _get_target(target_str):
    target_str = '<ekb>' + target_str + '</ekb>'
    tp = TripsProcessor(target_str)
    terms = tp.tree.findall('TERM')
    term_id = terms[0].attrib['id']
    agent = tp._get_agent_by_id(term_id, None)
    return agent

def _get_targets(target_arg):
    print target_arg
    target_args = target_arg.split('</TERM>')
    agent = []
    for i in range(len(target_args)-1):
        target_str = '<ekb>' + target_args[i] + '</TERM></ekb>'
        tp = TripsProcessor(target_str)
        terms = tp.tree.findall('TERM')
        term_id = terms[0].attrib['id']
        agent.append(tp._get_agent_by_id(term_id, None))
    return agent

def make_failure(reason):
    msg = KQMLList('FAILURE')
    msg.set('reason', reason)
    return msg

if __name__ == "__main__":
    TFTA_Module(['-name', 'TFTA'] + sys.argv[1:])

