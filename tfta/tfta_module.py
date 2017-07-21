"""TFTA module is to receive and decode messages and send responses from
and to other agents in the system"""

import sys
import logging
from kqml import KQMLModule, KQMLPerformative, KQMLList
from tfta import TFTA, TFNotFoundException, TargetNotFoundException, PathwayNotFoundException, GONotFoundException
from indra.trips.processor import TripsProcessor

logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
		    level=logging.INFO)
logger = logging.getLogger('TFTA')

class TFTA_Module(KQMLModule):
    """TFTA module is used to receive and decode messages and send
    responses from and to other agents in the system."""
    def __init__(self, argv):
	# Call the constructor of KQMLModule
        super(TFTA_Module, self).__init__(argv)
	#Instantiate a singleton TFTA agent
	self.tfta = TFTA()
	
        self.tasks = ['IS-TF-TARGET', 'FIND-TF-TARGET', 'FIND-TARGET-TF',
                      'FIND-PATHWAY-GENE', 'FIND-PATHWAY-DB-GENE',
                      'FIND-TF-PATHWAY', 'FIND-GENE-PATHWAY',
                      'FIND-PATHWAY-KEYWORD', 'FIND-TF-KEYWORD',
                      'FIND-COMMON-TF-GENES', 'FIND-OVERLAP-TARGETS-TF-GENES',
                      'FIND-COMMON-PATHWAY-GENES','FIND-PATHWAY-GENE-KEYWORD',
		      'FIND-COMMON-PATHWAY-GENES-KEYWORD', 'FIND-GENE-GO-TF',
                      'IS-TF-TARGET-TISSUE', 'FIND-TF-TARGET-TISSUE',
                      'FIND-TARGET-TF-TISSUE']

        #Send subscribe messages
        for task in self.tasks:
            #msg_txt = '(subscribe :content (request &key :content (%s . *)))' % task
            #self.send(KQMLPerformative.from_string(msg_txt))
	    self.subscribe_request(task)
	
        #send ready message
        self.ready()
        #super(TFTA_Module, self).start()
	self.start()

    def receive_request(self, msg, content):
        """If a "request" message is received, decode the task and
        the content and call the appropriate function to prepare the
        response. A reply message is then sent back.
        """
        task_str = content.head().upper()
        if task_str == 'IS-TF-TARGET':
            reply_content = self.respond_is_tf_target(content)
        elif task_str == 'FIND-TF-TARGET':
            reply_content = self.respond_find_tf_targets(content)
        elif task_str == 'FIND-TARGET-TF':
            reply_content = self.respond_find_target_tfs(content)
        elif task_str == 'FIND-PATHWAY-GENE':
            reply_content = self.respond_find_pathway_gene(content)
	elif task_str == 'FIND-PATHWAY-GENE-KEYWORD':
            reply_content = self.respond_find_pathway_gene_keyword(content)
        elif task_str == 'FIND-PATHWAY-DB-GENE':
            reply_content = self.respond_find_pathway_db_gene(content)
        elif task_str == 'FIND-TF-PATHWAY':
            reply_content = self.respond_find_tf_pathway(content)
        elif task_str == 'FIND-GENE-PATHWAY':
            reply_content = self.respond_find_gene_pathway(content)
        elif task_str == 'FIND-PATHWAY-KEYWORD':
            reply_content = self.respond_find_pathway_chemical(content)
        elif task_str == 'FIND-TF-KEYWORD':
            reply_content = self.respond_find_tf_chemical(content)
        elif task_str == 'FIND-COMMON-TF-GENES':
            reply_content = self.respond_find_common_tfs_genes(content)
        elif task_str == 'FIND-OVERLAP-TARGETS-TF-GENES':
            reply_content = self.respond_find_overlap_targets_tf_genes(content)
        elif task_str == 'FIND-COMMON-PATHWAY-GENES':
	    reply_content = self.respond_find_common_pathway_genes(content)
	elif task_str == 'FIND-COMMON-PATHWAY-GENES-KEYWORD':
	    reply_content = self.respond_find_common_pathway_genes_keyword(content)
	elif task_str == 'FIND-GENE-GO-TF':
            reply_content = self.respond_find_genes_go_tf(content)
        elif task_str == 'IS-TF-TARGET-TISSUE':
            reply_content = self.respond_is_tf_target_tissue(content)
        elif task_str == 'FIND-TF-TARGET-TISSUE':
            reply_content = self.respond_find_tf_targets_tissue(content)
        elif task_str == 'FIND-TARGET-TF-TISSUE':
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
        #print 'tf=' + tf.name

        target_arg = content.gets('target')
        target = _get_target(target_arg)
        target_name = target.name

        tissue_arg = content.get('tissue')
        tissue_name = tissue_arg.head()
	tissue_name = trim_quotes(tissue_name)
	tissue_name = tissue_name.lower()
	tissue_name = tissue_name.replace(' ', '_')

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
	try:
            tfs = _get_targets(tf_arg)
            tf_names = []
            for tf in tfs:
                tf_names.append(tf.name)
	except Exception as e:
	    print 'received message:' + content
	    reply = make_failure('TF_NOT_FOUND')
	    return reply
            
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
	
	print "reply=",reply
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
	tissue_name = trim_quotes(tissue_name)
	tissue_name = tissue_name.lower()
	tissue_name = tissue_name.replace(' ', '_')

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
	try:
            targets = _get_targets(target_arg)
            target_names = []
            for target in targets:
                target_names.append(target.name)
	except Exception as e:
	    reply = make_failure('TARGET_NOT_FOUND')
	    return reply
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
	print "reply=", reply
        return reply

    def respond_find_target_tfs_tissue(self, content):
        """Response content to find-target-tf request
        For a target list, reply the tfs found within a given tissue"""
        target_arg = content.gets('target')
        targets = _get_targets(target_arg)
	target_names = []
	for target in targets:
            target_names.append(target.name)

        tissue_arg = content.get('tissue')
        tissue_name = tissue_arg.head()
	tissue_name = trim_quotes(tissue_name)
	tissue_name = tissue_name.lower()
	tissue_name = tissue_name.replace(' ', '_')

        if tissue_name not in tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        
        try:
            tf_names = self.tfta.find_tfs_tissue(target_names, tissue_name)
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
	try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
	except Exception as e:
	    reply = make_failure('GENE_NOT_FOUND')
	    return reply

        try:
            pathwayId, pathwayName, externalId, source,dblink = \
                self.tfta.find_pathways_from_genelist(gene_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName, externalId, source, dblink):
	    pnslash = '"' + pn + '"'
	    eidslash = '"' + eid + '"'
	    dbl = '"' + dbl + '"'
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pnslash, eidslash ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_gene_keyword(self,content):
        """Response content to find_pathway_gene_name request
        For a given gene list and keyword, reply the related pathways information"""
        keyword_arg = content.get('keyword')
        keyword_name = keyword_arg.head()
        keyword = trim_quotes(keyword_name)
        
        gene_arg = content.gets('gene')
	try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
	except Exception as e:
	    reply = make_failure('GENE_NOT_FOUND')
	    return reply

        try:
            pathwayId, pathwayName, externalId, source,dblink = \
                self.tfta.find_pathways_from_genelist_keyword(gene_names, keyword)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName, externalId, source, dblink):
	    pnslash = '"' + pn +'"'
	    eidslash= '"' + eid +'"'
	    dbl = '"' + dbl +'"'
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pnslash, eidslash ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_db_gene(self, content):
        """Response content to FIND-PATHWAY-DB-GENE request
        For a given gene list and certain db source, reply the related
        pathways information"""
        db_arg = content.get('database')
        db_name = db_arg.head()
	print db_name
	db_name = trim_quotes(db_name)

        gene_arg = content.gets('gene')
	try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
	except Exception as e:
	    reply = make_failure('GENE_NOT_FOUND')
	    return reply

        try:
            pathwayId, pathwayName, externalId, source, dblink = \
                self.tfta.find_pathways_from_dbsource_geneName(db_name,gene_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName, externalId, source, dblink):
            pnslash = '"' + pn + '"'
	    eidslash= '"' + eid + '"'
	    dbl = '"' + dbl + '"'
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pnslash, eidslash ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tf_pathway(self, content):
        """Response content to FIND_TF_PATHWAY request
        For a given pathway name, reply the tfs within the pathway"""
        pathway_arg = content.get('pathway')
        #pathway_name = pathway_arg.head()
	pathway_name = _get_targets(pathway_arg)

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
	    pnslash = '"' + pn +'"'
	    #eidslash= '"' + eid +'"'
            pathway_list_str += '(:name %s :tfs %s) ' % (pnslash, tf_list_str)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_gene_pathway(self, content):
        """Response content to FIND-GENE-PATHWAY request
        For a given pathway name, reply the genes within the pathway"""
        pathway_arg = content.gets('pathway')
        #pathway_name = pathway_arg.head()
        target = _get_target(pathway_arg)
        pathway_name = target.name

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
            pnslash = '"' + pn +'"'
            pathway_list_str += '(:name %s :genes %s) ' % (pnslash, gene_list_str)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_chemical(self, content):
        """Response content to FIND_PATHWAY_KEYWORD request
        For a given chemical name, reply the pathways involving the chemical"""
        chemical_arg = content.get('keyword')
        chemical_name = chemical_arg.head()
	chemical_name = trim_quotes(chemical_name)

        try:
            pathwayId, pathwayName, externalId, source, dblink = \
                self.tfta.find_pathways_from_chemical(chemical_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply

        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName,externalId,source,dblink):
            pnslash = '"' + pn +'"'
	    eidslash= '"' + eid +'"'
	    dbl = '"' + dbl +'"'
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pnslash, eidslash ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tf_chemical(self, content):
        """Response content to FIND_TF_KEYWORD request
        For a given chemical name, reply the tfs within the pathways
        involving the chemical"""
        chemical_arg = content.get('keyword')
        chemical_name = chemical_arg.head()
	chemical_name = trim_quotes(chemical_name)

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
            pnslash = '"' + pn + '"'
            pathway_list_str += '(:name %s :tfs %s) ' % (pnslash, tf_list_str)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_common_tfs_genes(self, content):
        """Response content to FIND-COMMON-TF-GENES request
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

    def respond_find_common_tfs_genes2(self, content):
        """Response content to FIND-COMMON-TF-GENES request
        For a given target list, reply the tfs regulating all these genes
        ,same as find_target_tf"""
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
        """response content to FIND-COMMON-PATHWAY-GENES request"""
	target_arg = content.gets('target')
	targets = _get_targets(target_arg)
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
            pnslash = '"' + pn +'"'
	    eidslash= '"' + eid +'"'
	    dbl = '"' + dbl +'"'
            path_list_str += '(:name %s :externalId %s :source %s :dblink %s :count %s) ' % (pnslash, eidslash ,src, dbl, ct)

	reply = KQMLList.from_string(
               '(SUCCESS :pathways (' + path_list_str + '))')
        return reply

    def respond_find_common_pathway_genes2(self, content):
        """response content to FIND-COMMON-PATHWAY-GENES request,same as find-pathway-gene
	For a given gene list, reply the related pathways information"""
        gene_arg = content.gets('target')
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
	    pnslash = '"' + pn +'"'
	    eidslash= '"' + eid +'"'
	    dbl = '"' + dbl +'"'
            pathway_list_str += \
                '(:name %s :externalId %s :source %s :dblink %s) ' % (pnslash, eidslash ,src, dbl)

        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_common_pathway_genes_keyword(self, content):
        """response content to FIND-COMMON-PATHWAY-GENES-KEYWORD request"""
	keyword_arg = content.get('keyword')
        keyword_name = keyword_arg.head()
        keyword = trim_quotes(keyword_name)
	
	target_arg = content.gets('target')
	targets = _get_targets(target_arg)
	target_names = []
	for tg in targets:
            target_names.append(tg.name)
		
        try:
            pathwayName,externalId,source,dblink,counts = self.tfta.find_pathway_count_genes_keyword(target_names, keyword)
	except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
			
        path_list_str = ''
	for pn,eid,src,dbl,ct in zip(pathwayName,externalId,source,dblink,counts):
            pnslash = '"' + pn +'"'
	    eidslash= '"' + eid +'"'
	    dbl = '"' + dbl +'"'
            path_list_str += '(:name %s :externalId %s :source %s :dblink %s :count %s) ' % (pnslash, eidslash ,src, dbl, ct)

	reply = KQMLList.from_string(
               '(SUCCESS :pathways (' + path_list_str + '))')
        return reply

    def respond_find_genes_go_tf(self, content):
	"""
        Respond content to FIND-GENE-GO-TF request
        """
        keyword_arg = content.get('keyword')
        keyword_name = keyword_arg.head()
        keyword = trim_quotes(keyword_name)
        keyword = keyword.lower()
        
        tf_arg = content.gets('tf')
        tfs = _get_targets(tf_arg)
        tf_names = []
        for tf in tfs:
            tf_names.append(tf.name)
            
        try:
            go_ids,go_types,go_names,go_genes = \
                self.tfta.find_genes_GO_tf(keyword, tf_names)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        except GONotFoundException:
            reply = make_failure('GOTERM_NOT_FOUND')
            return reply
            
        go_list_str = ''
        for gid,gn in zip(go_ids, go_names):
            gene_list_str = ''
            for gene in go_genes[gid]:
                gene_list_str += '(:name %s) ' % gene.encode('ascii', 'ignore')
                
            gene_list_str = '(' + gene_list_str + ')'
            gn_str = '"' + gn + '"'
            go_list_str += '(:name %s :genes %s) ' % (gn_str, gene_list_str)
            
        reply = KQMLList.from_string(
            '(SUCCESS :GOTERMs (' + go_list_str + '))')
            
        return reply

def _get_target(target_str):
    tp = TripsProcessor(target_str)
    terms = tp.tree.findall('TERM')
    term_id = terms[0].attrib['id']
    agent = tp._get_agent_by_id(term_id, None)
    return agent

def _get_targets(target_arg):
    tp = TripsProcessor(target_arg)
    agent = []
    for term in tp.tree.findall('TERM'):
        term_id = term.attrib['id']
        agent.append(tp._get_agent_by_id(term_id, None))
    return agent

def make_failure(reason):
    msg = KQMLList('FAILURE')
    msg.set('reason', reason)
    return msg

def trim_quotes(descr):
    if descr[0] == '"':
	descr = descr[1:]
    if descr[-1] == '"':
	descr = descr[:-1]
    return descr


if __name__ == "__main__":
    TFTA_Module(['-name', 'TFTA'] + sys.argv[1:])

