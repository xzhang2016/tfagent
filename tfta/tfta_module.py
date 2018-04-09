"""TFTA module is to receive and decode messages and send responses from
and to other agents in the system"""

import sys
import re
import time
import xml.etree.ElementTree as ET
import logging
logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA')
from kqml import KQMLModule, KQMLPerformative, KQMLList
from .tfta import TFTA, TFNotFoundException, TargetNotFoundException, PathwayNotFoundException 
from .tfta import GONotFoundException, miRNANotFoundException, TissueNotFoundException
from indra.sources.trips.processor import TripsProcessor
from collections import defaultdict
from bioagents import Bioagent

class TFTA_Module(Bioagent):
    """TFTA module is used to receive and decode messages and send
    responses from and to other agents in the system."""
    name = 'TFTA'
    tasks = ['IS-TF-TARGET', 'FIND-TF-TARGET', 'FIND-TARGET-TF',
             'FIND-PATHWAY-GENE', 'FIND-PATHWAY-DB-GENE',
             'FIND-TF-PATHWAY', 'FIND-GENE-PATHWAY',
             'FIND-PATHWAY-KEYWORD', 'FIND-TF-KEYWORD',
             'FIND-COMMON-TF-GENES', 'FIND-OVERLAP-TARGETS-TF-GENES',
             'FIND-COMMON-PATHWAY-GENES','FIND-PATHWAY-GENE-KEYWORD',
             'FIND-COMMON-PATHWAY-GENES-KEYWORD', 'FIND-GENE-GO-TF',
             'IS-TF-TARGET-TISSUE', 'FIND-TF-TARGET-TISSUE',
             'FIND-TARGET-TF-TISSUE', 'IS-PATHWAY-GENE','IS-MIRNA-TARGET', 
             'FIND-MIRNA-TARGET', 'FIND-TARGET-MIRNA', 'FIND-EVIDENCE-MIRNA-TARGET',
             'FIND-MIRNA-COUNT-GENE','FIND-GENE-COUNT-MIRNA',
             'FIND-PATHWAY-DB-KEYWORD', 'FIND-TISSUE', 'IS-REGULATION',
             'FIND-TF', 'FIND-PATHWAY', 'FIND-TARGET', 'FIND-GENE', 'FIND-MIRNA',
             'IS-GENE-ONTO', 'FIND-GENE-ONTO']

    #keep the genes from the most recent previous call, which are used to input 
    #find-gene-onto if there's no gene input 
    #gene_list = ['STAT3', 'JAK1', 'JAK2', 'ELK1', 'FOS', 'SMAD2', 'KDM4B']
    gene_list = []
    #keep track of TFs
    #global_tf_list = []

    def __init__(self, **kwargs):
        #Instantiate a singleton TFTA agent
        self.tfta = TFTA()
        # Call the constructor of KQMLModule
        super(TFTA_Module, self).__init__(**kwargs)
        
    def receive_tell(self, msg, content):
        #handle tell broadcast
        #now just do nothing here, but to avoid error message senting out
        pass

    def receive_request(self, msg, content):
        """If a "request" message is received, decode the task and
        the content and call the appropriate function to prepare the
        response. A reply message is then sent back.
        """
        task_str = content.head().upper()
        if task_str == 'IS-REGULATION':
            reply_content = self.respond_is_regulation(content)
        elif task_str == 'FIND-TF':
            reply_content = self.respond_find_tf(content)
        elif task_str == 'FIND-PATHWAY':
            reply_content = self.respond_find_pathway(content)
        elif task_str == 'FIND-TARGET':
            reply_content = self.respond_find_target(content)
        elif task_str == 'FIND-GENE':
            reply_content = self.respond_find_gene(content)
        elif task_str == 'FIND-MIRNA':
            reply_content = self.respond_find_miRNA(content)
        elif task_str == 'FIND-COMMON-PATHWAY-GENES':
            reply_content = self.respond_find_common_pathway_genes2keyword(content)
        elif task_str == 'FIND-MIRNA-COUNT-GENE':
            reply_content = self.respond_find_miRNA_count_target(content)
        elif task_str == 'FIND-GENE-COUNT-MIRNA':
            reply_content = self.respond_find_target_count_miRNA(content)
        elif task_str == 'IS-MIRNA-TARGET':
            reply_content = self.respond_is_miRNA_target(content)
        elif task_str == 'FIND-COMMON-TF-GENES':
            reply_content = self.respond_find_common_tfs_genes2(content)
        elif task_str == 'IS-PATHWAY-GENE':
            reply_content = self.respond_is_pathway_gene(content)
        elif task_str == 'FIND-TARGET-MIRNA':
            reply_content = self.respond_find_target_miRNA(content)
        elif task_str == 'FIND-MIRNA-TARGET':
            reply_content = self.respond_find_miRNA_target(content)
        elif task_str == 'IS-GENE-ONTO':
            reply_content = self.respond_is_gene_onto(content)
        elif task_str == 'FIND-GENE-ONTO':
            reply_content = self.respond_find_gene_onto(content)
        elif task_str == 'IS-TF-TARGET':
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
        elif task_str == 'FIND-OVERLAP-TARGETS-TF-GENES':
            reply_content = self.respond_find_overlap_targets_tf_genes(content)
        elif task_str == 'FIND-COMMON-PATHWAY-GENES-KEYWORD':
            reply_content = self.respond_find_common_pathway_genes_keyword(content)
        elif task_str == 'FIND-GENE-GO-TF':
            reply_content = self.respond_find_genes_go_tf2(content)
        elif task_str == 'IS-TF-TARGET-TISSUE':
            reply_content = self.respond_is_tf_target_tissue(content)
        elif task_str == 'FIND-TF-TARGET-TISSUE':
            reply_content = self.respond_find_tf_targets_tissue(content)
        elif task_str == 'FIND-TARGET-TF-TISSUE':
            reply_content = self.respond_find_target_tfs_tissue(content)
        elif task_str == 'FIND-EVIDENCE-MIRNA-TARGET':
            reply_content = self.respond_find_evidence_miRNA_target(content)
        elif task_str == 'FIND-PATHWAY-DB-KEYWORD':
            reply_content = self.respond_find_pathway_db_keyword(content)
        elif task_str == 'FIND-TISSUE':
            reply_content = self.respond_find_tissue_gene(content)
        else:
            self.error_reply(msg, 'unknown request task ' + task_str)
            return

        reply_msg = KQMLPerformative('reply')
        reply_msg.set('content', reply_content)
        self.reply(msg, reply_msg)
        
    def respond_is_regulation(self, content):
        """
        Response content to is-regulation request which includes:
        is-tf-target
        is-tf-target-tissue
        is-mirna-target (not covered here)
        """
        tf_arg = content.gets('tf')
        target_arg = content.gets('target')
        tissue_arg = content.get('tissue')
        mirna_arg = content.gets('mirna')
        if all([tf_arg,target_arg,tissue_arg]):
            reply = self.respond_is_tf_target_tissue(content)
        elif all([tf_arg,target_arg]):
            reply = self.respond_is_tf_target(content)
        #elif all([mirna_arg,target_arg]):
            #reply = self.respond_is_miRNA_target(content)
        else:
            reply = make_failure('UNKNOW_TASK')
        return reply
    
    def respond_find_tf(self, content):
        """
        Response content to find-tf request which includes:
        find-target-tf, find-target-tf-tissue, find-tf-keyword, find-tf-pathway,
        find-common-tf-genes (not covered)
        """
        target_arg = content.gets('target')
        pathway_arg = content.gets('pathway')
        keyword_arg = content.get('keyword')
        tissue_arg = content.get('tissue')
        #count_arg = content.get('count')
        if all([target_arg,tissue_arg]):
            reply = self.respond_find_target_tfs_tissue(content)
        #elif all([target_arg,count_arg]):
            #reply = self.respond_find_common_tfs_genes(content)
        elif target_arg:
            reply = self.respond_find_target_tfs(content)
        elif pathway_arg:
            reply = self.respond_find_tf_pathway(content)
        elif keyword_arg:
            reply = self.respond_find_tf_chemical(content)
        else:
            reply = make_failure('UNKNOW_TASK')
        return reply
    
    def respond_find_pathway(self, content):
        """
        Response content to find-pathway request, which includes:
        find-pathway-gene, find-pathway-db-gene, find-pathway-gene-keyword,
        find-pathway-keyword, find-pathway-db-keyword.
        """
        gene_arg = content.gets('gene')
        pathway_arg = content.gets('pathway')
        db_arg = content.get('database')
        keyword_arg = content.get('keyword')
        #count_arg = content.get('count')
        if all([keyword_arg,gene_arg]):
            reply = self.respond_find_pathway_gene_keyword(content)
        elif all([keyword_arg,db_arg]):
            reply = self.respond_find_pathway_db_keyword(content)
        elif keyword_arg:
            reply = self.respond_find_pathway_chemical(content)
        elif all([gene_arg,db_arg]):
            reply = self.respond_find_pathway_db_gene(content)
        elif gene_arg:
            reply = self.respond_find_pathway_gene(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
    
    def respond_find_target(self, content):
        """
        Response content to find-target request, which includes these cases:
        find-tf-target, find-tf-target-tissue, find-target-miran (not covered) 
        """
        tf_arg = content.gets('tf')
        target_arg = content.gets('target')
        miRNA_arg = content.gets('miRNA')
        tissue_arg = content.get('tissue')
        #count_arg = content.get('count')
        if all([tf_arg,tissue_arg]):
            reply = self.respond_find_tf_targets_tissue(content)
        elif tf_arg:
            reply = self.respond_find_tf_targets(content)
        #elif miRNA_arg:
            #reply = self.respond_find_target_miRNA(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
    
    def respond_find_gene(self, content):
        """
        Response content to find-gene request, which includes:
        find-gene-pathway, find-gene-go-tf
        """
        pathway_arg = content.gets('pathway')
        goid_arg = content.get('goid')
        tf_arg = content.gets('tf')
        if all([goid_arg,tf_arg]):
            reply = self.respond_find_genes_go_tf2(content)
        elif pathway_arg:
            reply = self.respond_find_gene_pathway(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
      
    def respond_find_miRNA(self, content):
        """
        Response content to find-mirna request, which includes:
        find-mirna-target, find-mirna-count-gene
        """
        target_arg = content.gets('target')
        count_arg = content.get('count')
        if all([target_arg,count_arg]):
            reply = self.respond_find_miRNA_count_target(content)
        elif target_arg:
            reply = self.respond_find_miRNA_target(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
        
    def respond_find_common_pathway_genes2keyword(self, content):
        """
        Respond to find-common-pathway-genes, which covers:
        find-common-pathway-genes, find-common-pathway-genes-keyword
        """
        target_arg = content.gets('target')
        keyword_arg = content.get('keyword')
        if all([target_arg, keyword_arg]):
            reply = self.respond_find_common_pathway_genes_keyword2(content)
        elif target_arg:
            reply = self.respond_find_common_pathway_genes2(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
        
    def respond_is_gene_onto(self, content):
        """
        Respond to is-gene-onto
        """
        try:
            keyword_arg = content.get('keyword')
            keyword_name = keyword_arg.data
            keyword_name = keyword_name.lower()
        except Exception as e:
            reply = make_failure('NO_KEYWORD')
            return reply
        try:
            gene_arg = content.gets('gene')
            gene = _get_target(gene_arg)
            gene_name = gene.name
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            is_onto = self.tfta.Is_gene_onto(keyword_name, gene_name)
        except GONotFoundException:
            reply = make_failure('GO_NOT_FOUND')
            return reply
        reply = KQMLList('SUCCESS')
        is_onto_str = 'TRUE' if is_onto else 'FALSE'
        reply.set('result', is_onto_str)
        return reply
        
    def respond_find_gene_onto(self, content):
        """
        Respond to find-gene-onto
        """
        try:
            keyword_arg = content.get('keyword')
            keyword_name = keyword_arg.data
            keyword_name = keyword_name.lower()
        except Exception as e:
            reply = make_failure('NO_KEYWORD')
            return reply
        try:
            gene_arg = content.get('gene')
            gene_arg_str = gene_arg.data
            gene_arg_str = gene_arg_str.replace(' ', '')
            gene_arg_str = gene_arg_str.upper()
            gene_names = gene_arg_str.split(',')
            #genes = _get_targets(gene_arg)
            #gene_names = []
            #for gene in genes:
                #gene_names.append(gene.name)
        except Exception as e:
            if self.gene_list:
                gene_names = self.gene_list
            else:
                reply = make_failure('NO_GENE_NAME')
                return reply
        try:
            results = self.tfta.find_gene_onto(keyword_name, gene_names)
        except GONotFoundException:
            reply = make_failure('GO_NOT_FOUND')
            return reply
        if len(results):
            gene_list_str = ''
            for res in results:
                gene_list_str += '(:name %s) ' % res
            reply = KQMLList.from_string(
                '(SUCCESS :genes (' + gene_list_str + '))')
        else:
            reply = make_failure('NO_GENE_FOUND')
        return reply
        
    def respond_is_tf_target(self, content):
        """
        Response content to is-tf-target request.
        """
        tf_arg = content.gets('tf')
        try:
            tf = _get_target(tf_arg)
            tf_name = tf.name
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply
        target_arg = content.gets('target')
        try:
            target = _get_target(target_arg)
            target_name = target.name
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
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
        reply.set('result', is_target_str)
        return reply

    tissue_list = ['bladder','blood','bone','bone_marrow','brain','cervix',
                   'colon','eye','heart','kidney','larynx','liver','lung',
                   'lymph_node','mammary_gland','muscle','ovary','pancreas',
                   'peripheral_nervous_system','placenta','prostate','skin',
                   'small_intestine','soft_tissue','spleen','stomach',
                   'testis','thymus','tongue','uterus']


    def respond_is_tf_target_tissue(self, content):
        """
        Response content to is-tf-target-tissue request.
        """
        tf_arg = content.gets('tf')
        try:
            tf = _get_target(tf_arg)
            tf_name = tf.name
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply
        target_arg = content.gets('target')
        try:
            target = _get_target(target_arg)
            target_name = target.name
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
        try:
            tissue_arg = content.get('tissue')
            #tissue_name = tissue_arg.head()
            tissue_name = tissue_arg.data
            #tissue_name = trim_quotes(tissue_name)
            tissue_name = tissue_name.lower()
            tissue_name = tissue_name.replace(' ', '_')
            tissue_name = tissue_name.replace('-', '_')
        except Exception as e:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
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
        reply.set('result', is_target_str)
        return reply

    def respond_find_tf_targets(self, content):
        """Response content to find-tf-target request
        For a tf list, reply with the targets found
        """
        t0 = time.clock()
        tf_arg = content.gets('tf')
        try:
            tfs = _get_targets(tf_arg)
            tf_names = []
            for tf in tfs:
                tf_names.append(tf.name)
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply
        #consider an optional parameter for sequencing query
        of_targets_names = []
        of_targets_arg = content.get('of-targets')
        if of_targets_arg:
            of_targets_data = of_targets_arg.data
            of_targets_data = of_targets_data.replace(' ', '')
            of_targets_data = of_targets_data.upper()
            of_targets_names = of_targets_data.split(',')  
        try:
            target_names = self.tfta.find_targets(tf_names)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        #check if it's a sequencing query
        if of_targets_names:
            target_names = list(set(of_targets_names) & set(target_names))
        if len(target_names):
            target_list_str = ''
            for tg in target_names:
                target_list_str += '(:name %s ) ' % tg
            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        t1 = time.clock()
        f = open('time-test.txt','a')
        t01 = t1 - t0
        f.write('This query used time(s):' + str(t01) + '\n')
        f.write('=============test end==================\n')
        f.close()
        self.gene_list = target_names
        return reply

    def respond_find_tf_targets_tissue(self, content):
        """
        Response content to find-tf-target-tissue request
        For tf list, reply with the targets found within given tissue
        """
        tf_arg = content.gets('tf')
        try:
            tfs = _get_targets(tf_arg)
            tf_names = []
            for tf in tfs:
                tf_names.append(tf.name)
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply
        try:
            tissue_arg = content.get('tissue')
            #tissue_name = tissue_arg.head()
            tissue_name = tissue_arg.data
            #tissue_name = trim_quotes(tissue_name)
            tissue_name = tissue_name.lower()
            tissue_name = tissue_name.replace(' ', '_')
            tissue_name = tissue_name.replace('-', '_')
        except Exception as e:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        #consider an optional parameter for sequencing query
        of_targets_names = []
        of_targets_arg = content.get('of-targets')
        if of_targets_arg:
            of_targets_data = of_targets_arg.data
            of_targets_data = of_targets_data.replace(' ', '')
            of_targets_data = of_targets_data.upper()
            of_targets_names = of_targets_data.split(',')    
        try:
            target_names = self.tfta.find_targets_tissue(tf_names, tissue_name)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        #check if it's a sequencing query
        if of_targets_names:
            target_names = list(set(of_targets_names) & set(target_names))  
        if len(target_names):
            target_list_str = ''
            for tg in target_names:
                target_list_str += '(:name %s) ' % tg
            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        self.gene_list = target_names
        return reply

    def respond_find_target_tfs(self, content):
        """
        Response content to find-target-tf request
        For a target list, reply the tfs found
        """
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for target in targets:
                target_names.append(target.name)
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
        #consider an optional parameter for sequencing query
        of_tfs_names = []
        of_tfs_arg = content.get('of-tfs')
        if of_tfs_arg:
            of_tfs_data = of_tfs_arg.data
            of_tfs_data = of_tfs_data.replace(' ', '')
            of_tfs_data = of_tfs_data.upper()
            of_tfs_names = of_tfs_data.split(',')
        try:
            tf_names = self.tfta.find_tfs(target_names)
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply
        #check if it's a sequencing query
        if of_tfs_names:
            tf_names = list(set(of_tfs_names) & set(tf_names))   
        if len(tf_names):    
            tf_list_str = ''
            for tf in tf_names:
                tf_list_str += '(:name %s) ' % tf
            reply = KQMLList.from_string(
                '(SUCCESS :tfs (' + tf_list_str + '))')
        else:
            reply = make_failure('NO_TF_FOUND')
        self.gene_list = tf_names
        return reply

    def respond_find_target_tfs_tissue(self, content):
        """
        Response content to find-target-tf-tissue request
        For a target list, reply the tfs found within a given tissue
        """
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for target in targets:
                target_names.append(target.name)
            #debug: tract the target_names
            f = open('targets-find-tfs-tissue.txt', 'a')
            f.write(';'.join(target_names))
            f.write('\n' + '================' + '\n')
            f.close()
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
        try:
            tissue_arg = content.get('tissue')
            #tissue_name = tissue_arg.head()
            tissue_name = tissue_arg.data
            #tissue_name = trim_quotes(tissue_name)
            tissue_name = tissue_name.lower()
            tissue_name = tissue_name.replace(' ', '_')
            tissue_name = tissue_name.replace('-', '_')
        except Exception as e:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        #consider an optional parameter for sequencing query
        of_tfs_names = []
        of_tfs_arg = content.get('of-tfs')
        if of_tfs_arg:
            of_tfs_data = of_tfs_arg.data
            of_tfs_data = of_tfs_data.replace(' ', '')
            of_tfs_data = of_tfs_data.upper()
            of_tfs_names = of_tfs_data.split(',')
        try:
            tf_names = self.tfta.find_tfs_tissue(target_names, tissue_name)
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply
        #check if it's a sequencing query
        if of_tfs_names:
            tf_names = list(set(of_tfs_names) & set(tf_names))         
        if len(tf_names):
            tf_list_str = ''
            for tf in tf_names:
                tf_list_str += '(:name %s) ' % tf
            reply = KQMLList.from_string(
                '(SUCCESS :tfs (' + tf_list_str + '))')
        else:
            reply = make_failure('NO_TF_FOUND')
        self.gene_list = tf_names
        return reply

    def respond_find_pathway_gene(self,content):
        """
        Response content to find_pathway_gene request
        For a given gene list, reply the related pathways information
        """
        gene_arg = content.gets('gene')
        try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            pathwayName, dblink = \
                self.tfta.find_pathways_from_genelist(gene_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        pathway_list_str = ''
        for pn, dbl in zip(pathwayName, dblink):
            pnslash = '"' + pn + '"'
            dbl = '"' + dbl + '"'
            pathway_list_str += \
                '(:name %s :dblink %s) ' % (pnslash, dbl)
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_gene_keyword(self,content):
        """
        Response content to find_pathway_gene_name request
        For a given gene list and keyword, reply the related pathways information
        """
        keyword_arg = content.get('keyword')
        try:
            #keyword_name = keyword_arg.head()
            keyword_name = keyword_arg.data
        except Exception as e:
            reply = make_failure('NO_KEYWORD')
            return reply        
        gene_arg = content.gets('gene')
        try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            pathwayId, pathwayName, externalId, source,dblink = \
                self.tfta.find_pathways_from_genelist_keyword(gene_names, keyword_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        pathway_list_str = ''
        for pn, eid, src, dbl in zip(pathwayName, externalId, source, dblink):
            pnslash = '"' + pn +'"'
            dbl = '"' + dbl +'"'
            pathway_list_str += \
                '(:name %s :dblink %s) ' % (pnslash, dbl)
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_pathway_db_gene(self, content):
        """Response content to FIND-PATHWAY-DB-GENE request
        For a given gene list and certain db source, reply the related
        pathways information"""
        db_arg = content.get('database')
        try:
            #db_name = db_arg.head()
            db_name = db_arg.data
            #db_name = trim_quotes(db_name)
        except Exception as e:
            reply = make_failure('NO_DB_NAME')
            return reply
        gene_arg = content.gets('gene')
        try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
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
            dbl = '"' + dbl + '"'
            pathway_list_str += \
                '(:name %s :dblink %s) ' % (pnslash, dbl)
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tf_pathway(self, content):
        """
        Response content to FIND_TF_PATHWAY request
        For a given pathway name, reply the tfs within the pathway
        """
        pathway_arg = content.gets('pathway')
        pathway_names = _get_pathway_name(pathway_arg)
        pathway_names = trim_word(pathway_names, 'pathway')
        if not len(pathway_names):
            reply = make_failure('NO_PATHWAY_NAME')
            return reply    
        try:
            pathwayId,pathwayName,tflist,dblink = \
                self.tfta.find_tfs_from_pathwayName(pathway_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        pathway_list_str = ''
        keys = list(tflist.keys())
        temp_tfs = [] #track all the tfs
        if keys:
            for key in keys:
                temp_tfs += tflist[key]
                tf_list_str = ''
                for tf in tflist[key]:
                    tf_list_str += '(:name %s) ' % tf
                tf_list_str = '(' + tf_list_str + ')'
                pn = '"' + pathwayName[key] + '"'
                dl = '"' + dblink[key] + '"'
                pathway_list_str += '(:name %s :dblink %s :tfs %s) ' % (pn, dl, tf_list_str)
            reply = KQMLList.from_string(
                    '(SUCCESS :pathways (' + pathway_list_str + '))')
        else:
            reply = make_failure('PathwayNotFoundException')
            return reply
        self.gene_list = list(set(temp_tfs))
        return reply

    def respond_find_gene_pathway(self, content):
        """
        Response content to FIND-GENE-PATHWAY request
        For a given pathway name, reply the genes within the pathway
        """
        pathway_arg = content.gets('pathway')
        pathway_names = _get_pathway_name(pathway_arg)
        pathway_names = trim_word(pathway_names, 'pathway')
        if not len(pathway_names):
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        try:
            pathwayId,pathwayName,genelist,plink = \
                self.tfta.find_genes_from_pathwayName(pathway_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        pathway_list_str = ''
        temp_genes = []
        for pid in pathwayId:
            temp_genes += genelist[pid]
            gene_list_str = ''
            for gene in genelist[pid]:
                gene_list_str += '(:name %s) ' % gene
            gene_list_str = '(' + gene_list_str + ')'
            pnslash = '"' + pathwayName[pid] +'"'
            pwlink = '"' + plink[pid] + '"'
            pathway_list_str += '(:name %s :dblink %s :genes %s) ' % (pnslash, pwlink, gene_list_str)            
        reply = KQMLList.from_string(
                        '(SUCCESS :pathways (' + pathway_list_str + '))')
        self.gene_list = list(set(temp_genes))
        return reply

    def respond_find_pathway_chemical(self, content):
        """
        Response content to FIND_PATHWAY_KEYWORD request
        For a given chemical name, reply the pathways involving the chemical
        """
        chemical_arg = content.get('keyword')
        try:
            chemical_name = chemical_arg.data
        except Exception as e:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        try:
            pathwayId, pathwayName, externalId, source, dblink = \
                self.tfta.find_pathways_from_chemical(chemical_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        pathway_list_str = ''
        for pn, dbl in zip(pathwayName,dblink):
            pnslash = '"' + pn +'"'
            dbl = '"' + dbl +'"'
            pathway_list_str += \
                '(:name %s :dblink %s) ' % (pnslash, dbl)
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tf_chemical(self, content):
        """
        Response content to FIND_TF_KEYWORD request
        For a given chemical name, reply the tfs within the pathways
        involving the chemical
        """
        chemical_arg = content.get('keyword')
        try:
            chemical_name = chemical_arg.data
        except Exception as e:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        try:
            pathwayId, pathwayName, tflist, dblink = \
                self.tfta.find_tfs_from_pathwaysWithChemical(chemical_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        pathway_list_str = ''
        temp_tfs = []
        for pid, pn, dl in zip(pathwayId, pathwayName, dblink):
            tf_list_str = ''
            temp_tfs += tflist[pid]
            for tf in tflist[pid]:
                tf_list_str += '(:name %s) ' % tf
            tf_list_str = '(' + tf_list_str + ')'
            pnslash = '"' + pn + '"'
            dl = '"' + dl + '"'
            pathway_list_str += '(:name %s :dblink %s :tfs %s) ' % (pnslash, dl, tf_list_str)
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        self.gene_list = list(set(temp_tfs))
        return reply

    def respond_find_common_tfs_genes(self, content):
        """
        Response content to FIND-COMMON-TF-GENES request
        For a given target list, reply the tfs regulating these genes
        and the frequency of each TF
        """
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for target in targets:
                target_names.append(target.name)
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
        try:
            tf_counts = self.tfta.find_tfs_count(target_names)
        except TFNotFoundException:
            reply = make_failure('NO_TF_FOUND')
            return reply
        #cluster the tfs according to the count
        tf_clustered = cluster_dict_by_value(tf_counts)
        tf_list_str = ''
        counts = list(tf_clustered.keys())
        counts.reverse()
        for ct in counts:
            tf_list = ''
            for tf in tf_clustered[ct]:
                tf_list += '(:name %s) ' % tf
            tf_list = ':tf-list (' + tf_list + ')'
            tf_list += ' :count %d' % ct
            tf_list_str += '(' + tf_list + ') '
        reply = KQMLList.from_string(
            '(SUCCESS :tfs (' + tf_list_str + '))')
        return reply
        
    def respond_find_common_tfs_genes2(self, content):
        """
        Response content to FIND-COMMON-TF-GENES request
        For a given target list, reply the tfs regulating these genes
        and the regulated targets by each TF
        """
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for target in targets:
                target_names.append(target.name)
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
        try:
            tf_targets = self.tfta.find_common_tfs(target_names)
        except TFNotFoundException:
            reply = make_failure('NO_TF_FOUND')
            return reply
        #cluster the tfs according to the targets
        tf_clustered = cluster_dict_by_value2(tf_targets)
        tf_list_str = ''
        targets = list(tf_clustered.keys())
        for tg in targets:
            target_list = tg.split(',')
            target_str = ''
            for t in target_list:
                target_str += '(:name %s) ' % t
            target_str = ':target-list (' + target_str + ')'
            tf_list = ''
            for tf in tf_clustered[tg]:
                tf_list += '(:name %s) ' % tf
            tf_list = ':tf-list (' + tf_list + ') '
            tf_list += target_str
            tf_list_str += '(' + tf_list + ') '
        reply = KQMLList.from_string(
            '(SUCCESS :tfs (' + tf_list_str + '))')
        return reply

    def respond_find_overlap_targets_tf_genes(self, content):
        """
        For given tf list, find the targets which are regulated
        by all of them; then return the overlap between the targets
        and the given gene list
        """
        tf_arg = content.gets('tf')
        try:
            tfs = _get_targets(tf_arg)
            tf_names = []
            for tf in tfs:
                tf_names.append(tf.name)
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for tg in targets:
                target_names.append(tg.name)
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply        
        try:
            overlap_targets = \
                self.tfta.find_overlap_targets_tfs_genes(tf_names, target_names)
        except TFNotFoundException:
            reply = make_failure('TF_NOT_FOUND')
            return reply
        if len(overlap_targets):
            target_list_str = ''
            for tg in overlap_targets:
                target_list_str += '(:name %s ) ' % tg
            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        return reply
        
    def respond_find_common_pathway_genes(self, content):
        """
        response content to FIND-COMMON-PATHWAY-GENES request
        """
        target_arg = content.gets('target')
        if not target_arg:
            target_arg = content.gets('gene')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for tg in targets:
                target_names.append(tg.name)
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            pathwayName,externalId,source,dblink,counts = self.tfta.find_pathway_count_genes(target_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        path_list_str = ''
        for pn,eid,src,dbl,ct in zip(pathwayName,externalId,source,dblink,counts):
            pnslash = '"' + pn +'"'
            #eidslash= '"' + eid +'"'
            dbl = '"' + dbl +'"'
            path_list_str += '(:name %s :dblink %s :count %d) ' % (pnslash, dbl, ct)
        reply = KQMLList.from_string(
               '(SUCCESS :pathways (' + path_list_str + '))')
        return reply

    def respond_find_common_pathway_genes_keyword(self, content):
        """
        Response content to FIND-COMMON-PATHWAY-GENES-KEYWORD request
        """
        keyword_arg = content.get('keyword')
        try:
            #keyword_name = keyword_arg.head()
            keyword_name = keyword_arg.data
        except Exception as e:
            reply = make_failure('NO_KEYWORD')
            return reply
        target_arg = content.gets('gene')
        if not target_arg:
            target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for tg in targets:
                target_names.append(tg.name)
            target_names = list(set(target_names))
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            pathwayName,externalId,source,dblink,counts = \
               self.tfta.find_pathway_count_genes_keyword(target_names, keyword_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        path_list_str = ''
        for pn,dbl,ct in zip(pathwayName,dblink,counts):
            pnslash = '"' + pn +'"'
            dbl = '"' + dbl +'"'
            path_list_str += '(:name %s :dblink %s :count %d) ' % (pnslash, dbl, ct)
        reply = KQMLList.from_string(
               '(SUCCESS :pathways (' + path_list_str + '))')
        return reply
        
    def respond_find_common_pathway_genes2(self, content):
        """
        response content to FIND-COMMON-PATHWAY-GENES request
        """
        target_arg = content.gets('target')
        if not target_arg:
            target_arg = content.gets('gene')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for tg in targets:
                target_names.append(tg.name)
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            pathwayName, dblink, genes = self.tfta.find_common_pathway_genes(target_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        path_list_str = ''
        for pth in genes.keys():
            pnslash = '"' + pathwayName[pth] +'"'
            dbl = '"' + dblink[pth] +'"'
            gene_list_str = ''
            for g in genes[pth]:
                gene_list_str += '(:name %s) ' % g
            path_list_str += '(:name %s :dblink %s ' % (pnslash, dbl)
            path_list_str += ':gene-list (' + gene_list_str + ')) '
        reply = KQMLList.from_string(
               '(SUCCESS :pathways (' + path_list_str + '))')
        return reply

    def respond_find_common_pathway_genes_keyword2(self, content):
        """
        Response content to FIND-COMMON-PATHWAY-GENES-KEYWORD request
        """
        keyword_arg = content.get('keyword')
        try:
            #keyword_name = keyword_arg.head()
            keyword_name = keyword_arg.data
        except Exception as e:
            reply = make_failure('NO_KEYWORD')
            return reply
        target_arg = content.gets('gene')
        if not target_arg:
            target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for tg in targets:
                target_names.append(tg.name)
            target_names = list(set(target_names))
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            pathwayName,dblink,genes = \
               self.tfta.find_common_pathway_genes_keyword(target_names, keyword_name)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply
        path_list_str = ''
        for pth in genes.keys():
            pnslash = '"' + pathwayName[pth] +'"'
            dbl = '"' + dblink[pth] +'"'
            gene_list_str = ''
            for g in genes[pth]:
                gene_list_str += '(:name %s) ' % g
            path_list_str += '(:name %s :dblink %s ' % (pnslash, dbl)
            path_list_str += ':gene-list (' + gene_list_str + ')) '
        reply = KQMLList.from_string(
               '(SUCCESS :pathways (' + path_list_str + '))')
        return reply

    def respond_is_pathway_gene(self, content):
        """
        Respond to IS-PATHWAY-GENE request
        query like: Does the mTor pathway utilize SGK1? 
        """
        pathway_arg = content.gets('pathway')
        pathway_names = _get_pathway_name(pathway_arg)
        pathway_names = trim_word(pathway_names, 'pathway')
        if not len(pathway_names):
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        gene_arg = content.gets('gene')
        try:
            genes = _get_targets(gene_arg)
            gene_names = []
            for gene in genes:
                gene_names.append(gene.name)
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply  
        try:
            pids, pathwayName, dblink = \
                self.tfta.Is_pathway_gene(pathway_names, gene_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply            
        pathway_list_str = ''
        for pid in pids:
            pn = '"' + pathwayName[pid] + '"'
            dbl = '"' + dblink[pid] + '"'
            pathway_list_str += \
                '(:name %s :dblink %s) ' % (pn, dbl)                
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_genes_go_tf(self, content):
        """
        Respond content to FIND-GENE-GO-TF request
        The inputs here are keyword from GO names and tf list
        """
        keyword_arg = content.get('keyword')
        try:
            #keyword_name = keyword_arg.head()
            keyword_name = keyword_arg.data
            #keyword = trim_quotes(keyword_name)
            keyword = keyword_name.lower()
        except Exception as e:
            reply = make_failure('NO_GO_NAME')
            return reply        
        tf_arg = content.gets('tf')
        try:
            tfs = _get_targets(tf_arg)
            tf_names = []
            for tf in tfs:
                tf_names.append(tf.name)
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply            
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
                gene_list_str += '(:name %s) ' % gene
            gene_list_str = '(' + gene_list_str + ')'
            gn_str = '"' + gn + '"'
            gid_str = '"' + gid + '"'
            go_list_str += '(:id %s :name %s :genes %s) ' % (gid_str, gn_str, gene_list_str)            
        reply = KQMLList.from_string(
            '(SUCCESS :go-terms (' + go_list_str + '))')        
        return reply

    def respond_find_genes_go_tf2(self, content):
        """
        Respond content to FIND-GENE-GO-TF request
        The inputs here are GO id and tf list
        """
        go_arg = content.get('goid')
        try:
            #goid = go_arg.head()
            goid = go_arg.data
            #goid = trim_quotes(goid)
        except Exception as e:
            reply = make_failure('NO_GO_ID')
            return reply       
        tf_arg = content.gets('tf')
        try:
            tfs = _get_targets(tf_arg)
            tf_names = []
            for tf in tfs:
                tf_names.append(tf.name)
        except Exception as e:
            reply = make_failure('NO_TF_NAME')
            return reply           
        try:
            go_ids,go_types,go_names,go_genes = \
                self.tfta.find_genes_GO_tf2(goid, tf_names)
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
                gene_list_str += '(:name %s) ' % gene
            gene_list_str = '(' + gene_list_str + ')'
            gn_str = '"' + gn + '"'
            gid_str = '"' + gid + '"'
            go_list_str += '(:id %s :name %s :genes %s) ' % (gid_str, gn_str, gene_list_str)            
        reply = KQMLList.from_string(
            '(SUCCESS :go-terms (' + go_list_str + '))')            
        return reply

    def respond_is_miRNA_target(self, content):
        """
        Respond to IS-MIRNA-TARGET request
        """
        #assume the miRNA is also in EKB XML format
        miRNA_arg = content.gets('miRNA')
        miRNA_name = _get_miRNA_name(miRNA_arg)
        if not len(miRNA_name):
            reply = make_failure('NO_MIRNA_NAME')
            return reply       
        target_arg = content.gets('target')
        try:
            target = _get_target(target_arg)
            target_name = target.name
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply        
        try:
            is_target = self.tfta.Is_miRNA_target(miRNA_name[0], target_name)
        except miRNANotFoundException:
            reply = make_failure('MIRNA_NOT_FOUND')
            return reply
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply            
        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if is_target else 'FALSE'
        reply.set('is-miRNA-target', is_target_str)
        return reply
        
    def respond_find_miRNA_target(self, content):
        """
        Respond to FIND-MIRNA-TARGET request
        """
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for target in targets:
                target_names.append(target.name)
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply        
        try:
            miRNA_names = self.tfta.find_miRNA_target(target_names)
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply            
        if len(miRNA_names):
            miRNA_list_str = ''
            for m in miRNA_names:
                miRNA_list_str += '(:name %s ) ' % m
            reply = KQMLList.from_string(
                '(SUCCESS :miRNAs (' + miRNA_list_str + '))')
        else:
            reply = make_failure('NO_MIRNA_FOUND')            
        return reply
        
    def respond_find_target_miRNA(self, content):
        """
        Respond to FIND-TARGET-MIRNA request
        """
        #assume the miRNA is also in EKB XML format
        try:
            miRNA_arg = content.gets('miRNA')
            miRNA_names = _get_miRNA_name(miRNA_arg)
        except Exception as e:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        if not len(miRNA_names):
            reply = make_failure('NO_MIRNA_NAME')
            return reply 
        try:
            target_names = self.tfta.find_target_miRNA(miRNA_names)
        except miRNANotFoundException:
            #for user clarification
            if len(miRNA_names) == 1:
                try:
                    clari_mirna = self.tfta.get_similar_miRNAs(miRNA_names[0])
                    c_str = ''
                    for c in clari_mirna:
                        c_str += '(:name %s) ' % c
                    reply = make_failure_clarification('MIRNA_NOT_FOUND', '(' + c_str + ')')
                    return reply
                except miRNANotFoundException:
                    reply = make_failure('NO_SIMILAR_MIRNA')
                    return reply
            else:
                reply = make_failure('MIRNA_NOT_FOUND')
                return reply
        if len(target_names):
            target_list_str = ''
            for tg in target_names:
                target_list_str += '(:name %s ) ' % tg
            reply = KQMLList.from_string(
                '(SUCCESS :targets (' + target_list_str + '))')
        else:
            reply = make_failure('NO_TARGET_FOUND')
        self.gene_list = target_names  
        return reply
        
    def respond_find_evidence_miRNA_target(self, content):
        """
        Respond to FIND-EVIDENCE-MIRNA-TARGET request
        """
        #assume the miRNA is also in EKB XML format
        miRNA_arg = content.gets('miRNA')
        miRNA_name = _get_miRNA_name(miRNA_arg)
        if not len(miRNA_name):
            reply = make_failure('NO_MIRNA_NAME')
            return reply      
        target_arg = content.gets('target')
        try:
            target = _get_target(target_arg)
            target_name = target.name
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply           
        try:
            experiments,supportType,pmlink = \
                self.tfta.find_evidence_miRNA_target(miRNA_name[0], target_name)
        except miRNANotFoundException:
            reply = make_failure('MIRNA_NOT_FOUND')
            return reply
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply            
        if len(experiments):
            evidence_list_str = ''
            for e,s,l in zip(experiments,supportType,pmlink):
                en = '"' + e + '"'
                sn = '"' + s + '"'
                ln = '"' + l + '"'
                evidence_list_str += \
                    '(:experiment %s :supportType %s :pubmedLink %s) ' % (en, sn, ln)
            reply = KQMLList.from_string(
                '(SUCCESS :evidence (' + evidence_list_str + '))')
        else:
            reply = make_failure('NO_EVIDENCE_FOUND')   
        return reply

    def respond_find_target_count_miRNA(self, content):
        """
        Respond to FIND-GENE-COUNT-MIRNA request
        """
        #assume the miRNA is also in EKB XML format
        miRNA_arg = content.gets('miRNA')
        miRNA_names = _get_miRNA_name(miRNA_arg)
        if not len(miRNA_names):
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        try:
            targets,counts,mrna = self.tfta.find_gene_count_miRNA(miRNA_names)
        except miRNANotFoundException:
            reply = make_failure('MIRNA_NOT_FOUND')
            return reply
        except TargetNotFoundException:
            reply = make_failure('NO_TARGET_FOUND')
            return reply           
        target_str = ''
        for t,ct in zip(targets,counts):
            ms = mrna[t]
            m_str = ''
            for m in ms:
                m_str += '(:name %s)' % m
            m_str = '(' + m_str + ')'
            target_str += '(:name %s :count %d :miRNA %s) ' % (t, ct, m_str)
        reply = KQMLList.from_string(
            '(SUCCESS :targets (' + target_str + '))')
        self.gene_list = targets
        return reply
             
    def respond_find_miRNA_count_target(self, content):
        """
        Respond to FIND-MIRNA-COUNT-GENE request
        """
        target_arg = content.gets('target')
        try:
            targets = _get_targets(target_arg)
            target_names = []
            for tg in targets:
                target_names.append(tg.name)
        except Exception as e:
            reply = make_failure('NO_TARGET_NAME')
            return reply
        try:
            mirnas,counts,genes = self.tfta.find_miRNA_count_gene(target_names)
        except miRNANotFoundException:
            reply = make_failure('NO_MIRNA_FOUND')
            return reply
        except TargetNotFoundException:
            reply = make_failure('TARGET_NOT_FOUND')
            return reply            
        mirna_str = ''
        for m,ct in zip(mirnas,counts):
            gs = genes[m]
            g_str = ''
            for g in gs:
                g_str += '(:name %s) ' % g
            g_str = '(' + g_str + ')'
            mirna_str += '(:name %s :count %d :targets %s) ' % (m, ct, g_str)
        reply = KQMLList.from_string(
            '(SUCCESS :miRNAs (' + mirna_str + '))')
        return reply

    def respond_find_pathway_db_keyword(self, content):
        """
        Respond to FIND-PATHWAY-DB-KEYWORD request
        """
        try:
            db_arg = content.get('database')
            #db_name = db_arg.head()
            db_name = db_arg.data
        except Exception as e:
            reply = make_failure('NO_DB_NAME')
            return reply
        try:
            pathway_arg = content.get('keyword')
            pathway_names = pathway_arg.data
            pathway_names = pathway_names.replace('W::', '')
            pathway_names = pathway_names.replace('-', ' ')
            pathway_names = trim_word([pathway_names], 'pathway')
        except Exception as e:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        try:
            pathwayId,pathwayName,dblink = \
                    self.tfta.find_pathway_db_keyword(db_name, pathway_names)
        except PathwayNotFoundException:
            reply = make_failure('PathwayNotFoundException')
            return reply        
        pathway_list_str = ''
        for id in pathwayId:
            pn = '"' + pathwayName[id] + '"'
            dbl = '"' + dblink[id] + '"'
            pathway_list_str += \
                '(:name %s :dblink %s) ' % (pn, dbl)
        reply = KQMLList.from_string(
            '(SUCCESS :pathways (' + pathway_list_str + '))')
        return reply

    def respond_find_tissue_gene(self, content):
        """
        Response to FIND-TISSUE-GENE task
        """
        try:
            gene_arg = content.gets('gene')
            gene = _get_target(gene_arg)
            gene_name = gene.name
        except Exception as e:
            reply = make_failure('NO_GENE_NAME')
            return reply
        try:
            tissues = self.tfta.find_tissue_gene(gene_name)
        except TissueNotFoundException:
            reply = make_failure('NO_TISSUE_FOUND')
            return reply    
        tissue_str = ''
        for ts in tissues:
            tissue_str += '(:name %s) ' % ts
        reply = KQMLList.from_string(
            '(SUCCESS :tissue (' + tissue_str + '))')
        return reply

def _get_target(target_str):
    tp = TripsProcessor(target_str)
    terms = tp.tree.findall('TERM')
    term_id = terms[0].attrib['id']
    agent = tp._get_agent_by_id(term_id, None)
    return agent

def _get_targets(target_arg):
    agent = []
    ont1 = ['ONT::PROTEIN', 'ONT::GENE-PROTEIN', 'ONT::GENE']
    tp = TripsProcessor(target_arg)
    for term in tp.tree.findall('TERM'):
        if term.find('type').text in ont1:
            term_id = term.attrib['id']
            agent.append(tp._get_agent_by_id(term_id, None))
    return agent

def _get_pathway_name(target_str):
    #print('In _get_pathway_name')
    #tree = ET.XML(xml_string, parser=UTB())
    pathway_name = []
    try:
        #test the ekb xml format
        f = open('TFTA-test-pathway-ekb.txt', 'a')
        f.write(target_str + '\n')
        f.write('===============================\n')
        
        root = ET.fromstring(target_str)
    except Exception as e:
        return pathway_name
    try:
        for term in root.findall('TERM'):
            for t in term.find('drum-terms').findall('drum-term'):
                s = t.get('matched-name')
                if s and s not in ['PATHWAY', 'SIGNALING-PATHWAY']:
                    pathway_name = pathway_name + [t.get('matched-name')]
                s = t.get('name')
                if s and s not in ['PATHWAY', 'SIGNALING-PATHWAY']:
                    pathway_name = pathway_name + [t.get('name')]
        pathway_name = list(set(pathway_name))
    except Exception as e:
        try:
            for term in root.findall('TERM'):
                s = term.find('name')
                if s is not None:
                    s1 = s.text
                    if s1 not in ['PATHWAY', 'SIGNALING-PATHWAY']:
                        pathway_name = pathway_name + [s1.replace('-', ' ').lower()]
            pathway_name = list(set(pathway_name))
        except Exception as e:
            return pathway_name
    f.write('Extracted pathwayName=' + ';;;'.join(pathway_name) + '\n')
    f.write('===============================\n')
    f.close()
    return pathway_name

def _get_miRNA_name(xml_string):
    miRNA_names = []
    try:
        root = ET.fromstring(xml_string)
    except Exception as e:
        return miRNA_names
    try:
        for term in root.findall('TERM'):
            if term is not None:
                dts = term.find('drum-terms').findall('drum-term')
                for dt in dts:
                    if dt.get('matched-name') is not None:
                        #change miRNA name to standard name
                        s1 = dt.get('matched-name')
                        #matched_pattern = re.findall('([0-9]+-)[a-zA-Z]', s1)
                        s1 = rtrim_hyphen(s1)
                        miRNA_names.append(s1.upper())
        miRNA_names = list(set(miRNA_names))
    except Exception as e:
        try:
            for term in root.findall('TERM'):
                s = term.find('name')
                if s is not None:
                    s1 = s.text
                    #change miRNA name to standard name
                    #matched_pattern = re.findall('([0-9]+-)[a-zA-Z]', s1)
                    s1 = rtrim_hyphen(s1)
                    miRNA_names.append(s1.upper())
            miRNA_names = list(set(miRNA_names))
        except Exception as e:
            return miRNA_names
    return miRNA_names

def make_failure(reason):
    msg = KQMLList('FAILURE')
    msg.set('reason', reason)
    return msg

def make_failure_clarification(reason, clarification):
    msg = KQMLList('FAILURE')
    msg.set('reason', reason)
    msg.set('clarification', clarification)
    return msg

def trim_quotes(descr):
    if descr[0] == '(':
        descr = descr[1:]
    if descr[-1] == ')':
        descr = descr[:-1]
    if descr[0] == '"':
        descr = descr[1:]
    if descr[-1] == '"':
        descr = descr[:-1]
    return descr

def trim_hyphen(descr):
    if descr[0] == '-':
        descr = descr[1:]
    return descr

def trim_word(descr, word):
    #descr is a list
    ds = []
    if len(descr):
        for d in descr:
            if d[-len(word):] == word:
                if len(d[:-len(word)-1]):
                    ds.append(d[:-len(word)-1])
            else:
                ds.append(d)
    return ds

def rtrim_hyphen(str1):
    plist = re.findall('([0-9]+-[a-zA-Z])', str1)
    s = str1
    for p in plist:
        p1 = p.replace('-','')
        s = s.replace(p, p1)
    return s

def cluster_dict_by_value(d):
    #d is a list with tuple pair
    clusters = defaultdict(list)
    for key, val in d:
        clusters[val].append(key)
    return clusters
    
def cluster_dict_by_value2(d):
    #d is a defaultdict
    clusters = defaultdict(list)
    for key, val in d.items():
        clusters[','.join(val)].append(key)
    return clusters
    

if __name__ == "__main__":
    TFTA_Module(argv=sys.argv[1:])
