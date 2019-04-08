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
from .tfta import KinaseNotFoundException
from indra.sources.trips.processor import TripsProcessor
from collections import defaultdict
from bioagents import Bioagent

stmt_provenance_map = {'increase':'upregulates', 'decrease':'downregulates',
                 'regulate':'regulates'}
                 
stmt_provenance_map_target = {'increase':'upregulated', 'decrease':'downregulated',
                 'regulate':'regulated'}

stmt_type_map = {'increase':['IncreaseAmount'], 'decrease':['DecreaseAmount'],
                 'regulate':['RegulateAmount']}
                 
dbname_pmid_map = {'TRED':'17202159', 'ITFP':'18713790', 'ENCODE':'22955616',
                 'TRRUST':'26066708', 'Marbach2016':'26950747', 'Neph2012':'22959076'}


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
             'IS-GENE-ONTO', 'FIND-GENE-ONTO', 'FIND-KINASE-REGULATION',
             'FIND-TF-MIRNA', 'FIND-REGULATION', 'FIND-EVIDENCE', 'FIND-GENE-TISSUE',
             'IS-GENE-TISSUE', 'FIND-KINASE-PATHWAY','TEST-FUNCTION']
    #keep the genes from the most recent previous call, which are used to input 
    #find-gene-onto if there's no gene input 
    #gene_list = ['STAT3', 'JAK1', 'JAK2', 'ELK1', 'FOS', 'SMAD2', 'KDM4B']
    gene_list = []

    def __init__(self, **kwargs):
        #Instantiate a singleton TFTA agent
        self.tfta = TFTA()
        self.stmts_indra = dict()
        self.hgnc_info = dict()
    
        # Call the constructor of Bioagent
        super(TFTA_Module, self).__init__(**kwargs)
        
    def receive_tell(self, msg, content):
        #handle tell broadcast
        #now just do nothing here, but to avoid error message sending out
        pass

    def respond_is_regulation(self, content):
        """
        Response content to is-regulation request which includes:
        is-tf-target
        is-tf-target-tissue
        """
        tf_arg = content.gets('tf')
        target_arg = content.gets('target')
        tissue_arg = content.get('tissue')
        keyword_name = _get_keyword_name(content, descr='keyword')
        
        if all([tf_arg,target_arg,tissue_arg]):
            reply = self.respond_is_tf_target_tissue(content)
        elif all([tf_arg,target_arg, keyword_name]):
            if keyword_name in ['increase', 'decrease']:
                reply = self.respond_is_tf_target_literature(content)
            else:
                reply = self.respond_is_tf_target(content)
        elif all([tf_arg,target_arg]):
            reply = self.respond_is_tf_target(content)
        else:
            reply = make_failure('UNKNOW_TASK')
        return reply
    
    def respond_find_tf(self, content):
        """
        Response content to find-tf request which includes:
        find-target-tf, find-target-tf-tissue
        """
        target_arg = content.gets('target')
        tissue_arg = content.get('tissue')
        keyword_name = _get_keyword_name(content, descr='keyword')
        
        if all([target_arg,tissue_arg]):
            reply = self.respond_find_target_tfs_tissue(content)
        elif all([target_arg, keyword_name]):
            if keyword_name in ['increase', 'decrease']:
                reply = self.respond_find_tf_literature(content)
            else:
                reply = self.respond_find_target_tfs(content)
        elif target_arg:
            reply = self.respond_find_target_tfs(content)
        else:
            reply = make_failure('UNKNOW_TASK')
        return reply
    
    def respond_find_pathway(self, content):
        """
        Response content to find-pathway request, which includes:
        find-pathway-gene, find-pathway-db-gene, find-pathway-gene-keyword,
        find-pathway-keyword.
        """
        gene_arg = content.gets('gene')
        #pathway_arg = content.gets('pathway')
        db_arg = content.get('database')
        keyword_arg = content.get('keyword')
        #count_arg = content.get('count')
        if all([keyword_arg,gene_arg]):
            reply = self.respond_find_pathway_gene_keyword(content)
        elif all([gene_arg,db_arg]):
            reply = self.respond_find_pathway_db_gene(content)
        elif gene_arg:
            reply = self.respond_find_pathway_gene(content)
        elif keyword_arg:
            reply = self.respond_find_pathway_keyword(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
    
    def respond_find_target(self, content):
        """
        Response content to find-target request, which includes these cases:
        find-tf-target, find-tf-target-tissue. 
        """
        tf_arg = content.gets('tf')
        tissue_arg = content.get('tissue')
        keyword_name = _get_keyword_name(content, descr='keyword')
        
        if all([tf_arg,tissue_arg]):
            reply = self.respond_find_tf_targets_tissue(content)
        elif all([tf_arg, keyword_name]):
            if keyword_name in ['increase', 'decrease']:
                reply = self.respond_find_target_literature(content)
            else:
                reply = self.respond_find_tf_targets(content)
        elif tf_arg:
            reply = self.respond_find_tf_targets(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
    
    def respond_find_gene(self, content):
        """
        Response content to find-gene request, which includes:
        find-gene-pathway, find-tf-pathway, find-kinase-pathway
        """
        #pathway_arg = content.gets('pathway')
        subtype_name = _get_keyword_name(content, descr='subtype')
        if subtype_name == 'tf':
            reply = self.respond_find_tf_pathway(content)
        elif subtype_name == 'gene':
            reply = self.respond_find_gene_pathway(content)
        elif subtype_name == 'kinase':
            reply = self.respond_find_kinase_pathway(content)
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
        
    def respond_find_common_pathway_all(self, content):
        """
        Respond to find-common-pathway-genes, which covers:
        find-common-pathway-genes, find-common-pathway-genes-keyword,
        find-common-pathway-genes-database
        """
        target_arg = content.gets('target')
        keyword_arg = content.get('keyword')
        db_arg = content.get('database')
        if all([target_arg, keyword_arg]):
            reply = self.respond_find_common_pathway_genes_keyword(content)
        elif all([target_arg, db_arg]):
            reply = self.respond_find_common_pathway_genes_db(content)
        elif target_arg:
            reply = self.respond_find_common_pathway_genes(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
        
    def respond_find_regulation(self, content):
        """
        Respond to find-regulation
        """
        #target_arg = content.gets('target')
        agent_arg = content.get('regulator')
        source_arg = content.get('source')
        tissue_arg = content.get('tissue')
        if tissue_arg:
            temp_reply = self.respond_find_target_tfs_tissue(content)
            tf_str = temp_reply.get('tfs')
            if tf_str != 'NIL':
                reply = KQMLList.from_string(
                    '(SUCCESS :regulators (:tf-db ' + tf_str.to_string() + '))')
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        elif source_arg:
            reply = self.respond_find_regulation_source(content)
        elif agent_arg:
            reply = self.respond_find_regulation_agent(content)
        else:
            reply = self.respond_find_regulation_all(content)
        return reply
        
    def respond_is_gene_onto(self, content):
        """
        Respond to is-gene-onto
        """
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        gene_name = get_gene(content, descr='gene')
        gene_arg = content.gets('gene')
        if not gene_name:
            reply = self.wrap_family_message1(gene_arg, 'NO_GENE_NAME')
            return reply
            
        #check if keyword is protein or gene
        if keyword_name in ['gene', 'protein']:
            is_onto,is_ekb = self.is_protein_gene(gene_arg, keyword_name)
            if not is_ekb:
                reply = make_failure('ONLY_SUPPORT_EKB')
                return reply
        else:
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
        regulator_arg = content.gets('regulator')
        gene_arg = content.gets('gene')
        if regulator_arg:
            reply = self.respond_find_gene_onto_regulator(content)
        elif gene_arg:
            reply = self.respond_find_gene_onto_gene(content)
        else:
            reply = make_failure('UNKNOWN_TASK')
        return reply
        
    def respond_find_gene_onto_gene(self, content):
        """
        Respond to find-gene-onto
        """
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        gene_names = get_of_those_list(content, descr='gene')
        if not gene_names:
            gene_arg = content.gets('gene')
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            results = self.tfta.find_gene_onto(keyword_name, gene_names)
        except GONotFoundException:
            reply = make_failure('GO_NOT_FOUND')
            return reply
        if len(results):
            gene_list_str = self.wrap_message(':genes', results)
            reply = KQMLList.from_string(
                '(SUCCESS ' + gene_list_str + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :genes NIL)')
        return reply
        
    def respond_find_gene_onto_regulator(self, content):
        """
        Respond to find-gene-onto
        """
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
            
        regulator_names = get_of_those_list(content, descr='regulator')
        if not regulator_names:
            regulator_arg = content.gets('regulator')
            reply = _wrap_family_message(regulator_arg, 'NO_REGULATOR_NAME')
            return reply
        #get genes regulated by regulators in regulator_names
        try:
            gene_names,dbname = self.tfta.find_targets(regulator_names)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :genes NIL)')
            return reply
        if len(gene_names):
            try:
                results = self.tfta.find_gene_onto(keyword_name, gene_names)
            except GONotFoundException:
                reply = make_failure('GO_NOT_FOUND')
                return reply
            if len(results):
                gene_list_str = self.wrap_message(':genes', results)
                reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
            else:
                reply = KQMLList.from_string('(SUCCESS :genes NIL)')
        else:
            reply = KQMLList.from_string('(SUCCESS :genes NIL)')
        return reply
        
    def respond_is_tf_target(self, content):
        """
        Response content to is-tf-target request.
        """
        literature = False
        tf_name = get_gene(content, descr='tf')
        if not tf_name:
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
        
        target_name = get_gene(content, descr='target')
        if not target_name:
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        #check if it exists in literature
        term_tuple = (tf_name, target_name, 'REGULATE')
        if term_tuple not in self.stmts_indra:
            stmts,success = self.tfta.find_statement_indraDB(subj=tf_name, obj=target_name, stmt_types=['RegulateAmount'])
            if success:
                self.stmts_indra[term_tuple] = stmts
        else:
            stmts = self.stmts_indra[term_tuple]
        #provenance support
        self.send_background_support(stmts, tf_name, target_name, 'regulate')
        if len(stmts):
            literature = True
        #check the db
        try:
            is_target,dbname = self.tfta.Is_tf_target(tf_name, target_name)
            #provenance support
            self.send_background_support_db(tf_name, target_name, dbname)
        except Exception as e:
            #provenance support
            self.send_background_support_db(tf_name, target_name, '')
            if not literature:
                reply = KQMLList.from_string('(SUCCESS :result FALSE :db FALSE :literature FALSE)')
                return reply
        
        result = is_target or literature
        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if result else 'FALSE'
        reply.set('result', is_target_str)
        reply.set('db', 'TRUE' if is_target else 'FALSE')
        reply.set('literature', 'TRUE' if literature else 'FALSE')
        return reply
        
    def respond_is_tf_target_literature(self, content):
        """
        Response content to is-tf-target request with regulating direction
        """
        literature = False
        tf_name = get_gene(content, descr='tf')
        if not tf_name:
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
            
        target_name = get_gene(content, descr='target')
        if not target_name:
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
            
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        try:
            stmt_types = stmt_type_map[keyword_name]
        except KeyError as e:
            reply = make_failure('INVALID_KEYWORD')
            return reply
            
        term_tuple = (tf_name, target_name, keyword_name)
        if term_tuple not in self.stmts_indra:
            stmts,success = self.tfta.find_statement_indraDB(subj=tf_name, obj=target_name, stmt_types=stmt_types)
            if success:
                self.stmts_indra[term_tuple] = stmts
        else:
            stmts = self.stmts_indra[term_tuple]
        #provenance support
        self.send_background_support(stmts, tf_name, target_name, keyword_name)
        if len(stmts):
            literature = True
        #set is_target to False since it doesn't support 'increase' or 'decrease' query
        is_target = False
        result = literature
        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if result else 'FALSE'
        reply.set('result', is_target_str)
        reply.set('db', 'FALSE')
        reply.set('literature', 'TRUE' if literature else 'FALSE')
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
        tf_name = get_gene(content, descr='tf')
        if not tf_name:
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
            
        target_name = get_gene(content, descr='target')
        if not target_name:
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        #get tissue name
        tissue_name = _get_tissue_name(content)
        if not tissue_name:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        try:
            is_target = self.tfta.Is_tf_target_tissue(tf_name, target_name, tissue_name)
        except Exception as e:
            reply = KQMLList.from_string('(SUCCESS :result FALSE)')
            return reply
        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if is_target else 'FALSE'
        reply.set('result', is_target_str)
        return reply

    def respond_find_tf_targets(self, content):
        """
        Response content to find-tf-target request
        Covered by find-target
        For a tf list, reply with the targets found
        """
        tf_names = get_of_those_list(content, descr='tf')
        if not len(tf_names):
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply

        #tf_str = ', '.join(tf_names[:-1]) + ' and ' + tf_names[-1]
        #consider an optional parameter for subsequent query
        of_targets_names = get_of_those_list(content)
        
        #consider an optional parameter to only return given type of genes
        target_type_set = self.get_target_type_set(content)
        
        try:
            target_names,dbname = self.tfta.find_targets(tf_names)
        except TFNotFoundException:
            #provenance support
            self.send_background_support_db(tf_names, [], '')
            
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
            return reply
        #check if it's a subsequent query
        if of_targets_names:
            target_names = list(set(of_targets_names) & set(target_names))
        
        #check if it requires returning a type of genes
        if target_type_set:
            target_names = target_type_set & set(target_names)
        
        if len(target_names):
            gene_list_str = self.wrap_message(':targets', target_names)
            reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
            
        self.gene_list = target_names
        
        #provenance support
        self.send_background_support_db(tf_names, target_names, dbname, find_target=True)
        
        return reply
        
    def respond_find_target_literature(self, content):
        """
        Respond to find-target request with information in literatures
        """
        tf_names = get_of_those_list(content, descr='tf')
        if not len(tf_names):
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        #consider an optional parameter for sequencing query
        of_those_names = get_of_those_list(content)
        
        #consider an optional parameter to only return given types of genes
        target_type = _get_keyword_name(content, descr='target-type')
    
        try:
            stmt_types = stmt_type_map[keyword_name]
        except KeyError as e:
            reply = make_failure('INVALID_KEYWORD')
            return reply
        lit_messages = self.get_target_indra(tf_names, stmt_types, keyword_name, 
                                             of_those=of_those_names, target_type=target_type)
        if len(lit_messages):
            reply = KQMLList.from_string(
                    '(SUCCESS ' + lit_messages + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply

    def respond_find_tf_targets_tissue(self, content):
        """
        Response content to find-tf-target-tissue request
        For tf list, reply with the targets found within given tissue
        """
        tf_names = get_of_those_list(content, descr='tf')
        if not len(tf_names):
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
            
        #get tissue name
        tissue_name = _get_tissue_name(content)
        if not tissue_name:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
            
        #consider an optional parameter for sequencing query
        of_targets_names = get_of_those_list(content)
        
        #consider an optional parameter to only return given types of genes
        target_type_set = self.get_target_type_set(content)
        
        try:
            target_names = self.tfta.find_targets_tissue(tf_names, tissue_name)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
            return reply
            
        #check if it's a sequencing query
        if of_targets_names:
            target_names = list(set(of_targets_names) & set(target_names))
        #check if it requires returning a type of genes
        if target_type_set:
            target_names = target_type_set & set(target_names)
            
        if len(target_names):
            gene_list_str = self.wrap_message(':targets', target_names)
            reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        self.gene_list = target_names
        return reply

    def respond_find_target_tfs(self, content):
        """
        Response content to find-target-tf request
        For a target list, reply the tfs found
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
            
        #consider an optional parameter for subsequent query
        of_tfs_names = get_of_those_list(content)
        
        try:
            tf_names,dbname = self.tfta.find_tfs(target_names)
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
            #provenance support
            self.send_background_support_db([], target_names, '')
            return reply
        #check if it's a subsequent query
        if of_tfs_names:
            tf_names = list(set(of_tfs_names) & set(tf_names))
        if len(tf_names):
            gene_list_str = self.wrap_message(':tfs', tf_names)
            reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        #self.gene_list = tf_names
        #provenance support
        self.send_background_support_db(tf_names, target_names, dbname, find_tf=True)
        return reply
        
    def respond_find_tf_literature(self, content):
        """
        Respond to find-tf request with information in literatures
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        try:
            stmt_types = stmt_type_map[keyword_name]
        except KeyError as e:
            reply = make_failure('INVALID_KEYWORD')
            return reply
        lit_messages = self.get_tf_indra(target_names, stmt_types, keyword_name)
        if len(lit_messages):
            reply = KQMLList.from_string(
                    '(SUCCESS ' + lit_messages + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        return reply

    def respond_find_target_tfs_tissue(self, content):
        """
        Response content to find-target-tf-tissue request
        For a target list, reply the tfs found within a given tissue
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        #get tissue name
        tissue_name = _get_tissue_name(content)
        if not tissue_name:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
        #consider an optional parameter for subsequent query
        of_tfs_names = get_of_those_list(content)
        
        try:
            tf_names = self.tfta.find_tfs_tissue(target_names, tissue_name)
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
            return reply
        #check if it's a sequencing query
        if of_tfs_names:
            tf_names = list(set(of_tfs_names) & set(tf_names))         
        if len(tf_names):
            gene_list_str = self.wrap_message(':tfs', tf_names)
            reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        self.gene_list = tf_names
        return reply

    def respond_find_pathway_gene(self,content):
        """
        Response content to find_pathway_gene request
        For a given gene list, reply the related pathways information
        """
        gene_names = get_of_those_list(content, descr='gene')
        if not len(gene_names):
            gene_arg = content.gets('gene')
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName, dblink = \
                self.tfta.find_pathways_from_genelist(gene_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink)
        return reply

    def respond_find_pathway_gene_keyword(self,content):
        """
        Response content to find_pathway_gene_name request
        For a given gene list and keyword, reply the related pathways information
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name or keyword_name in ['pathway', 'signaling pathway']:
            reply = make_failure('NO_KEYWORD')
            return reply
            
        gene_names = get_of_those_list(content, descr='gene')
        if not len(gene_names):
            gene_arg = content.gets('gene')
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName, dblink = \
                self.tfta.find_pathway_gene_keyword(gene_names, keyword_name)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=[keyword_name])
        return reply

    def respond_find_pathway_db_gene(self, content):
        """Response content to FIND-PATHWAY-DB-GENE request
        For a given gene list and certain db source, reply the related
        pathways information"""
        db_name = _get_keyword_name(content, descr='database', low_case=False)
        if not db_name:
            reply = make_failure('NO_DB_NAME')
            return reply
            
        gene_names = get_of_those_list(content, descr='gene')
        if not gene_names:
            gene_arg = content.gets('gene')
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName,dblink = self.tfta.find_pathways_from_dbsource_geneName(db_name,gene_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink)
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
        
        #consider another optional parameter for subsequent query
        of_gene_names = get_of_those_list(content, descr='of-those')
        
        try:
            pathwayName,tflist,dblink = \
                self.tfta.find_tf_pathway(pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_genelist_message(pathwayName, dblink, tflist, pathway_names=pathway_names,
                                               gene_descr=':tfs', of_gene_names=of_gene_names)
        return reply
        
    def respond_find_kinase_pathway(self, content):
        """
        Response content to FIND_KINASE_PATHWAY request
        For a given pathway name, reply the kinases within the pathway
        """
        pathway_arg = content.gets('pathway')
        pathway_names = _get_pathway_name(pathway_arg)
        pathway_names = trim_word(pathway_names, 'pathway')
        if not len(pathway_names):
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
            
        #consider an optional parameter for subsequent query
        of_those_names = get_of_those_list(content)
        
        try:
            pathwayName, kinaselist, dblink = \
                self.tfta.find_kinase_pathway(pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        
        reply = _wrap_pathway_genelist_message(pathwayName, dblink, kinaselist, pathway_names=pathway_names,
                                               gene_descr=':kinase', of_gene_names=of_those_names)
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
            pathwayName,genelist,plink = self.tfta.find_genes_from_pathwayName(pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        #consider an optional parameter for subsequent query
        of_gene_names = get_of_those_list(content)
        
        reply = _wrap_pathway_genelist_message(pathwayName, plink, genelist, pathway_names=pathway_names,
                                       gene_descr=':genes', of_gene_names=of_gene_names)
        return reply

    def respond_find_pathway_keyword(self, content):
        """
        Response content to FIND-PATHWAY-KEYWORD request
        For a given keyword, reply the pathways involving the keyword
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name or keyword_name in ['pathway', 'signaling pathway']:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        else:
            keyword_name = trim_word([keyword_name], 'pathway')
        
        try:
            pathwayName, dblink = \
                self.tfta.find_pathway_keyword(keyword_name[0])
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=keyword_name)
        return reply

    def respond_find_tf_keyword(self, content):
        """
        Response content to FIND-TF-KEYWORD request
        For a given keyword, reply the tfs within the pathways
        involving the keyword
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        
        try:
            pathwayName, tflist, dblink = self.tfta.find_tf_keyword(keyword_name)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        
        #consider an optional parameter for subsequent query
        of_those_names = get_of_those_list(content)
        
        reply = _wrap_pathway_genelist_message(pathwayName, dblink, tflist, pathway_names=[keyword_name],
                                               gene_descr=':tfs', of_gene_names=of_those_names)
        return reply

    def respond_find_common_tfs_genes(self, content):
        """
        Response content to FIND-COMMON-TF-GENES request
        For a given target list, reply the tfs regulating these genes
        and the regulated targets by each TF
        """
        target_names = get_of_those_list(content, descr='target')
        if not target_names:
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        #consider another parameter for subsequent query
        of_those_names = get_of_those_list(content, descr='of-those')
        
        try:
            tf_targets = self.tfta.find_common_tfs(target_names)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
            return reply
        #cluster the tfs according to the targets
        tf_clustered = cluster_dict_by_value2(tf_targets)
        tf_list_str = ''
        targets = list(tf_clustered.keys())
        for tg in targets:
            tf_list = ''
            if of_those_names:
                tfs = set(of_those_names) & set(tf_clustered[tg])
            else:
                tfs = tf_clustered[tg]
            if tfs:
                target_list = tg.split(',')
                target_str = ''
                for t in target_list:
                    target_str += '(:name %s) ' % t
                target_str = ':target-list (' + target_str + ')'
                for tf in tfs:
                    tf_list += '(:name %s) ' % tf
                tf_list = ':tf-list (' + tf_list + ') '
                tf_list += target_str
                tf_list_str += '(' + tf_list + ') '
        if tf_list_str:
            reply = KQMLList.from_string(
                '(SUCCESS :tfs (' + tf_list_str + '))')
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        return reply
        
    def respond_find_common_pathway_genes(self, content):
        """
        response content to FIND-COMMON-PATHWAY-GENES request
        """
        gene_arg = content.gets('target')
        if not gene_arg:
            gene_arg = content.gets('gene')
        gene_names = get_of_those_list(content, descr='target')
        if not gene_names:
            gene_names = get_of_those_list(content, descr='gene')
        if not len(gene_names):
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName, dblink, genes = self.tfta.find_common_pathway_genes(gene_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = _wrap_pathway_genelist_message(pathwayName, dblink, genes, gene_descr=':gene-list')
        return reply

    def respond_find_common_pathway_genes_keyword(self, content):
        """
        Response content to FIND-COMMON-PATHWAY-GENES-KEYWORD request
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name or keyword_name in ['pathway', 'signaling pathway']:
            reply = make_failure('NO_KEYWORD')
            return reply
        else:
            keyword_name = trim_word([keyword_name], 'pathway')
            
        gene_arg = content.gets('target')
        if not gene_arg:
            gene_arg = content.gets('gene')
        gene_names = get_of_those_list(content, descr='target')
        if not gene_names:
            gene_names = get_of_those_list(content, descr='gene')
        if not len(gene_names):
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName,dblink,genes = \
               self.tfta.find_common_pathway_genes_keyword(gene_names, keyword_name[0])
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = _wrap_pathway_genelist_message(pathwayName, dblink, genes, pathway_names=keyword_name,
                                               gene_descr=':gene-list')
        return reply
        
    def respond_find_common_pathway_genes_db(self, content):
        """
        Response content to FIND-COMMON-PATHWAY-GENES-DB request
        """
        db_name = _get_keyword_name(content, descr='database', low_case=False)
        if not db_name:
            reply = make_failure('NO_DB_NAME')
            return reply
            
        gene_arg = content.gets('target')
        if not gene_arg:
            gene_arg = content.gets('gene')
        gene_names = get_of_those_list(content, descr='target')
        if not gene_names:
            gene_names = get_of_those_list(content, descr='gene')
        if not len(gene_names):
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName,dblink,genes = \
               self.tfta.find_common_pathway_genes_db(gene_names, db_name)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = _wrap_pathway_genelist_message(pathwayName, dblink, genes, gene_descr=':gene-list')
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
            
        gene_names = get_of_those_list(content, descr='gene')
        if not len(gene_names):
            gene_arg = content.gets('gene')
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            pathwayName, dblink = \
                self.tfta.Is_pathway_gene(pathway_names, gene_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=pathway_names)
        return reply

    def respond_find_genes_go_tf(self, content):
        """
        Respond content to FIND-GENE-GO-TF request
        The inputs here are keyword from GO names and tf list
        """
        keyword_name = _get_keyword_name(content)
        if not keyword_name:
            reply = make_failure('NO_GO_NAME')
            return reply
            
        tf_names = get_of_those_list(content, descr='tf')
        if not tf_names:
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
            
        try:
            go_ids,go_types,go_names,go_genes = \
                self.tfta.find_genes_GO_tf(keyword_name, tf_names)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :go-terms NIL)')
            return reply
        except GONotFoundException:
            reply = KQMLList.from_string('(SUCCESS :go-terms NIL)')
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
        goid = _get_keyword_name(content, descr='goid', low_case=False)
        if not goid:
            reply = make_failure('NO_GO_ID')
            return reply
        
        tf_names = get_of_those_list(content, descr='tf')
        if not tf_names:
            tf_arg = content.gets('tf')
            reply = _wrap_family_message(tf_arg, 'NO_TF_NAME')
            return reply
        
        try:
            go_ids,go_types,go_names,go_genes = \
                self.tfta.find_genes_GO_tf2(goid, tf_names)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :go-terms NIL)')
            return reply
        except GONotFoundException:
            reply = KQMLList.from_string('(SUCCESS :go-terms NIL)')
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
        miRNA_name = list(_get_miRNA_name(miRNA_arg).keys())
        if not len(miRNA_name):
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        
        target_name = get_gene(content, descr='target')
        if not len(target_name):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
            
        try:
            is_target = self.tfta.Is_miRNA_target(miRNA_name[0], target_name)
        except miRNANotFoundException:
            reply = KQMLList.from_string('(SUCCESS :is-miRNA-target FALSE)')
            return reply
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :is-miRNA-target FALSE)')
            return reply            
        reply = KQMLList('SUCCESS')
        is_target_str = 'TRUE' if is_target else 'FALSE'
        reply.set('is-miRNA-target', is_target_str)
        return reply
        
    def respond_is_miRNA_target2(self, content):
        """
        Respond to IS-MIRNA-TARGET request
        """
        miRNA_arg = content.gets('miRNA')
        miRNA_name = _get_miRNA_name(miRNA_arg)
        if not len(miRNA_name):
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        
        target_name = get_gene(content, descr='target')
        if not len(target_name):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        is_target,expr,supt,pmid,miRNA_mis = self.tfta.Is_miRNA_target2(miRNA_name, target_name)
        
        #provenance support
        self.send_background_support_mirna(list(miRNA_name.keys())[0], target_name, expr, supt, pmid)
        
        #respond to BA
        #check if it's necessary for user clarification
        if len(miRNA_mis):
            reply = self.get_reply_mirna_clarification(miRNA_mis)
            return reply
        else:
            reply = KQMLList('SUCCESS')
            is_target_str = 'TRUE' if is_target else 'FALSE'
            reply.set('is-miRNA-target', is_target_str)
        return reply
        
    def respond_find_miRNA_target(self, content):
        """
        Respond to FIND-MIRNA-TARGET request
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        ##consider another parameter for subsequent query
        of_those_names = get_of_those_mirna(content, descr='of-those')
        
        try:
            miRNA_names = self.tfta.find_miRNA_target(target_names)
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)') 
            return reply
            
        if of_those_names:
            res_db = set(','.join(miRNA_names).upper().split(','))
            res_of = set(','.join(of_those_names).upper().split(','))
            miRNA_names = res_db & res_of
            
        if len(miRNA_names):
            miRNA_list_str = ''
            for m in miRNA_names:
                miRNA_list_str += '(:name %s ) ' % m
            reply = KQMLList.from_string(
                '(SUCCESS :miRNAs (' + miRNA_list_str + '))')
        else:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)')
        return reply
        
    def respond_find_target_miRNA(self, content):
        """
        Respond to FIND-TARGET-MIRNA request
        """
        strength_arg = content.get('strength')
        if strength_arg:
            reply = self.respond_find_target_miRNA_strength(content)
        else:
            reply = self.respond_find_target_miRNA_all(content)
        return reply

    def respond_find_target_miRNA_strength(self, content):
        """
        Respond to FIND-TARGET-MIRNA request with strength parameter
        """
        try:
            miRNA_arg = content.gets('miRNA')
            miRNA_names = _get_miRNA_name(miRNA_arg)
        except Exception as e:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        if not len(miRNA_names):
            reply = make_failure('NO_MIRNA_NAME')
            return reply
            
        strength_name = _get_keyword_name(content, descr='strength')
        if not strength_name:
            reply = make_failure('NO_STRENGTH_NAME')
            return reply
            
        #consider another parameter for subsequent query
        of_those_names = get_of_those_list(content)
        
        #consider an optional parameter to only return given type of genes
        target_type_set = self.get_target_type_set(content)
        
        try:
            target_names,miRNA_mis = self.tfta.find_target_miRNA_strength(miRNA_names, strength_name)
        except miRNANotFoundException:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
            return reply
        #check if it's necessary for user clarification
        if len(miRNA_mis):
            reply = self.get_reply_mirna_clarification(miRNA_mis)
            return reply
        else:
            if of_those_names:
                target_names = set(of_those_names) & set(target_names)
            #check if it requires returning a type of genes
            if target_type_set:
                target_names = target_type_set & set(target_names)
            
            if len(target_names):
                gene_list_str = self.wrap_message(':targets', target_names)
                reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        #self.gene_list = target_names  
        return reply
        
    def respond_find_target_miRNA_all(self, content):
        """
        Respond to FIND-TARGET-MIRNA request without strength parameter
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
            
        #consider another parameter for subsequent query
        of_those_names = get_of_those_list(content)
        
        #consider an optional parameter to only return given type of genes
        target_type_set = self.get_target_type_set(content)
        
        target_names,miRNA_mis = self.tfta.find_target_miRNA(miRNA_names)
        #check if it's necessary for user clarification
        if len(miRNA_mis):
            reply = self.get_reply_mirna_clarification(miRNA_mis)
            return reply
        else:
            if of_those_names:
                target_names = set(of_those_names) & set(target_names)
                
            #check if it requires returning a type of genes
            if target_type_set:
                target_names = target_type_set & set(target_names)
            
            if len(target_names):
                gene_list_str = self.wrap_message(':targets', target_names)
                reply = KQMLList.from_string('(SUCCESS ' + gene_list_str + ')')
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        #self.gene_list = target_names  
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
        
        target_name = get_gene(content, descr='target')
        if not len(target_name):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
                   
        try:
            experiments,supportType,pmlink,miRNA_mis = \
                self.tfta.find_evidence_miRNA_target(miRNA_name, target_name)
        except Exception as e:
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)')
            return reply
        if miRNA_mis:
            reply = self.get_reply_mirna_clarification(miRNA_mis)
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
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)')   
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
        
        #consider another parameter for subsequent query
        of_those_names = get_of_those_list(content)
        
        targets,counts,mrna,miRNA_mis = self.tfta.find_gene_count_miRNA(miRNA_names)
        if not len(targets):
            #clarification
            if miRNA_mis:
                reply = self.get_reply_mirna_clarification(miRNA_mis)
                return reply
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
                return reply
                
        target_str = ''
        if of_those_names:
            ftargets = set(of_those_names) & set(targets)
        else:
            ftargets = targets
        if ftargets:
            for t in ftargets:
                ct = counts[t]
                ms = mrna[t]
                m_str = ''
                for m in ms:
                    m_str += '(:name %s)' % m
                m_str = '(' + m_str + ')'
                target_str += '(:name %s :count %d :miRNA %s) ' % (t, ct, m_str)
            reply = KQMLList.from_string('(SUCCESS :targets (' + target_str + '))')
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        #self.gene_list = targets
        return reply
             
    def respond_find_miRNA_count_target(self, content):
        """
        Respond to FIND-MIRNA-COUNT-GENE request
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
            
        try:
            mirnas,counts,genes = self.tfta.find_miRNA_count_gene(target_names)
        except miRNANotFoundException:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)')
            return reply
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)')
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
        db_name = _get_keyword_name(content, descr='database', low_case=False)
        if not db_name:
            reply = make_failure('NO_DB_NAME')
            return reply
        
        pathway_names = _get_keyword_name(content, hyphen=True)
        if not pathway_names:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        else:
            pathway_names = trim_word([pathway_names], 'pathway')
        
        try:
            pathwayName,dblink = \
                    self.tfta.find_pathway_db_keyword(db_name, pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=pathway_names)
        return reply

    def respond_find_tissue_gene(self, content):
        """
        Response to FIND-TISSUE-GENE task
        """
        gene_arg = content.gets('gene')
        #without gene parameter, then return all the tissues
        if not gene_arg:
            tissue_str = ''
            for ts in self.tissue_list:
                tissue_str += '(:name %s) ' % ts
            reply = KQMLList.from_string(
                    '(SUCCESS :tissue (' + tissue_str + '))')
            return reply
        gene_name = get_gene(content, descr='gene')
        if not len(gene_name):
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        try:
            tissues = self.tfta.find_tissue_gene(gene_name)
        except TissueNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tissue NIL)')
            return reply    
        tissue_str = ''
        for ts in tissues:
            tissue_str += '(:name %s) ' % ts
        reply = KQMLList.from_string(
            '(SUCCESS :tissue (' + tissue_str + '))')
        return reply
        
    def respond_find_kinase_target(self, content):
        """
        Response to find-kinase-regulation with only target parameter
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
            
        try:
            kinases = self.tfta.find_kinase_target(target_names)
        except KinaseNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :kinase NIL)')
            return reply
        kinase_str = self.wrap_message(':kinase ', kinases)
        reply = KQMLList.from_string('(SUCCESS ' + kinase_str + ')')
        return reply
        
    def respond_find_kinase_target_keyword(self, content):
        """
        Response to find-kinase-regulation with target and keyword parameter
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        try:
            kinases = self.tfta.find_kinase_target_keyword(target_names, keyword_name)
        except KinaseNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :kinase NIL)')
            return reply
        kinase_str = self.wrap_message(':kinase ', kinases)
        reply = KQMLList.from_string('(SUCCESS ' + kinase_str + ')')
        return reply
        
    def respond_find_kinase_regulation(self, content):
        """
        Response to FIND-KINASE-REGULATION
        """
        target_arg = content.gets('target')
        keyword_arg = content.get('keyword')
        if all([target_arg, keyword_arg]):
            reply = self.respond_find_kinase_target_keyword(content)
        elif target_arg:
            reply = self.respond_find_kinase_target(content)
        else:
            reply = make_failure('UNKNOW_TASK')
        return reply
        
    def respond_find_tf_miRNA(self, content):
        """
        Response to FIND-TF-MIRNA
        """
        try:
            miRNA_arg = content.gets('miRNA')
            miRNA_names = _get_miRNA_name(miRNA_arg)
        except Exception as e:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        if not len(miRNA_names):
            reply = make_failure('NO_MIRNA_NAME')
            return reply 
            
        #consider an additional parameter for subsequent query
        of_those_names = get_of_those_list(content, descr='of-those')
        
        try:
            tf_names,miRNA_mis = self.tfta.find_tf_miRNA(miRNA_names)
            #check if it's necessary for user clarification
            if len(miRNA_mis):
                reply = self.get_reply_mirna_clarification(miRNA_mis)
                return reply
            else:
                if of_those_names:
                    tf_names = set(of_those_names) & set(tf_names)
                tf_list_str = self.wrap_message(':tfs ', tf_names)
                reply = KQMLList.from_string('(SUCCESS ' + tf_list_str + ')')
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
            return reply
        #self.gene_list = tf_names  
        return reply
    
    def respond_find_regulation_source(self, content):
        """
        Respond to find-regulation
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        target_names = list(set(target_names))
        
        keyword_name = _get_keyword_name(content)
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
            
        source_name = _get_keyword_name(content, descr='source')
        if not source_name:
            reply = make_failure('INVALID_SOURCE')
            return reply
            
        #consider an optional parameter for subsequent query
        of_those_names = set(get_of_those_list(content))
            
        if source_name == 'literature':
            #literature result
            try:
                stmt_types = stmt_type_map[keyword_name]
            except KeyError as e:
                reply = make_failure('INVALID_KEYWORD')
                return reply
            if of_those_names:
                lit_messages = self.get_regulator_indra_those(target_names, stmt_types, keyword_name,
                                                        of_those=of_those_names)
            else:
                lit_messages = self.get_regulator_indra(target_names, stmt_types, keyword_name)
            if len(lit_messages):
                reply = KQMLList.from_string(
                        '(SUCCESS :regulators (' + lit_messages + '))')
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        elif source_name == 'geo rnai':
            #kinase regualtion
            kin_messages = self.get_kinase_regulation(target_names, keyword_name, of_those=of_those_names)
            
            if len(kin_messages):
                reply = KQMLList.from_string(
                        '(SUCCESS :regulators (' + kin_messages + '))')
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        else:
            reply = make_failure('INVALID_SOURCE')
        return reply
        
    def respond_find_regulation_agent(self, content):
        """
        Response to find-regulation request
        For example: Which kinases regulate the cfos gene?
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        target_names = list(set(target_names))
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        agent_name = _get_keyword_name(content, descr='regulator')
        if not agent_name:
            reply = make_failure('NO_REGULATOR_NAME')
            return reply
            
        #consider an optional parameter for subsequent query
        of_those_names = set(get_of_those_list(content))
            
        if agent_name == 'kinase':
            #kinase regualtion
            kin_messages = self.get_kinase_regulation(target_names, keyword_name, of_those=of_those_names)
            if len(kin_messages):
                reply = KQMLList.from_string(
                        '(SUCCESS :regulators (' + kin_messages + '))')
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        elif agent_name in ['transcription factor', 'tf']:
            temp_reply = self.respond_find_target_tfs(content)
            tf_str = temp_reply.get('tfs')
            if tf_str != 'NIL':
                reply = KQMLList.from_string(
                    '(SUCCESS :regulators (:tf-db ' + tf_str.to_string() + '))')
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        else:
            reply = make_failure('INVALID_REGULATOR')
        return reply
        
    def respond_find_regulation_all(self, content):
        """
        Response to find-regulation request
        For example: what regulate MYC?
        """
        target_names = get_of_those_list(content, descr='target')
        if not len(target_names):
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        target_names = list(set(target_names))
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        #consider an optional parameter for subsequent query
        of_those_names = set(get_of_those_list(content))
        
        #literature result
        try:
            stmt_types = stmt_type_map[keyword_name.lower()]
        except KeyError as e:
            reply = make_failure('INVALID_KEYWORD')
            return reply
            
        if of_those_names:
            lit_messages = self.get_regulator_indra_those(target_names, stmt_types, keyword_name,
                                                          of_those=of_those_names)
        else:
            lit_messages = self.get_regulator_indra(target_names, stmt_types, keyword_name)
        
        #kinase regulation
        kin_messages = self.get_kinase_regulation(target_names, keyword_name, of_those=of_those_names)
        
        #db result
        #only considering the regulation case
        if keyword_name == 'regulate':
            try:
                tf_names,dbname = self.tfta.find_tfs(target_names)
                if of_those_names:
                    tf_names = list(set(tf_names).intersection(of_those_names))
                #provenance support
                self.send_background_support_db(tf_names, target_names, dbname, find_tf=True)
            except TargetNotFoundException:
                #provenance support
                self.send_background_support_db([], target_names, '')
                
                messages = _combine_messages([kin_messages, lit_messages])
                if len(messages):
                    reply = KQMLList.from_string(
                        '(SUCCESS :regulators (' + messages + '))')
                else:
                    reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
                return reply
        else:
            tf_names = []
        if len(tf_names):
            tf_messages = self.wrap_message(':tf-db', tf_names)
            messages = _combine_messages([tf_messages, kin_messages, lit_messages])
            reply = KQMLList.from_string(
                '(SUCCESS :regulators (' + messages + '))')
        else:
            messages = _combine_messages([kin_messages, lit_messages])
            if len(messages):
                reply = KQMLList.from_string(
                        '(SUCCESS :regulators (' + messages + '))')
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        #self.gene_list = tf_names
        return reply
        
    def respond_find_evidence(self, content):
        """
        response to find-evidence request
        """
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        if keyword_name == 'bind':
            reply = self.respond_find_evidence_tfdb(content)
        else:
            reply = self.respond_find_evidence_literature(content)
        return reply
    
    def respond_find_evidence_literature(self, content):
        """
        response to find-evidence request from literature
        for example: show me evidence that IL6 increases the amount of SOCS1.
        Only consider one-one case
        """
        regulator_name = get_gene(content, descr='regulator')
        if not regulator_name:
            regulator_arg = content.gets('regulator')
            reply = _wrap_family_message(regulator_arg, 'NO_REGULATOR_NAME')
            return reply
        
        target_name = get_gene(content, descr='target')
        if not target_name:
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        try:
            stmt_types = stmt_type_map[keyword_name]
        except KeyError as e:
            reply = make_failure('INVALID_KEYWORD')
            return reply
        
        evi_message = ''
        db_message = ''
        term_tuple = (regulator_name, target_name, keyword_name)
        if term_tuple not in self.stmts_indra:
            stmts,success = self.tfta.find_statement_indraDB(subj=regulator_name, obj=target_name, stmt_types=stmt_types)
            if success:
                self.stmts_indra[term_tuple] = stmts
        else:
            stmts = self.stmts_indra[term_tuple]
        self.send_background_support(stmts, regulator_name, target_name, keyword_name)
        if len(stmts):
            evidences = self.tfta.find_evidence_indra(stmts)
            if len(evidences):
                evi_message = _wrap_evidence_message(':literature', evidences)
        if keyword_name == 'regulate':
            db_names = self.tfta.find_evidence_dbname(regulator_name, target_name)
            if len(db_names):
                for db in db_names:
                    try:
                        pmid = dbname_pmid_map[db]
                    except KeyError as e:
                        pmid = ''
                    db_message += '(:name %s :pmid %s) ' % (db, pmid)
                db_message = ':tf-db (' + db_message + ') '
        message = _combine_messages([db_message, evi_message])
        if message:
            reply = KQMLList.from_string(
                             '(SUCCESS :evidence (' + message + '))')
        else:
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)')
        return reply
    
    def respond_find_evidence_tfdb(self, content):
        """
        Respond to find-evidence request from tfdb
        """
        regulator_name = get_gene(content, descr='regulator')
        if not regulator_name:
            regulator_arg = content.gets('regulator')
            reply = _wrap_family_message(regulator_arg, 'NO_REGULATOR_NAME')
            return reply
        
        target_name = get_gene(content, descr='target')
        if not target_name:
            target_arg = content.gets('target')
            reply = _wrap_family_message(target_arg, 'NO_TARGET_NAME')
            return reply
            
        db_names = self.tfta.find_evidence_dbname(regulator_name, target_name)
        if len(db_names):
            db_str = ''
            for db in db_names:
                try:
                    pmid = dbname_pmid_map[db]
                except KeyError as e:
                    pmid = ''
                db_str += '(:name %s :pmid %s) ' % (db, pmid)
            reply = KQMLList.from_string(
                    '(SUCCESS :evidence (:tf-db (' + db_str + ')))')
        else:
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)')
        return reply
            
    def respond_find_gene_tissue(self, content):
        """
        respond to find-gene-tissue request
        for example: which of these are expressed in liver?(subsequent query)
        what genes are expressed in liver?
        """
        #get tissue name
        tissue_name = _get_tissue_name(content)
        if not tissue_name:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
            
        #optional keyword
        keyword_name = _get_keyword_name(content, descr='keyword')
        if len(keyword_name):
            try:
                gene_list = self.tfta.find_gene_tissue_exclusive(tissue_name)
            except TissueNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :genes NIL)')
                return reply
        else:
            try:
                gene_list = self.tfta.find_gene_tissue(tissue_name)
            except TissueNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :genes NIL)')
                return reply
        
        #check if it's a subsequent query with optional parameter
        ogenes = get_of_those_list(content, descr='gene')
        if len(ogenes):
            results = list(set(gene_list) & set(ogenes))
        else:
            results = gene_list
        if len(results):
            results.sort()
            gene_str = self.wrap_message(':genes ', results)
            reply = KQMLList.from_string('(SUCCESS ' + gene_str + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :genes NIL)')
        return reply
        
    def respond_is_gene_tissue(self, content):
        """
        Respond to is-tissue-gene request
        Is stat3 expressed in liver?
        """
        #get tissue name
        tissue_name = _get_tissue_name(content)
        if not tissue_name:
            reply = make_failure('NO_TISSUE_NAME')
            return reply
        if tissue_name not in self.tissue_list:
            reply = make_failure('INVALID_TISSUE')
            return reply
            
        gene_name = get_gene(content, descr='gene')
        if not gene_name:
            gene_arg = content.gets('gene')
            reply = _wrap_family_message(gene_arg, 'NO_GENE_NAME')
            return reply
            
        #optional keyword
        keyword_name = _get_keyword_name(content, descr='keyword')
        if keyword_name == 'exclusive':
            try:
                is_express = self.tfta.is_tissue_gene_exclusive(tissue_name, gene_name)
            except Exception as e:
                reply = KQMLList.from_string('(SUCCESS :result FALSE)')
                return reply
        else:
            try:
                is_express = self.tfta.is_tissue_gene(tissue_name, gene_name)
            except Exception as e:
                reply = KQMLList.from_string('(SUCCESS :result FALSE)')
                return reply
        reply = KQMLList('SUCCESS')
        is_express_str = 'TRUE' if is_express else 'FALSE'
        reply.set('result', is_express_str)
        return reply
        
#
    def respond_get_family_name(self, content):
        family_arg = content.gets('family')
        
        agent = _get_family_name(family_arg)
        if agent:
            res_str = ' (:for ' + agent[0].name + ' :error FAMILY_NAME_NOT_ALLOWED)'
            reply = make_failure(res_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :family NIL)')
        return reply
        
    
        
    task_func = {'IS-REGULATION':respond_is_regulation, 'FIND-TF':respond_find_tf,
                 'FIND-PATHWAY':respond_find_pathway, 'FIND-TARGET':respond_find_target,
                 'FIND-GENE':respond_find_gene, 'FIND-MIRNA':respond_find_miRNA,
                 'FIND-COMMON-PATHWAY-GENES':respond_find_common_pathway_all,
                 'FIND-MIRNA-COUNT-GENE':respond_find_miRNA_count_target,
                 'FIND-GENE-COUNT-MIRNA':respond_find_target_count_miRNA,
                 'IS-MIRNA-TARGET':respond_is_miRNA_target2,
                 'FIND-COMMON-TF-GENES':respond_find_common_tfs_genes,
                 'IS-PATHWAY-GENE':respond_is_pathway_gene,
                 'FIND-TARGET-MIRNA':respond_find_target_miRNA,
                 'FIND-MIRNA-TARGET':respond_find_miRNA_target,
                 'IS-GENE-ONTO':respond_is_gene_onto,
                 'FIND-GENE-ONTO':respond_find_gene_onto,
                 'IS-TF-TARGET':respond_is_tf_target,
                 'FIND-TF-TARGET':respond_find_tf_targets,
                 'FIND-TARGET-TF':respond_find_target_tfs,
                 'FIND-PATHWAY-GENE':respond_find_pathway_gene,
                 'FIND-PATHWAY-GENE-KEYWORD':respond_find_pathway_gene_keyword,
                 'FIND-PATHWAY-DB-GENE':respond_find_pathway_db_gene,
                 'FIND-TF-PATHWAY':respond_find_tf_pathway,
                 'FIND-GENE-PATHWAY':respond_find_gene_pathway,
                 'FIND-PATHWAY-KEYWORD':respond_find_pathway_keyword,
                 'FIND-TF-KEYWORD':respond_find_tf_keyword,
                 'FIND-COMMON-PATHWAY-GENES-KEYWORD':respond_find_common_pathway_genes_keyword,
                 'FIND-GENE-GO-TF':respond_find_genes_go_tf2,
                 'IS-TF-TARGET-TISSUE':respond_is_tf_target_tissue,
                 'FIND-TF-TARGET-TISSUE':respond_find_tf_targets_tissue,
                 'FIND-TARGET-TF-TISSUE':respond_find_target_tfs_tissue,
                 'FIND-EVIDENCE-MIRNA-TARGET':respond_find_evidence_miRNA_target,
                 'FIND-PATHWAY-DB-KEYWORD':respond_find_pathway_db_keyword,
                 'FIND-TISSUE':respond_find_tissue_gene,
                 'FIND-KINASE-REGULATION':respond_find_kinase_regulation,
                 'FIND-TF-MIRNA':respond_find_tf_miRNA,
                 'FIND-REGULATION':respond_find_regulation,
                 'FIND-EVIDENCE':respond_find_evidence,
                 'FIND-GENE-TISSUE':respond_find_gene_tissue,
                 'IS-GENE-TISSUE':respond_is_gene_tissue,
                 'FIND-KINASE-PATHWAY':respond_find_kinase_pathway,
                 'TEST-FUNCTION':respond_get_family_name}
    
    def receive_request(self, msg, content):
        """If a "request" message is received, decode the task and
        the content and call the appropriate function to prepare the
        response. A reply message is then sent back.
        """
        task_str = content.head().upper()
        try:
            reply_content = self.task_func[task_str](self, content)
        except KeyError:
            self.error_reply(msg, 'unknown request task ' + task_str)
            return
        reply_msg = KQMLPerformative('reply')
        reply_msg.set('content', reply_content)
        self.reply(msg, reply_msg)
        
    def get_regulator_indra0(self, target_names, stmt_types, keyword_name):
        """
        wrap message for multiple targets case
        target_names: list
        stmt_types: indra statement type
        """
        lit_messages = ''
        tfs = defaultdict(set)
        others = defaultdict(set)
        mirnas = defaultdict(set)
        genes = defaultdict(set)
        for target in target_names:
            term_tuple = (target, 'target', keyword_name)
            if term_tuple not in self.stmts_indra:
                stmts,success = self.tfta.find_statement_indraDB(obj=target, stmt_types=stmt_types)
                if success:
                    self.stmts_indra[term_tuple] = stmts
            else:
                stmts = self.stmts_indra[term_tuple]
            #provenance support
            self.send_background_support(stmts, 'what', target, keyword_name)
            if len(stmts):
                tfs[target], genes[target], mirnas[target], others[target] = self.tfta.find_regulator_indra(stmts)
        #take the intersection
        ftfs = tfs[target_names[0]]
        fmirnas = mirnas[target_names[0]]
        fothers = others[target_names[0]]
        fgenes = genes[target_names[0]]
        if len(target_names)>1:
            for i in range(1, len(target_names)):
                ftfs = ftfs.intersection(tfs[target_names[i]])
                fmirnas = fmirnas.intersection(mirnas[target_names[i]])
                fothers = fothers.intersection(others[target_names[i]])
                fgenes = fgenes.intersection(genes[target_names[i]])
        if len(ftfs):
            lit_messages += self.wrap_message(':tf-literature', ftfs)
        if len(fgenes):
            lit_messages += self.wrap_message(':gene-literature', fgenes)
        if len(fmirnas):
            lit_messages += self.wrap_message(':miRNA-literature', fmirnas, hgnc_id=False)
        if len(fothers):
            lit_messages += self.wrap_message(':other-literature', fothers, hgnc_id=False)
        return lit_messages
        
    def get_regulator_indra(self, target_names, stmt_types, keyword_name):
        """
        wrap message for multiple targets case
        target_names: list
        stmt_types: indra statement type
        """
        lit_messages = ''
        tfs = defaultdict(set)
        others = defaultdict(set)
        mirnas = defaultdict(set)
        genes = defaultdict(set)
        stmts_d = defaultdict(list)
        if len(target_names) > 1:
            target_str = ', '.join(target_names[:-1]) + ' and ' + target_names[-1]
        else:
            target_str = target_names[0]
        for target in target_names:
            term_tuple = (target, 'target', keyword_name)
            if term_tuple not in self.stmts_indra:
                stmts,success = self.tfta.find_statement_indraDB(obj=target, stmt_types=stmt_types)
                if success:
                    self.stmts_indra[term_tuple] = stmts
            else:
                stmts = self.stmts_indra[term_tuple]
            if len(stmts):
                stmts_d[target] = stmts
            else:
                self.send_background_support([], 'what', target_str, keyword_name)
                return lit_messages
        tfs, genes, mirnas, others, stmt_f = self.tfta.find_regulators_indra(stmts_d)
        
        if len(tfs):
            lit_messages += self.wrap_message(':tf-literature', tfs)
        if len(genes):
            lit_messages += self.wrap_message(':gene-literature', genes)
        if len(mirnas):
            lit_messages += self.wrap_message(':miRNA-literature', mirnas, hgnc_id=False)
        if len(others):
            lit_messages += self.wrap_message(':other-literature', others, hgnc_id=False)
            
        #provenance support
        self.send_background_support(stmt_f, 'what', target_str, keyword_name)
        
        return lit_messages
        
    def get_regulator_indra_those(self, target_names, stmt_types, keyword_name, of_those=None, target_type=None):
        """
        wrap message for multiple targets case
        target_names: list
        stmt_types: indra statement type
        of_those: set
        """
        lit_messages = ''
        stmts_d = defaultdict(list)
        if len(target_names) > 1:
            target_str = ', '.join(target_names[:-1]) + ' and ' + target_names[-1]
        else:
            target_str = target_names[0]
        for target in target_names:
            term_tuple = (target, 'target', keyword_name)
            if term_tuple not in self.stmts_indra:
                stmts,success = self.tfta.find_statement_indraDB(obj=target, stmt_types=stmt_types)
                if success:
                    self.stmts_indra[term_tuple] = stmts
            else:
                stmts = self.stmts_indra[term_tuple]
            if len(stmts):
                stmts_d[target] = stmts
            else:
                self.send_background_support([], 'what', target_str, keyword_name)
                return lit_messages
        genes, stmt_f = self.tfta.find_regulator_indra_targets(stmts_d, of_those=of_those, target_type=target_type)
        
        if len(genes):
            lit_messages += self.wrap_message(':gene-literature', genes)
        #provenance support
        self.send_background_support(stmt_f, 'what', target_str, keyword_name)
        
        return lit_messages
        
    def get_target_indra(self, regulator_names, stmt_types, keyword_name, of_those=None, target_type=None):
        """
        wrap message for multiple regulators case
        regulator_names: list
        stmt_types: indra statement type
        keyword_name: str
        """
        lit_messages = ''
        stmts_d = defaultdict(list)
        if len(regulator_names) > 1:
            regulator_str = ', '.join(regulator_names[:-1]) + ' and ' + regulator_names[-1]
        else:
            regulator_str = regulator_names[0]
        for regulator in regulator_names:
            term_tuple = (regulator,'regulator',keyword_name)
            if term_tuple not in self.stmts_indra:
                stmts,success = self.tfta.find_statement_indraDB(subj=regulator, stmt_types=stmt_types)
                if success:
                    self.stmts_indra[term_tuple] = stmts
            else:
                stmts = self.stmts_indra[term_tuple]
            if len(stmts):
                stmts_d[regulator] = stmts
            else:
                self.send_background_support([], regulator_str, 'what', keyword_name)
                return lit_messages
        genes, stmt_f = self.tfta.find_target_indra_regulators(stmts_d, of_those=of_those, target_type=target_type)
        
        if len(genes):
            lit_messages += self.wrap_message(':targets', genes)
        #provenance support
        self.send_background_support(stmt_f, regulator_str, 'what', keyword_name)
        
        return lit_messages
        
    def get_tf_indra(self, target_names, stmt_types, keyword_name):
        """
        wrap message for multiple regulators case
        target_names: list
        stmt_types: indra statement type
        keyword_name: str
        """
        lit_messages = ''
        stmts_d = defaultdict(list)
        if len(target_names) > 1:
            target_str = ', '.join(target_names[:-1]) + ' and ' + target_names[-1]
        else:
            target_str = target_names[0]
        for target in target_names:
            term_tuple = (target, 'target', keyword_name)
            if term_tuple not in self.stmts_indra:
                stmts,success = self.tfta.find_statement_indraDB(obj=target, stmt_types=stmt_types)
                if success:
                    self.stmts_indra[term_tuple] = stmts
            else:
                stmts = self.stmts_indra[term_tuple]
            if len(stmts):
                stmts_d[target] = stmts
            else:
                self.send_background_support([], target_str, 'what', keyword_name)
                return lit_messages
        tfs, stmt_f = self.tfta.find_tfs_indra(stmts_d)        
        
        if len(tfs):
            lit_messages += self.wrap_message(':tfs', tfs)
        #provenance support
        self.send_background_support(stmt_f, target_str, 'what', keyword_name)
        
        return lit_messages


    def send_table_to_provenance_mirna(self, mirna_name, target_name, experiment, support_type, pmid, nl_question):
        """
        Post a concise table listing evidence found for mirna-target relationship.
        mirna_name: list
        target_name: list
        experiment: list
        support_type: list
        pmid: list
        nl_question: str
        """
        publink = "https://www.ncbi.nlm.nih.gov/pubmed/"
        head_str = '<head><style>.table-borders{border:1px solid black;padding:8px;}</style></head>'
        html_str = head_str + '<h4>Supporting information from TFTA: %s</h4>\n' % nl_question
        #html_str = '<h4>Supporting information from TFTA: %s</h4>\n' % nl_question
        html_str += '<table class = table-borders>\n'
        row_list = ['<th class = table-borders>MiRNA</th><th class = table-borders>Target</th> \
                   <th class = table-borders>Experiment</th><th class = table-borders>Support Type</th> \
                   <th class = table-borders>PMID</th>']
        for mirna,target,expe,st,pd in zip(mirna_name, target_name, experiment, support_type, pmid):
            pd_list = pd.split(',')
            pd_str = ''
            for p in pd_list:
                pd_str += '<a href=' + publink + p + ' target="_blank">' + p + '</a>;\n'
            row_list.append('<td class = table-borders>%s</td><td class = table-borders>%s</td> \
                             <td class = table-borders>%s</td><td class = table-borders>%s</td> \
                             <td class = table-borders>%s</td>' % (mirna, target, expe, st, pd_str[:-2]))
        html_str += '\n'.join(['  <tr>%s</tr>\n' % row_str
                               for row_str in row_list])
        html_str += '</table> <hr>'
        content = KQMLList('add-provenance')
        content.sets('html', html_str)
        return self.tell(content)
        
    def send_background_support_mirna(self, mirna_name, target_name, experiment, support_type, pmid, find_mirna=False, find_target=False):
        """
        Send the evidence from the MiRNA-target database
        mirna_name: list or str
        target_name: list or str
        experiment: dict
        support_type: dict
        pmid: dict
        """
        if pmid:
            if all([type(mirna_name).__name__ == 'str', type(target_name).__name__ == 'str']):
                nl = 'does ' + mirna_name + ' regulate ' + target_name + '?'
                exp_str = ';\n'.join(experiment[target_name])
                suptype_str = ';\n'.join(support_type[target_name])
                pmid_str = ','.join(pmid[target_name])
                self.send_table_to_provenance_mirna([mirna_name], [target_name], [exp_str], [suptype_str], [pmid_str], nl)
            elif find_mirna:
                nl = 'what miRNAs regulate ' + target_name + '?'
                target_list = []
                exp_list = []
                sup_list = []
                pmid_list = []
                for i in range(len(mirna_name)):
                    target_list.append(target_name)
                    exp_list.append(';\n'.join(experiment[mirna_name[i]]))
                    sup_list.append(';\n'.join(support_type[mirna_name[i]]))
                    pmid_list.append(','.join(pmid[mirna_name[i]]))
                self.send_table_to_provenance_mirna(mirna_name, target_list, exp_list, sup_list, pmid_list, nl)
            elif find_target:
                nl = 'what genes are regulated by ' + mirna_name + '?'
                mirna_list = []
                exp_list = []
                sup_list = []
                pmid_list = []
                for i in range(len(target_name)):
                    mirna_list.append(mirna_name)
                    exp_list.append(';\n'.join(experiment[target_name[i]]))
                    sup_list.append(';\n'.join(support_type[target_name[i]]))
                    pmid_list.append(','.join(pmid[target_name[i]]))
                self.send_table_to_provenance_mirna(mirna_list, target_name, exp_list, sup_list, pmid_list, nl)
            else:
                return
        else:
            for_what = 'your query'
            cause_txt = 'MiRNA-target db'
            reason_txt = ''
            self.send_null_provenance(stmt=for_what, for_what=cause_txt, reason=reason_txt)
        
    def send_table_to_provenance(self, tf_name, target_name, dbname, nl_question):
        """
        Post a concise table listing evidence found for tf-target relationship.
        tf_name: list
        target_name: list
        dbname: list
        nl_question: str
        """
        publink = "https://www.ncbi.nlm.nih.gov/pubmed/"
        #head_str = '<head><style>table,th,td{border:1px solid black;padding:8px}</style></head>'
        head_str = '<head><style>.table-borders{border:1px solid black;padding:8px;}</style></head>'
        html_str = head_str + '<h4>Supporting information from TFTA: %s</h4>\n' % nl_question
        #html_str = '<h4>Supporting information from TFTA: %s</h4>\n' % nl_question
        html_str += '<table class = table-borders>\n'
        row_list = ['<th class = table-borders>TF</th><th class = table-borders>Target</th> \
                    <th class = table-borders>Source</th>']
        for tf,target,db in zip(tf_name, target_name, dbname):
            db_list = db.split(',')
            db_str = ''
            for d in db_list:
                db_str += '<a href=' + publink + dbname_pmid_map[d] + ' target="_blank">' + d + '</a>,'
            row_list.append('<td class = table-borders>%s</td><td class = table-borders>%s</td> \
                             <td class = table-borders>%s</td>' % (tf, target, db_str[:-1]))
        html_str += '\n'.join(['  <tr>%s</tr>\n' % row_str
                               for row_str in row_list])
        html_str += '</table> <hr>'
        content = KQMLList('add-provenance')
        content.sets('html', html_str)
        return self.tell(content)
        
    def send_background_support_db(self, tf_name, target_name, dbname, tissue=None, find_tf=False, find_target=False):
        """
        Send the evidence from the tf-target db
        """
        if dbname:
            if all([type(tf_name).__name__ == 'str', type(target_name).__name__ == 'str']):
                if tissue:
                    nl = 'does ' + tf_name + ' regulate ' + target_name + ' in ' + tissue + '?'
                else:
                    nl = 'does ' + tf_name + ' regulate ' + target_name + '?'
                self.send_table_to_provenance([tf_name], [target_name], [dbname], nl)
            elif find_tf:
                if len(target_name) > 1:
                    target_str = ', '.join(target_name[:-1]) + ' and ' + target_name[-1]
                else:
                    target_str = target_name[0]
                if tissue:
                    nl = 'what transcription factors regulate ' + target_str + ' in ' + tissue + '?'
                else:
                    nl = 'what transcription factors regulate ' + target_str + '?'
                target_list = []
                db_list = []
                tf_list = []
                for tf in tf_name:
                    for target in target_name:
                        tf_list.append(tf)
                        target_list.append(target)
                        db_list.append(dbname[(tf,target)])
                self.send_table_to_provenance(tf_list, target_list, db_list, nl)
            elif find_target:
                if len(tf_name) > 1:
                    tf_str = ', '.join(tf_name[:-1]) + ' and ' + tf_name[-1]
                else:
                    tf_str = tf_name[0]
                if tissue:
                    nl = 'what genes are regulated by ' + tf_str + ' in ' + tissue + '?'
                else:
                    nl = 'what genes are regulated by ' + tf_str + '?'
                tf_list = []
                target_list = []
                db_list = []
                for target in target_name:
                    for tf in tf_name:
                        tf_list.append(tf)
                        target_list.append(target)
                        db_list.append(dbname[(tf,target)])
                self.send_table_to_provenance(tf_list, target_list, db_list, nl)
            else:
                return
        else:
            for_what = 'your query'
            cause_txt = 'tf-db'
            reason_txt = ''
            self.send_null_provenance(stmt=for_what, for_what=cause_txt, reason=reason_txt)
        
    def get_kinase_regulation(self, target_names, keyword_name, of_those=None):
        """
        target_names: list
        keyword_name: string
        """
        kin_messages = ''
        kinase_names = []
        if keyword_name == 'regulate':
            try:
                kinase_names = self.tfta.find_kinase_target(target_names)
            except KinaseNotFoundException:
                return kin_messages
        else:
            try:
                kinase_names = self.tfta.find_kinase_target_keyword(target_names, keyword_name)
            except KinaseNotFoundException:
                return kin_messages
        
        if of_those:
            kinase_names = list(set(kinase_names).intersection(set(of_those)))
        if len(kinase_names):
            kin_messages += self.wrap_message(':kinase-db', kinase_names)
        return kin_messages
        
    def send_background_support(self, stmts, regulator_name, target_name, keyword_name):
        logger.info('Sending support for %d statements' % len(stmts))
        if target_name == 'what':
            interaction = stmt_provenance_map_target[keyword_name.lower()]
            for_what = 'what genes are ' + interaction + ' by ' + regulator_name
        else:
            interaction = stmt_provenance_map[keyword_name.lower()]
            for_what = regulator_name + ' ' + interaction + ' ' + target_name
        if len(stmts):
            self.send_provenance_for_stmts(stmts, for_what, limit = 50)
        else:
            cause_txt = 'literature'
            reason_txt = ''
            self.send_null_provenance(stmt=for_what, for_what=cause_txt, reason=reason_txt)
    
    def send_null_provenance(self, stmt, for_what, reason=''):
        """Send out that no provenance could be found for a given Statement."""
        content_fmt = ('<h4>No supporting evidence found for {statement} from '
                        '{cause}{reason}.</h4>')
        content = KQMLList('add-provenance')
        content.sets('html', content_fmt.format(statement=stmt,
                                                cause=for_what, reason=reason))
        return self.tell(content)
        
    def wrap_message(self, descr, data, hgnc_id=True):
        if type(data) is set:
            data = list(data)
        if not self.hgnc_info:
            self.hgnc_info = self.tfta.get_hgnc_mapping()
        if not data:
            return None
        data.sort()  
        tf_list_str = ''
        #don't consider strings from literature containing double quote
        dquote = '"'
        if hgnc_id:
            for tf in data:
                if dquote not in tf:
                    tf_str = '"' + tf + '"'
                    try:
                        id = self.hgnc_info[tf]
                    except KeyError:
                        id = None
                    if id:
                        tf_list_str += '(:name %s :hgnc %s) ' % (tf_str, str(id))
                    else:
                        tf_list_str += '(:name %s) ' % tf_str
        else:
            for tf in data:
                if dquote not in tf:
                    tf_str = '"' + tf + '"'
                    tf_list_str += '(:name %s) ' % tf_str
        return descr + ' (' + tf_list_str + ') '
    
    def get_target_type_set(self, content):
        target_type = _get_keyword_name(content, descr='target-type')
        target_type_set = set()
        if target_type:
            try:
                target_type_set = self.tfta.get_onto_set(target_type)
            except GONotFoundException:
                return target_type_set
        return target_type_set
        
    def wrap_family_message1(self, target_arg, msg):
        term_id = _get_term_id(target_arg)
        if len(term_id):
            members = self.tfta.find_members(term_id)
        if members:
            res_str = ''
            for id in members:
                for a in members[id]:
                    res_str += '(:name %s)' % a.name
                res_str = '(:term %s :as (%s))' % (id, res_str)
            res_str = '(resolve :agents (' + res_str + '))'
            reply = make_failure_clarification(msg, res_str)
        else:
            reply = make_failure(msg)
        return reply
    
    def get_reply_mirna_clarification(self, miRNA_mis):
        try:
            clari_mirna = self.tfta.get_similar_miRNAs(list(miRNA_mis.keys())[0])
            c_str = _wrap_mirna_clarification(miRNA_mis, clari_mirna)
            reply = make_failure_clarification('MIRNA_NOT_FOUND', c_str)
            return reply
        except miRNANotFoundException:
            reply = make_failure('NO_SIMILAR_MIRNA')
            return reply
    
    def is_protein_gene(self, gene_arg, keyword):
        gene_map = {'gene':['ONT::GENE-PROTEIN', 'ONT:GENE'], 'protein': ['ONT::GENE-PROTEIN', 'ONT:PROTEIN']}
        #gene_arg = content.gets('gene')
        is_onto = False
        is_ekb = True
        if '<ekb' in gene_arg or '<EKB' in gene_arg:
            tp = TripsProcessor(gene_arg)
            for term in tp.tree.findall('TERM'):
                if term.find('type').text in gene_map[keyword]:
                    is_onto = True
                    break
        else:
            is_ekb = False
        return is_onto, is_ekb
         
#------------------------------------------------------------------------#######
def _get_target(target_str):
    agent = None
    ont1 = ['ONT::PROTEIN', 'ONT::GENE-PROTEIN', 'ONT::GENE']
    tp = TripsProcessor(target_str)
    for term in tp.tree.findall('TERM'):
        if term.find('type').text in ont1:
            term_id = term.attrib['id']
            agent = tp._get_agent_by_id(term_id, None)
            break
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
    
def _get_family_name(target_arg):
    agent = []
    ont1 = ['ONT::PROTEIN-FAMILY', 'ONT::GENE-FAMILY']
    tp = TripsProcessor(target_arg)
    for term in tp.tree.findall('TERM'):
        if term.find('type').text in ont1:
            term_id = term.attrib['id']
            agent.append(tp._get_agent_by_id(term_id, None))
    return agent
    
def _get_term_id(target_arg, otype=1):
    term_id = dict()
    if '<ekb' in target_arg or '<EKB' in target_arg:
        if otype == 1:
            ont1 = ['ONT::PROTEIN-FAMILY', 'ONT::GENE-FAMILY']
        else:
            ont1 = ['ONT::RNA', 'ONT::GENE']
        tp = TripsProcessor(target_arg)
        for term in tp.tree.findall('TERM'):
            if term.find('type').text in ont1:
                id = term.attrib['id']
                term_id[id] = tp._get_agent_by_id(id, None)
    #print('term_id = ' + ','.join(term_id))
    return term_id

def _wrap_family_message(target_arg, msg):
    term_id = _get_term_id(target_arg)
    if len(term_id):
        res_str = ''
        for t in term_id:
            res_str += '(:for ' + t + ' :error FAMILY_NAME_NOT_ALLOWED) '
        res_str = '(' + res_str + ')'
        reply = make_failure(res_str)
    else:
        reply = make_failure(msg)
    return reply
    
def _wrap_family_message2(target_arg, msg):
    family = _get_family_name(target_arg)
    family_name = []
    for f in family:
        family_name.append(f.name)
    if family_name:
        res_str = ''
        for f in family_name:
            res_str += '(:for ' + f + ' :error FAMILY_NAME_NOT_ALLOWED) '
        res_str = '(' + res_str + ')'
        reply = make_failure(res_str)
    else:
        reply = make_failure(msg)
    return reply

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
    miRNA_names = dict()
    try:
        root = ET.fromstring(xml_string)
    except Exception as e:
        return miRNA_names
    ont1 = ['ONT::RNA', 'ONT::GENE']
    try:
        for term in root.findall('TERM'):
            if term.find('type').text in ont1:
                s = term.find('name')
                if s is not None:
                    s1 = s.text
                    s1 = rtrim_hyphen(s1)
                    miRNA_names[s1.upper()] = term.attrib['id']
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
    

def _wrap_evidence_message(descr, evids, limit = 10):
    """
    descr: descriptor, for example: ':evidence'
    evids: set of tuple(source_api, pmid, text)
    """
    evi_message = ''
    ind = 0
    for ev in evids:
        if all(ev) and ind < limit:
            evi_message += '('
            evi_message += ':source_api ' + ev[0]
            evi_message += ' :pmid ' + ev[1]
            evi_message += ' :text ' + '"' + ev[2] + '"'
            evi_message += ') '
            ind += 1
        else:
            break
    evi_message = descr + ' (' + evi_message + ') '
    return evi_message
    
def _combine_messages(mess_list):
    messages = ''
    for mess in mess_list:
        if len(mess):
            messages += mess
    return messages
    
def _wrap_mirna_clarification(miRNA_mis, clari_mirna):
    """
    miRNA_mis: dict
    clari_mirna: list
    """
    c_str = ''
    for c in clari_mirna:
        c_str += '(:name %s) ' % c
    c_str = '(resolve :term ' + list(miRNA_mis.values())[0] + ' :as (' + c_str + '))'
    return c_str

def _wrap_mirna_clarification2(miRNA_mis, clari_mirna):
    """
    miRNA_mis: dict
    clari_mirna: dict
    """
    mir_str = ''
    for mir in clari_mirna.keys():
        c_str = ''
        for c in clari_mirna[mir]:
            c_str += '(:name %s) ' % c
        c_str = '(:term ' + miRNA_mis[mir] + ' :as (' + c_str + '))'
        mir_str += c_str
    mir_str = '(resolve (' + mir_str + '))'
    return mir_str

def _wrap_pathway_message(pathwayName, dblink, keyword=None):
    """
    pathwayName: list
    dblink: list
    keyword: list or None
    """
    pathway_list_str = ''
    #limit the number of pathways to return
    limit = 50
    num = 0
    if keyword:
        if type(keyword).__name__ == 'str':
            keyword = [keyword]
        for pn, dbl in zip(pathwayName, dblink):
            num += 1
            if num > limit:
                break
            if _filter_subword(pn, keyword):
                pnslash = '"' + pn +'"'
                dbl = '"' + dbl +'"'
                pathway_list_str += '(:name %s :dblink %s) ' % (pnslash, dbl)
        if pathway_list_str:
            pathway_list_str = '(' + pathway_list_str + ')'
            reply = KQMLList('SUCCESS')
            reply.set('pathways', pathway_list_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
    else:
        for pn, dbl in zip(pathwayName, dblink):
            num += 1
            if num > limit:
                break
            pnslash = '"' + pn + '"'
            dbl = '"' + dbl + '"'
            pathway_list_str += '(:name %s :dblink %s) ' % (pnslash, dbl)
        pathway_list_str = '(' + pathway_list_str + ')'
        reply = KQMLList('SUCCESS')
        reply.set('pathways', pathway_list_str)
    return reply
    
def _wrap_pathway_genelist_message(pathwayName, dblink, genelist, pathway_names=None,
                                   gene_descr=':tfs', of_gene_names=None):
    """
    parameters
    -------------
    pathwayName: dict[pthid]
    dblink: dict[pthid]
    genelist: dict[pthid]
    pathway_names: list or None
    gene_descr: string
    of_gene_names: list
    """
    limit = 50
    num = 0
    pathway_list_str = ''
    keys = list(genelist.keys())
    if keys:
        if pathway_names:
            for key in keys:
                if _filter_subword(pathwayName[key], pathway_names):
                    gene_list_str = ''
                    if of_gene_names:
                        genes = set(of_gene_names) & set(genelist[key])
                    else:
                        genes = genelist[key]
                    if genes:
                    #check the limit
                        num += 1
                        if num > limit:
                            break
                        for gene in genes:
                            gene_list_str += '(:name %s) ' % gene
                        gene_list_str = '(' + gene_list_str + ')'
                        pn = '"' + pathwayName[key] + '"'
                        dl = '"' + dblink[key] + '"'
                        pathway_list_str += '(:name %s :dblink %s %s %s) ' % (pn, dl, gene_descr, gene_list_str)
        else:
            for key in keys:
                gene_list_str = ''
                if of_gene_names:
                    genes = set(of_gene_names) & set(genelist[key])
                else:
                    genes = genelist[key]
                if genes:
                #check the limit
                    num += 1
                    if num > limit:
                        break
                    for gene in genes:
                        gene_list_str += '(:name %s) ' % gene
                    gene_list_str = '(' + gene_list_str + ')'
                    pn = '"' + pathwayName[key] + '"'
                    dl = '"' + dblink[key] + '"'
                    pathway_list_str += '(:name %s :dblink %s %s %s) ' % (pn, dl, gene_descr, gene_list_str)
        if pathway_list_str:
            reply = KQMLList.from_string('(SUCCESS :pathways (' + pathway_list_str + '))')
        else:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
    else:
        reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
    return reply
    
def _filter_subword(sentence, pattern_list):
    """
    If the sentence contains any pattern in the pattern_list as a single word or phrase, 
    return true, else false
    """
    word = False
    sen_list = sentence.lower().split(' ')
    for p in pattern_list:
        ps = p.split(' ')
        ind = True
        for s in ps:
            ind = ind & (s.lower() in sen_list)
        if ind:
            word = True
            break
    return word
    
def get_of_those_list(content, descr='of-those'):
    """
    return a list of genes by parsing of_those_arg
    of_those_arg: str or EKB xml
    """
    #check if it's using ekb xml format
    gene_names = []
    gene_arg = content.gets(descr)
    if gene_arg:
        if '<ekb' in gene_arg or '<EKB' in gene_arg:
            genes = _get_targets(gene_arg)
            for gene in genes:
                gene_names.append(gene.name)
        else:
            gene_arg = content.get(descr)
            gene_arg_str = gene_arg.data
            gene_arg_str = gene_arg_str.replace(' ', '')
            gene_arg_str = gene_arg_str.upper()
            gene_names = gene_arg_str.split(',')
    return gene_names
    
def get_of_those_mirna(content, descr='of-those'):
    mirna_names = []
    mirna_arg = content.gets(descr)
    if mirna_arg:
        if '<ekb' in mirna_arg or '<EKB' in mirna_arg:
            mirna_names = set(_get_miRNA_name(mirna_arg).keys())
        else:
            mirna_arg = content.get(descr)
            mirna_str = mirna_arg.data
            mirna_str = mirna_str.replace(' ', '')
            mirna_names = mirna_str.split(',')
    return mirna_names
    
def get_gene(content, descr='gene'):
    gene_name = ''
    gene_arg = content.gets(descr)
    if gene_arg:
        if '<ekb' in gene_arg or '<EKB' in gene_arg:
            gene = _get_target(gene_arg)
            if gene:
                gene_name = gene.name
        else:
            gene_arg = content.get(descr)
            gene_arg_str = gene_arg.data
            gene_arg_str = gene_arg_str.replace(' ', '')
            gene_arg_str = gene_arg_str.upper()
            gene_name = gene_arg_str.split(',')[0]
    return gene_name
    
def _get_tissue_name(content):
    tissue_name = ''
    try:
        tissue_arg = content.get('tissue')
        #tissue_name = tissue_arg.head()
        tissue_name = tissue_arg.data
        #tissue_name = trim_quotes(tissue_name)
        tissue_name = tissue_name.lower()
        tissue_name = tissue_name.replace(' ', '_')
        tissue_name = tissue_name.replace('-', '_')
    except Exception as e:
        return tissue_name
    return tissue_name

def _get_keyword_name(content, descr='keyword', hyphen=False, low_case=True):
    keyword_name = ''
    try:
        keyword_arg = content.get(descr)
        keyword_name = keyword_arg.data
        if low_case:
            keyword_name = keyword_name.lower()
        if hyphen:
            keyword_name = keyword_name.replace('-', ' ')
    except Exception as e:
        return keyword_name
    return keyword_name
            

if __name__ == "__main__":
    TFTA_Module(argv=sys.argv[1:])
