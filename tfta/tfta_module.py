"""TFTA module is to receive and decode messages and send responses from
and to other agents in the system"""

import sys
import re
import time
#import xml.etree.ElementTree as ET
import logging
logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA')
from kqml import KQMLModule, KQMLPerformative, KQMLList
from .tfta import TFTA, TFNotFoundException, TargetNotFoundException, PathwayNotFoundException 
from .tfta import GONotFoundException, miRNANotFoundException, TissueNotFoundException
from .tfta import KinaseNotFoundException
from .mirDisease import mirDisease
from enrichment.GO import GOEnrich
from enrichment.pathway import PathwayEnrich
from utils.heatmap import generate_heatmap
#from indra.sources.trips.processor import TripsProcessor
from collections import defaultdict
from bioagents import Bioagent
from indra.statements import Agent
from indra.databases import hgnc_client

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
             'FIND-TF-PATHWAY', 'FIND-GENE-PATHWAY',
             'FIND-PATHWAY-KEYWORD', 'FIND-TF-KEYWORD',
             'FIND-COMMON-TF-GENES', 'FIND-GENE-GO-TF',
             'FIND-COMMON-PATHWAY-GENES', 'FIND-GENE-TISSUE',
             'IS-TF-TARGET-TISSUE', 'FIND-EVIDENCE-MIRNA-TARGET',
             'FIND-TARGET-TF-TISSUE', 'IS-PATHWAY-GENE','IS-MIRNA-TARGET', 
             'FIND-MIRNA-TARGET', 'FIND-TARGET-MIRNA', 'GO-ANNOTATION',
             'FIND-MIRNA-COUNT-GENE','FIND-GENE-COUNT-MIRNA',
             'FIND-PATHWAY-DB-KEYWORD', 'FIND-TISSUE', 'IS-REGULATION',
             'FIND-TF', 'FIND-PATHWAY', 'FIND-TARGET', 'FIND-GENE', 'FIND-MIRNA',
             'IS-GENE-ONTO', 'FIND-GENE-ONTO', 'FIND-KINASE-REGULATION',
             'FIND-TF-MIRNA', 'FIND-REGULATION', 'FIND-EVIDENCE', 
             'IS-GENE-TISSUE', 'FIND-KINASE-PATHWAY', 'GO-ENRICHMENT', 
             'IS-MIRNA-DISEASE', 'FIND-MIRNA-DISEASE', 'FIND-DISEASE-MIRNA',
             'MAKE-HEATMAP', 'PATHWAY-ENRICHMENT', 'DISEASE-ENRICHMENT',
             'MIRNA-DISEASE-ENRICHMENT', 'FIND-EVIDENCE-MIRNA-EXP']
    #keep the genes from the most recent previous call, which are used to input 
    #find-gene-onto if there's no gene input 
    #gene_list = ['STAT3', 'JAK1', 'JAK2', 'ELK1', 'FOS', 'SMAD2', 'KDM4B']
    gene_list = []

    def __init__(self, **kwargs):
        #Instantiate a singleton TFTA agent and other agents
        self.tfta = TFTA()
        self.go = GOEnrich()
        self.md = mirDisease()
        self.mirna_pop = self.md.get_mirna_pop()
        self.pop = self.tfta.get_hgnc_symbols()
        self.pw = PathwayEnrich(self.pop, self.mirna_pop, self.tfta)
        
        self.stmts_indra = dict()
        #self.hgnc_info = dict()
    
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
        tf_arg = content.get('tf')
        target_arg = content.get('target')
        tissue_arg = content.get('tissue')
        keyword_name = _get_keyword_name(content, descr='keyword')
        
        if all([tf_arg,target_arg,tissue_arg]):
            reply = self.respond_is_tf_target_tissue(content)
        elif all([tf_arg,target_arg, keyword_name]):
            if keyword_name in ['increase', 'decrease']:
                reply = self.respond_is_tf_target_literature(content)
            else:
                reply = self.respond_is_tf_target(content)
        else:
            reply = self.respond_is_tf_target(content)
        
        return reply
    
    def respond_find_tf(self, content):
        """
        Response content to find-tf request which includes:
        find-target-tf, find-target-tf-tissue
        """
        target_arg = content.get('target')
        tissue_arg = content.get('tissue')
        keyword_name = _get_keyword_name(content, descr='keyword')
        
        if all([target_arg,tissue_arg]):
            reply = self.respond_find_target_tfs_tissue(content)
        elif all([target_arg, keyword_name]):
            if keyword_name in ['increase', 'decrease']:
                reply = self.respond_find_tf_literature(content)
            else:
                reply = self.respond_find_target_tfs(content)
        else:
            reply = self.respond_find_target_tfs(content)
        
        return reply
    
    def respond_find_pathway(self, content):
        """
        Response content to find-pathway request, which includes:
        find-pathway-gene, find-pathway-db-gene, find-pathway-gene-keyword,
        find-pathway-keyword.
        """
        gene_arg = content.get('gene')
        regulator_arg = content.get('regulator')
        db_arg = content.get('database')
        keyword_arg = content.get('keyword')
        if all([keyword_arg,gene_arg]):
            reply = self.find_pathway_gene_keyword(content)
        elif all([gene_arg,db_arg]):
            reply = self.find_pathway_db_gene(content)
        elif all([keyword_arg,regulator_arg]):
            reply = self.find_pathway_keyword_regulator(content)
        elif all([regulator_arg,db_arg]):
            reply = self.find_pathway_db_regulator(content)
        elif regulator_arg:
            reply = self.find_pathway_regulator(content)
        elif keyword_arg:
            reply = self.respond_find_pathway_keyword(content)
        else:
            reply = self.find_pathway_gene(content)
        return reply
        
    def respond_find_common_pathway_genes(self, content):
        """
        Respond to find-common-pathway-genes, which covers:
        find-common-pathway-genes, find-common-pathway-genes-keyword,
        find-common-pathway-genes-database
        """
        target_arg = content.get('target')
        keyword_arg = content.get('keyword')
        db_arg = content.get('database')
        if all([target_arg, keyword_arg]):
            reply = self.get_common_pathway_genes_keyword(content)
        elif all([target_arg, db_arg]):
            reply = self.get_common_pathway_genes_db(content)
        else:
            reply = self.get_common_pathway_genes(content)
        
        return reply
    
    def respond_find_target(self, content):
        """
        Response content to find-target request, which includes these cases:
        find-tf-target, find-tf-target-tissue. 
        """
        tissue_arg = content.get('tissue')
        source_arg = content.get('source')
        
        if tissue_arg:
            reply = self.find_tf_targets_tissue(content)
        elif source_arg:
            reply = self.find_target_source(content)
        else:
            reply = self.find_target_all(content)
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
            reply = self.find_regulation_source(content)
        elif agent_arg:
            reply = self.find_regulation_agent(content)
        else:
            reply = self.find_regulation_all(content)
        return reply
        
    def respond_is_gene_onto(self, content):
        """
        Respond to is-gene-onto
        """
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        gene_name,term_id = self._get_targets(content, descr='gene')
        if not gene_name:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
            
        #check if keyword is protein or gene
        if keyword_name in ['gene', 'protein']:
            if gene_name:
                is_onto = True
            else:
                is_onto = False
        else:
            try:
                is_onto = self.tfta.Is_gene_onto(keyword_name, gene_name[0])
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
        regulator_arg = content.get('regulator')
        #gene_arg = content.get('gene')
        if regulator_arg:
            reply = self.respond_find_gene_onto_regulator(content)
        else:
            reply = self.respond_find_gene_onto_gene(content)
        
        return reply
        
    def respond_find_gene_onto_gene(self, content):
        """
        Respond to find-gene-onto
        """
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        gene_names,term_id = self._get_targets(content, descr='gene')
        if not gene_names:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
            
        try:
            results = self.tfta.find_gene_onto(keyword_name, gene_names)
        except GONotFoundException:
            reply = make_failure('GO_NOT_FOUND')
            return reply
        if len(results):
            gene_json = self._get_genes_json(results)
            reply = KQMLList('SUCCESS')
            reply.set('genes', gene_json)
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
            
        regulator_names,term_id = self._get_targets(content, descr='regulator')
        if not regulator_names:
            reply = self.wrap_family_message(term_id, 'NO_REGULATOR_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
                gene_json = self._get_genes_json(results)
                reply = KQMLList('SUCCESS')
                reply.set('genes', gene_json)
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
        tf_name,term_id = self._get_targets(content, descr='tf')
        if not tf_name:
            #tf_arg = content.gets('tf')
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
            return reply
        
        target_name,term_id = self._get_targets(content, descr='target')
        if not target_name:
            #target_arg = content.gets('target')
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        
        #check if it exists in literature
        term_tuple = (tf_name[0], target_name[0], 'REGULATE')
        if term_tuple not in self.stmts_indra:
            stmts,success = self.tfta.find_statement_indraDB(subj=tf_name[0], obj=target_name[0], stmt_types=['RegulateAmount'])
            if success:
                self.stmts_indra[term_tuple] = stmts
        else:
            stmts = self.stmts_indra[term_tuple]
        #provenance support
        self.send_background_support(stmts, tf_name[0], target_name[0], 'regulate')
        if len(stmts):
            literature = True
        #check the db
        try:
            is_target,dbname = self.tfta.Is_tf_target(tf_name[0], target_name[0])
            #provenance support
            self.send_background_support_db(tf_name[0], target_name[0], dbname)
        except Exception:
            #provenance support
            self.send_background_support_db(tf_name[0], target_name[0], '')
            if not literature:
                reply = KQMLList('SUCCESS')
                reply.set('result', 'FALSE')
                reply.set('db', 'FALSE')
                reply.set('literature', 'FALSE')
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
        tf_name,term_id = self._get_targets(content, descr='tf')
        if not tf_name:
            #tf_arg = content.gets('tf')
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
            return reply
            
        target_name,term_id = self._get_targets(content, descr='target')
        if not target_name:
            #target_arg = content.gets('target')
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
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
            
        term_tuple = (tf_name[0], target_name[0], keyword_name)
        if term_tuple not in self.stmts_indra:
            stmts,success = self.tfta.find_statement_indraDB(subj=tf_name[0], obj=target_name[0], stmt_types=stmt_types)
            if success:
                self.stmts_indra[term_tuple] = stmts
        else:
            stmts = self.stmts_indra[term_tuple]
        #provenance support
        self.send_background_support(stmts, tf_name[0], target_name[0], keyword_name)
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
        tf_name,term_id = self._get_targets(content, descr='tf')
        if not tf_name:
            #tf_arg = content.gets('tf')
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
            return reply
            
        target_name,term_id = self._get_targets(content, descr='target')
        if not target_name:
            #target_arg = content.gets('target')
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
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
            is_target = self.tfta.Is_tf_target_tissue(tf_name[0], target_name[0], tissue_name)
        except Exception as e:
            reply = KQMLList('SUCCESS')
            reply.set('result', 'FALSE')
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
        tf_names,term_id = self._get_targets(content, descr='tf')
        if not tf_names:
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply

        of_targets_names,nouse = self._get_targets(content, descr='of-those')
        
        #consider an optional parameter to only return given type of genes
        target_type_set,target_type = self.get_target_type_set(content)
        
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
        
        if target_names:
            reply = self.wrap_message('targets', target_names)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        
        #provenance support
        self.send_background_support_db(tf_names, target_names, dbname, find_target=True, target_type=target_type)
        
        return reply
        
    def respond_find_target_literature(self, content):
        """
        Respond to find-target request with information in literatures
        """
        tf_names,term_id = self._get_targets(content, descr='tf')
        if not tf_names:
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        #consider an optional parameter for sequencing query
        of_those_names,nouse = self._get_targets(content, descr='of-those')
        
        #consider an optional parameter to only return given types of genes
        target_type = _get_keyword_name(content, descr='target-type')
    
        try:
            stmt_types = stmt_type_map[keyword_name]
        except KeyError:
            reply = make_failure('INVALID_KEYWORD')
            return reply
        lit_messages = self.get_target_indra(tf_names, stmt_types, keyword_name, 
                                             of_those=of_those_names, target_type=target_type)
        if len(lit_messages):
            reply = KQMLList('SUCCESS')
            reply.set('targets', lit_messages)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply

    def find_tf_targets_tissue(self, content):
        """
        Response content to find-tf-target-tissue request
        For tf list, reply with the targets found within given tissue
        """
        tf_names,term_id = self._get_targets(content, descr='regulator')
        if not tf_names:
            reply = self.wrap_family_message(term_id, 'NO_REGULATOR_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
        of_targets_names,_ = self._get_targets(content, descr='of-those')
        
        #consider an optional parameter to only return given types of genes
        target_type_set,target_type = self.get_target_type_set(content)
        
        target_names = self.tfta.find_targets_tissue(tf_names, tissue_name)
            
        #check if it's a sequencing query
        if of_targets_names and target_names:
            target_names = list(set(of_targets_names) & set(target_names))
        #check if it requires returning a type of genes
        if target_type_set and target_names:
            target_names = target_type_set & set(target_names)
            
        if len(target_names):
            reply = KQMLList('SUCCESS')
            target_json = self._get_genes_json(target_names)
            tmess = KQMLList()
            tmess.set('target-db', target_json)
            reply.set('targets', tmess)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply

    def respond_find_target_tfs(self, content):
        """
        Response content to find-target-tf request
        For a target list, reply the tfs found
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
            
        #consider an optional parameter for subsequent query
        of_tfs_names,nouse = self._get_targets(content, descr='of-those')
        
        try:
            tf_names,dbname = self.tfta.find_tfs(target_names)
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
            #provenance support
            self.send_background_support_db([], target_names, '')
            return reply
        #check if it's a subsequent query
        if of_tfs_names:
            tf_names = sorted(list(set(of_tfs_names) & set(tf_names)))
        if len(tf_names):
            reply = self.wrap_message('tfs', tf_names)
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
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
            reply = KQMLList('SUCCESS')
            reply.set('tfs', lit_messages)
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        return reply

    def respond_find_target_tfs_tissue(self, content):
        """
        Response content to find-target-tf-tissue request
        For a target list, reply the tfs found within a given tissue
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
        of_tfs_names,nouse = self._get_targets(content, descr='of-those')
        
        tf_names = self.tfta.find_tfs_tissue(target_names, tissue_name)
            
        #check if it's a sequencing query
        if of_tfs_names and tf_names:
            tf_names = sorted(list(set(of_tfs_names) & set(tf_names)))
        if tf_names:
            reply = self.wrap_message('tfs', tf_names)
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        return reply

    def find_pathway_gene(self,content):
        """
        Response content to find-pathway request
        For a given gene list, return the related pathways information
        """
        gene_names,fmembers = self._get_targets2(content, descr='gene')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
            
        try:
            pathwayName,dblink = self.tfta.find_pathway_genes(gene_names, fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            #reply = self.wrap_family_message_pathway(term_id, descr='pathways', msg="PATHWAY_NOT_FOUND")
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink, of_those=of_those)
        return reply
        
    def find_pathway_regulator(self,content):
        """
        Response content to find-pathway request
        For a given regulator list, reply the pathways containing their targets
        """
        regulator_names,fmembers = self._get_targets2(content, descr='regulator')
        if not regulator_names and not fmembers:
            reply = make_failure('NO_REGULATOR_NAME')
            return reply
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
        
        #get genes regulated by regulators in regulator_names and fmembers
        try:
            gene_names,dbname = self.tfta.find_targets(regulator_names, fmembers)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        
        if len(gene_names) < 3:
            try:
                pathwayName, dblink = self.tfta.find_pathway_genes(gene_names)
            except PathwayNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
                return reply
            reply = _wrap_pathway_message(pathwayName, dblink, of_those=of_those)
            return reply
        else:
            try:
                pathwayName, dblink, genes = self.tfta.find_common_pathway_genes(gene_names)
            except PathwayNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
                return reply
            
            reply = self._wrap_pathway_genelist_message(pathwayName, dblink, genes, 
                    gene_descr='gene-list', of_those=of_those)
            return reply

    def find_pathway_gene_keyword(self,content):
        """
        Response content to find_pathway_gene_name request
        For a given gene list and keyword, reply the related pathways information
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name or keyword_name in ['pathway', 'signaling pathway']:
            reply = make_failure('NO_KEYWORD')
            return reply
            
        gene_names,fmembers = self._get_targets2(content, descr='gene')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
            
        try:
            pathwayName, dblink = \
                self.tfta.find_pathway_gene_keyword(gene_names, keyword_name, fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            #reply = self.wrap_family_message_pathway(term_id, descr='pathways', msg="PATHWAY_NOT_FOUND")
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=[keyword_name], of_those=of_those)
        return reply
        
    def find_pathway_keyword_regulator(self,content):
        """
        Response content to find_pathway_gene_name request
        For a given gene list and keyword, reply the related pathways information
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name or keyword_name in ['pathway', 'signaling pathway']:
            reply = make_failure('NO_KEYWORD')
            return reply
            
        regulator_names,fmembers = self._get_targets2(content, descr='regulator')
        if not regulator_names and not fmembers:
            reply = make_failure('NO_REGULATOR_NAME')
            return reply
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
            
        #get genes regulated by regulators in regulator_names
        try:
            gene_names,dbname = self.tfta.find_targets(regulator_names, fmembers)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        
        if len(gene_names) < 3:
            try:
                pathwayName, dblink = \
                    self.tfta.find_pathway_gene_keyword(gene_names, keyword_name)
            except PathwayNotFoundException:
                reply = self.KQMLList.from_string('(SUCCESS :pathways NIL)')
                return reply
            reply = _wrap_pathway_message(pathwayName, dblink, keyword=[keyword_name], of_those=of_those)
            return reply
        else:
            #go to find-common-pathway-genes
            try:
                pathwayName,dblink,genes = \
                   self.tfta.find_common_pathway_genes_keyword(gene_names, keyword_name)
            except PathwayNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
                return reply
            
            reply = self._wrap_pathway_genelist_message(pathwayName, dblink, genes, pathway_names=[keyword_name],
                                               gene_descr='gene-list', of_those=of_those)
            return reply

    def find_pathway_db_gene(self, content):
        """
        Response content to FIND-PATHWAY-DB-GENE request
        For a given gene list and certain db source, reply the related
        pathways information
        """
        db_name = _get_keyword_name(content, descr='database', low_case=False)
        if not db_name:
            reply = make_failure('NO_DB_NAME')
            return reply
            
        gene_names,fmembers = self._get_targets2(content, descr='gene')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
            
        try:
            pathwayName,dblink = self.tfta.find_pathway_db_gene(db_name,gene_names,fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink, of_those=of_those)
        return reply
        
    def find_pathway_db_regulator(self, content):
        """Response content to FIND-PATHWAY-DB-GENE request
        For a given gene list and certain db source, reply the related
        pathways information"""
        db_name = _get_keyword_name(content, descr='database', low_case=False)
        if not db_name:
            reply = make_failure('NO_DB_NAME')
            return reply
            
        regulator_names,fmembers = self._get_targets2(content, descr='regulator')
        if not regulator_names and not fmembers:
            reply = make_failure('NO_REGULATOR_NAME')
            return reply
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
        
        #get genes regulated by regulators in regulator_names
        try:
            gene_names,dbname = self.tfta.find_targets(regulator_names, fmembers)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        if len(gene_names) < 3:
            try:
                pathwayName,dblink = self.tfta.find_pathway_db_gene(db_name,gene_names)
            except PathwayNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
                return reply
            reply = _wrap_pathway_message(pathwayName, dblink, of_those=of_those)
            return reply
        else:
            try:
                pathwayName,dblink,genes = \
                   self.tfta.find_common_pathway_genes_db(gene_names, db_name)
            except PathwayNotFoundException:
                reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
                return reply
            
            reply = self._wrap_pathway_genelist_message(pathwayName, dblink, genes, 
                          gene_descr='gene-list', of_those=of_those)
            return reply

    def respond_find_tf_pathway(self, content):
        """
        Response content to FIND_TF_PATHWAY request
        For a given pathway name, reply the tfs within the pathway
        """
        logger.info('TFTA is processing the task: FIND-TF-PATHWAY.')
        pathway_names = self._get_pathway_name(content)
        if not pathway_names:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        
        #consider another optional parameter for subsequent query
        of_gene_names,nouse = self._get_targets(content, descr='of-those')
        
        try:
            pathwayName,tflist,dblink = self.tfta.find_tf_pathway(pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        logger.info('FIND-TF-PATHWAY: wrapping message...')
        reply = self._wrap_pathway_genelist_message(pathwayName, dblink, tflist, pathway_names=pathway_names,
                                               gene_descr='tfs', of_gene_names=of_gene_names)
        return reply
        
    def respond_find_kinase_pathway(self, content):
        """
        Response content to FIND_KINASE_PATHWAY request
        For a given pathway name, reply the kinases within the pathway
        """
        pathway_names = self._get_pathway_name(content)
        if not pathway_names:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
            
        #consider an optional parameter for subsequent query
        of_those_names,nouse = self._get_targets(content, descr='of-those')
        
        try:
            pathwayName, kinaselist, dblink = \
                self.tfta.find_kinase_pathway(pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        
        reply = self._wrap_pathway_genelist_message(pathwayName, dblink, kinaselist, pathway_names=pathway_names,
                                               gene_descr='kinase', of_gene_names=of_those_names)
        return reply

    def respond_find_gene_pathway(self, content):
        """
        Response content to FIND-GENE-PATHWAY request
        For a given pathway name, reply the genes within the pathway
        """
        logger.info('TFTA is processing the task: FIND-GENE-PATHWAY.')
        pathway_names = self._get_pathway_name(content)
        if not pathway_names:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
        try:
            pathwayName,genelist,plink = self.tfta.find_gene_pathway(pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        #consider an optional parameter for subsequent query
        of_gene_names,nouse = self._get_targets(content, descr='of-those')
        
        logger.info('FIND-GENE-PATHWAY: wrapping message...')
        reply = self._wrap_pathway_genelist_message(pathwayName, plink, genelist, pathway_names=pathway_names,
                                       gene_descr='genes', of_gene_names=of_gene_names)
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
            
        #Consider an optional parameter
        of_those = self._get_of_those_pathway(content)
        
        try:
            pathwayName, dblink = \
                self.tfta.find_pathway_keyword(keyword_name[0])
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=keyword_name, of_those=of_those)
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
        of_gene_names,nouse = self._get_targets(content, descr='of-those')
        
        reply = self._wrap_pathway_genelist_message(pathwayName, dblink, tflist, pathway_names=[keyword_name],
                                               gene_descr='tfs', of_gene_names=of_gene_names)
        return reply

    def respond_find_common_tf_genes(self, content):
        """
        Response content to FIND-COMMON-TF-GENES request
        For a given target list, reply the tfs regulating these genes
        and the regulated targets by each TF
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names or len(target_names) < 2:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
       
        #consider another parameter for subsequent query
        of_those_names,nouse = self._get_targets(content, descr='of-those')
        
        try:
            tf_targets = self.tfta.find_common_tfs(target_names)
        except TFNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
            return reply
        #cluster the tfs according to the targets
        tf_clustered = cluster_dict_by_value2(tf_targets)
        targets = list(tf_clustered.keys())
        tf_target_mes = []
        for tg in targets:
            if of_those_names:
                tfs = set(of_those_names) & set(tf_clustered[tg])
            else:
                tfs = tf_clustered[tg]
            if tfs:
                target_list = tg.split(',')
                target_json = self._get_genes_json(target_list)
                tf_json = self._get_genes_json(tfs)
                mes = KQMLList()
                mes.set('tf-list', tf_json)
                mes.set('target-list', target_json)
                tf_target_mes.append(mes.to_string())
        if tf_target_mes:
            mes_str = '(' + ' '.join(tf_target_mes) + ')'
            reply = KQMLList('SUCCESS')
            reply.set('tfs', mes_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        return reply
        
    def get_common_pathway_genes(self, content):
        """
        response content to FIND-COMMON-PATHWAY-GENES request
        """
        gene_names,fmembers = self._get_targets2(content, descr='target')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        #logger.info('FIND-COMMON-PATHWAY-GENES:genes=' + str(gene_names))
        #logger.info('FIND-COMMON-PATHWAY-GENES:fmembers=' + str(fmembers))
        
        try:
            pathwayName, dblink, genes = self.tfta.find_common_pathway_genes(gene_names, fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
        
        #logger.debug('FIND-COMMON-PATHWAY-GENES: starting wrapping message.')
        
        reply = self._wrap_pathway_genelist_message(pathwayName, dblink, genes, gene_descr='gene-list')
        
        logger.debug('FIND-COMMON-PATHWAY-GENES: sending message...')
        return reply

    def get_common_pathway_genes_keyword(self, content):
        """
        Response content to FIND-COMMON-PATHWAY-GENES-KEYWORD request
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name or keyword_name in ['pathway', 'signaling pathway']:
            reply = make_failure('NO_KEYWORD')
            return reply
        else:
            keyword_name = trim_word([keyword_name], 'pathway')
        
        gene_names,fmembers = self._get_targets2(content, descr='target')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        try:
            pathwayName,dblink,genes = self.tfta.find_common_pathway_genes_keyword(gene_names, keyword_name[0], fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = self._wrap_pathway_genelist_message(pathwayName, dblink, genes, pathway_names=keyword_name,
                                               gene_descr='gene-list')
        return reply
        
    def get_common_pathway_genes_db(self, content):
        """
        Response content to FIND-COMMON-PATHWAY-GENES-DB request
        """
        db_name = _get_keyword_name(content, descr='database', low_case=False)
        if not db_name:
            reply = make_failure('NO_DB_NAME')
            return reply
            
        gene_names,fmembers = self._get_targets2(content, descr='target')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        try:
            pathwayName,dblink,genes = \
               self.tfta.find_common_pathway_genes_db(gene_names, db_name, fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = self._wrap_pathway_genelist_message(pathwayName, dblink, genes, gene_descr='gene-list')
        return reply

    def respond_is_pathway_gene(self, content):
        """
        Respond to IS-PATHWAY-GENE request
        query like: Does the mTor pathway utilize SGK1? 
        """
        pathway_names = self._get_pathway_name(content)
        if not pathway_names:
            reply = make_failure('NO_PATHWAY_NAME')
            return reply
            
        gene_names,fmembers = self._get_targets2(content, descr='gene')
        if not gene_names and not fmembers:
            reply = make_failure('NO_GENE_NAME')
            return reply
            
        try:
            pathwayName, dblink = self.tfta.Is_pathway_gene(pathway_names, gene_names, fmembers)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        logger.info('IS-PATHWAY-GENE: wrapping message...')
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
            
        tf_names,term_id = self._get_targets(content, descr='tf')
        if not tf_names:
            #tf_arg = content.gets('tf')
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
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
        
        tf_names,term_id = self._get_targets(content, descr='tf')
        if not tf_names:
            #tf_arg = content.gets('tf')
            reply = self.wrap_family_message(term_id, 'NO_TF_NAME')
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
        miRNA_name = self._get_mirnas(content)
        if not miRNA_name:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        
        target_name,term_id = self._get_targets(content, descr='target')
        if not target_name:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        
        strength_name = _get_keyword_name(content, descr='strength')
        
        if strength_name:
            is_target,expr,supt,pmid,miRNA_mis = self.tfta.Is_miRNA_target_strength(miRNA_name, target_name[0], strength_name)
        else:
            is_target,expr,supt,pmid,miRNA_mis = self.tfta.Is_miRNA_target(miRNA_name, target_name[0])
        
        #provenance support
        self.send_background_support_mirna(list(miRNA_name.keys())[0], target_name[0], expr, supt, pmid, strength=strength_name)
        
        #respond to BA
        #check if it's necessary for user clarification
        if miRNA_mis:
            reply = self._get_mirna_clarification(miRNA_mis)
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
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
        
        strength_name = _get_keyword_name(content, descr='strength')
        
        ##consider another parameter for subsequent query
        of_those_names = self._get_mirnas(content, descr='of-those')
        
        try:
            if strength_name:
                miRNA_names,expr,supt,pmid = self.tfta.find_miRNA_target_strength(target_names, strength_name)
            else:
                miRNA_names,expr,supt,pmid = self.tfta.find_miRNA_target(target_names)
        except TargetNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)') 
            return reply
            
        if of_those_names:
            res_db = set(','.join(miRNA_names).upper().split(','))
            res_of = set(','.join(of_those_names).upper().split(','))
            miRNA_names = res_db & res_of
        
        #provenance support
        self.send_background_support_mirna(miRNA_names, target_names, expr, supt, pmid, strength=strength_name, find_mirna=True)
        
        #respond to BA
        if len(miRNA_names):
            mir_agent = [Agent(mir, db_refs={'type':'MIRNA'}) for mir in miRNA_names]
            mir_json = self.make_cljson(mir_agent)
            reply = KQMLList('SUCCESS')
            reply.set('miRNAs', mir_json)
        else:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)')
        return reply
        
    def respond_find_target_miRNA(self, content):
        """
        Respond to FIND-TARGET-MIRNA request
        """
        miRNA_names = self._get_mirnas(content)
        if not miRNA_names:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
            
        strength_name = _get_keyword_name(content, descr='strength')
            
        #consider another parameter for subsequent query
        of_those_names,nouse = self._get_targets(content, descr='of-those')
        
        #consider an optional parameter to only return given type of genes
        target_type_set,target_type = self.get_target_type_set(content)
        
        if strength_name:
            target_names,miRNA_mis,expr,supt,pmid = self.tfta.find_target_miRNA_strength(miRNA_names, strength_name)
        else:
            target_names,miRNA_mis,expr,supt,pmid = self.tfta.find_target_miRNA(miRNA_names)
        #check if it's necessary for user clarification
        if miRNA_mis:
            reply = self._get_mirna_clarification(miRNA_mis)
            return reply
        else:
            if of_those_names:
                target_names = set(of_those_names) & target_names
            #check if it requires returning a type of genes
            if target_type_set:
                target_names = target_type_set & target_names
            
            #provenance message
            self.send_background_support_mirna(list(miRNA_names.keys()), target_names, expr, supt, pmid, strength=strength_name, find_target=True)
            
            #Response to BA
            if len(target_names):
                reply = self.wrap_message('targets', target_names)
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply
        
    def respond_find_evidence_miRNA_target(self, content):
        """
        Respond to FIND-EVIDENCE-MIRNA-TARGET request
        """
        miRNA_name = self._get_mirnas(content)
        target_name,term_id = self._get_targets(content, descr='target')
        strength_name = _get_keyword_name(content, descr='strength')
        
        if all([miRNA_name, target_name, strength_name]):
            expe,sType,pmlink,miRNA_mis = self.tfta.find_evidence_miRNA_target_strength(miRNA_name, target_name[0], strength_name)
        elif all([miRNA_name, target_name]):
            expe,sType,pmlink,miRNA_mis = self.tfta.find_evidence_miRNA_target(miRNA_name, target_name[0])
        elif target_name:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        else:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        
        if miRNA_mis:
            reply = self._get_mirna_clarification(miRNA_mis)
            return reply
        
        mes_list = []
        if expe:
            for e,s,l in zip(expe,sType,pmlink):
                en = '"' + e + '"'
                sn = '"' + s + '"'
                ln = '"' + l + '"'
                mes = KQMLList()
                mes.set('experiment', en)
                mes.set('supportType', sn)
                mes.set('pubmedLink', ln)
                mes_list.append(mes.to_string())
        if mes_list:
            evi_str = '(' + ' '.join(mes_list) + ')'
            reply = KQMLList('SUCCESS')
            reply.set('evidence', evi_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)')   
        return reply
        
    def respond_find_evidence_miRNA_exp(self, content):
        """
        Respond to FIND-EVIDENCE-MIRNA-EXP request
        """
        strength = _get_keyword_name(content, descr='strength')
        res = self.tfta.find_evidence_strength(strength)
        if res:
            reply = KQMLList('SUCCESS')
            exp_str = ''
            for e in res:
                exp_str += '(:name %s )' % e
                exp_str = '( ' + exp_str + ') '
            reply.set('evidence', exp_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)') 
        return reply

    def respond_find_gene_count_miRNA(self, content):
        """
        Respond to FIND-GENE-COUNT-MIRNA request
        """
        miRNA_names = self._get_mirnas(content)
        if not miRNA_names:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        
        #consider another parameter for subsequent query
        of_those_names,nouse = self._get_targets(content, descr='of-those')
        
        target_count,target_mirna,miRNA_mis = self.tfta.find_gene_count_miRNA(miRNA_names, of_those=of_those_names)
        if not target_count:
            #clarification
            if miRNA_mis:
                reply = self._get_mirna_clarification(miRNA_mis)
                return reply
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
                return reply
        
        target_mir_mes = []
        if target_count:
            for t,c in target_count.items():
                mes = KQMLList()
                mes.set('target', self.make_cljson(Agent(t, db_refs={'TYPE':'ONT::GENE-PROTEIN'})))
                mes.set('count', str(c))
                mes.set('miRNA', self.make_cljson([Agent(mir, db_refs={'TYPE':'MIRNA'}) for mir in target_mirna[t]]))
                target_mir_mes.append(mes.to_string())
        if target_mir_mes:
            mes_str = '(' + ' '.join(target_mir_mes) + ')'
            reply = KQMLList('SUCCESS')
            reply.set('targets', mes_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply
             
    def respond_find_miRNA_count_gene(self, content):
        """
        Respond to FIND-MIRNA-COUNT-GENE request
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names or len(target_names) < 2:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        #consider another parameter for subsequent query
        of_those_names = self._get_mirnas(content, descr='of-those')
        
        mirna_count,mir_targets = self.tfta.find_miRNA_count_gene(target_names, of_those=of_those_names)
        
        mir_target_mes = []
        if mirna_count:
            for mir,count in mirna_count.items():
                mes = KQMLList()
                mes.set('miRNA', self.make_cljson(Agent(mir, db_refs={'type': 'MIRNA'})))
                mes.set('count', str(count))
                mes.set('target', self.make_cljson([Agent(t, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for t in mir_targets[mir]]))
                mir_target_mes.append(mes.to_string())
        if mir_target_mes:
            mes_str = '(' + ' '.join(mir_target_mes) + ')'
            reply = KQMLList('SUCCESS')
            reply.set('miRNAs', mes_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)')
        return reply
        
    def respond_is_mirna_disease(self, content):
        """
        Respond to IS-MIRNA-DISEASE request
        """
        miRNA_name = self._get_mirnas(content)
        if not miRNA_name:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        
        disease = _get_keyword_name(content, descr='disease', hyphen=True)
        if not disease:
            reply = self.wrap_family_message(term_id, 'NO_DISEASE_NAME')
            return reply
        
        res = self.md.is_mirna_disease(miRNA_name, disease)
        
        reply = KQMLList('SUCCESS')
        res_str = 'TRUE' if res else 'FALSE'
        reply.set('result', res_str)
        return reply
        
    def respond_find_mirna_disease(self, content):
        """
        Respond to find-mirna-disease request
        """
        disease = _get_keyword_name(content, descr='disease', hyphen=True)
        if not disease:
            reply = self.wrap_family_message(term_id, 'NO_DISEASE_NAME')
            return reply
            
        mirnas = self.md.find_mirna_disease(disease)
        
        ##consider another parameter for subsequent query
        of_those_names = self._get_mirnas(content, descr='of-those')
        
        if of_those_names:
            res_db = set(','.join(mirnas).upper().split(','))
            res_of = set(','.join(of_those_names).upper().split(','))
            mirnas = res_db & res_of
            
        if len(mirnas):
            #mir_agent = [Agent(mir, db_refs={'type':'MIRNA'}) for mir in miRNA_names]
            #mir_json = self.make_cljson(mir_agent)
            mir_json = self._to_json_mirna(mirnas)
            reply = KQMLList('SUCCESS')
            reply.set('miRNAs', mir_json)
        else:
            reply = KQMLList.from_string('(SUCCESS :miRNAs NIL)')
        return reply
        
    def respond_find_disease_mirna(self, content):
        """
        Respond to fidn-disease-mirna request.
        """
        miRNA_names = self._get_mirnas(content)
        if not miRNA_names:
            #return all the diseases
            diseases = self.md.find_disease()
        else:
            diseases = self.md.find_disease_mirna(miRNA_names)
            
        if diseases:
            dagent = [Agent(d, db_refs={'type':'disease'}) for d in diseases]
            djson = self.make_cljson(dagent)
            reply = KQMLList('SUCCESS')
            reply.set('disease', djson)
        else:
            reply = KQMLList.from_string('(SUCCESS :disease NIL)')
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
            pathwayName,dblink = self.tfta.find_pathway_db_keyword(db_name, pathway_names)
        except PathwayNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
            return reply
            
        reply = _wrap_pathway_message(pathwayName, dblink, keyword=pathway_names)
        return reply

    def respond_find_tissue_gene(self, content):
        """
        Response to FIND-TISSUE-GENE task
        """
        gene_arg = content.get('gene')
        #without gene parameter, then return all the tissues
        if not gene_arg:
            tissue_agent = [Agent(ts) for ts in self.tissue_list]
            tissue_json = self.make_cljson(tissue_agent)
            reply = KQMLList('SUCCESS')
            reply.set('tissue', tissue_json)
            return reply
            
        gene_name,term_id = self._get_targets(content, descr='gene')
        if not gene_name:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
            
        try:
            tissues = self.tfta.find_tissue_gene(gene_name[0])
        except TissueNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :tissue NIL)')
            return reply
            
        if tissues:
            tissue_agent = [Agent(ts, db_refs={'TYPE':'tissue'}) for ts in tissues]
            tissue_json = self.make_cljson(tissue_agent)
            reply = KQMLList('SUCCESS')
            reply.set('tissue', tissue_json)
        else:
            reply = KQMLList.from_string('(SUCCESS :tissue NIL)')
        return reply
        
    def get_target_kinase(self, kinase_names, keyword=None, of_those=None, target_type=None, out_gene=True):
        """
        Response to find-target with kinase parameter
        """
        if keyword in ['increase', 'decrease']:
            try:
                targets = self.tfta.find_target_kinase_keyword(kinase_names, keyword)
            except TargetNotFoundException:
                targets = set()
        else:
            try:
                targets = self.tfta.find_target_kinase(kinase_names)
            except TargetNotFoundException:
                targets = set()
        
        if targets and of_those:
            targets = set(targets).intersection(set(of_those))
        if targets and target_type:
            targets = set(targets).intersection(target_type)
        
        if out_gene:
            return targets
            
        if targets:
            reply = KQMLList('SUCCESS')
            target_json = self._get_genes_json(targets)
            tmess = KQMLList()
            tmess.set('target-db', target_json)
            reply.set('targets', tmess)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply
        
    def respond_find_kinase_target(self, content):
        """
        Response to find-kinase-regulation with only target parameter
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
            
        try:
            kinases = self.tfta.find_kinase_target(target_names)
        except KinaseNotFoundException:
            reply = KQMLList.from_string('(SUCCESS :kinase NIL)')
            return reply
        
        if kinases:
            kin_json = self._get_genes_json(kinases)
            reply = KQMLList('SUCCESS')
            reply.set('kinase', kin_json)
        else:
            reply = KQMLList.from_string('(SUCCESS :kinase NIL)')
        return reply
        
    def respond_find_kinase_target_keyword(self, content):
        """
        Response to find-kinase-regulation with target and keyword parameter
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
        if kinases:
            kin_json = self._get_genes_json(kinases)
            reply = KQMLList('SUCCESS')
            reply.set('kinase', kin_json)
        else:
            reply = KQMLList.from_string('(SUCCESS :kinase NIL)')
        return reply
        
    def respond_find_kinase_regulation(self, content):
        """
        Response to FIND-KINASE-REGULATION
        """
        target_arg = content.get('target')
        keyword_arg = content.get('keyword')
        if all([target_arg, keyword_arg]):
            reply = self.respond_find_kinase_target_keyword(content)
        else:
            reply = self.respond_find_kinase_target(content)
        
        return reply
        
    def respond_find_tf_miRNA(self, content):
        """
        Response to FIND-TF-MIRNA
        """
        miRNA_names = self._get_mirnas(content)
        if not miRNA_names:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
            
        #consider an additional parameter for subsequent query
        of_those_names,nouse = self._get_targets(content, descr='of-those')
        
        tf_names,miRNA_mis = self.tfta.find_tf_miRNA(miRNA_names)
        #check if it's necessary for user clarification
        if miRNA_mis:
            reply = self._get_mirna_clarification(miRNA_mis)
        else:
            if of_those_names:
                tf_names = set(of_those_names) & set(tf_names)
            if tf_names:
                reply = self.wrap_message(':tfs ', tf_names)
            else:
                reply = KQMLList.from_string('(SUCCESS :tfs NIL)')
        return reply
        
    def find_target_source(self, content):
        """
        Response to find-target with :source parameter
        """
        reg_names,term_id = self._get_targets(content, descr='regulator')
        if not reg_names:
            reply = self.wrap_family_message(term_id, 'NO_REGULATOR_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
        reg_names = list(set(reg_names))
        
        keyword_name = _get_keyword_name(content)
        if not keyword_name:
            keyword_name = 'regulate'
            
        source_name = _get_keyword_name(content, descr='source')
        if not source_name:
            reply = make_failure('INVALID_SOURCE')
            return reply
            
        #consider an optional parameter for subsequent query
        of_those,_ = self._get_targets(content, 'of-those')
        #consider an optional parameter to only return given type of genes
        target_type_set,target_type = self.get_target_type_set(content)
        
        if source_name == 'literature':
            try:
                stmt_types = stmt_type_map[keyword_name]
            except KeyError:
                reply = make_failure('INVALID_KEYWORD')
                return reply
            lit_json = self.get_target_indra(reg_names, stmt_types, keyword_name,
                                             of_those=of_those, target_type=target_type_set)
            if len(lit_json):
                reply = KQMLList('SUCCESS')
                tmess = KQMLList()
                tmess.set('target-literature', lit_json)
                reply.set('targets', tmess)
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        elif source_name == 'kinase':
            reply = self.get_target_kinase(reg_names, keyword=keyword_name, of_those=of_those,
                                             target_type=target_type_set, out_gene=False)
        else:
            #by default go to tf-db
            reply = self.find_target_tf(reg_names, of_those=of_those, target_type=(target_type_set, target_type),
                                          out_gene=False)
        return reply
    
    def find_target_all(self, content):
        """
        Response to find-target with results from all sources
        By default, it will only return results from our dbs. Only if :literature is set to true, 
        it will consider literature search in order to avoid long delay caused by literature search.
        """
        reg_names,term_id = self._get_targets(content, descr='regulator')
        if not reg_names:
            reply = self.wrap_family_message(term_id, 'NO_REGULATOR_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
        reg_names = list(set(reg_names))
        
        keyword_name = _get_keyword_name(content)
        if not keyword_name:
            keyword_name = 'regulate'
            
        #consider an optional parameter for subsequent query
        of_those,_ = self._get_targets(content, 'of-those')
        #consider an optional parameter to only return given type of genes
        target_type_set,target_type = self.get_target_type_set(content)
        
        #results from tf-db
        if keyword_name == 'regulate':
            targets_tfdb = self.find_target_tf(reg_names, of_those=of_those, target_type=(target_type_set,target_type))
        else:
            targets_tfdb = set()
        
        #results from kinase-db
        targets_kdb = self.get_target_kinase(reg_names, keyword=keyword_name, of_those=of_those,
                                                   target_type=target_type_set, out_gene=True)
        
        #results from literature
        lit = _get_keyword_name(content, descr='literature')
        if lit and lit.upper() == 'TRUE':
            try:
                stmt_types = stmt_type_map[keyword_name]
            except KeyError:
                stmt_types = stmt_type_map['regulate']
            lit_json = self.get_target_indra(reg_names, stmt_types, keyword_name,
                                        of_those=of_those, target_type=target_type_set)
        else:
            lit_json = []
        
        targets = set(targets_tfdb).union(set(targets_kdb))
        if targets or lit_json:
            reply = KQMLList('SUCCESS')
            tmess = KQMLList()
            if targets:
                target_json = self._get_genes_json(targets)
                tmess.set('target-db', target_json)
            if lit_json:
                tmess.set('target-literature', lit_json)
            reply.set('targets', tmess)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        
        return reply
    
    def find_target_tf(self, tf_names, of_those=None, target_type=None, out_gene=True):
        target_names = set()
        try:
            target_names,dbname = self.tfta.find_targets(tf_names)
        except TFNotFoundException:
            #provenance support
            self.send_background_support_db(tf_names, [], '')
            if out_gene:
                return target_names
            else:
                reply = KQMLList.from_string('(SUCCESS :targets NIL)')
                return reply
    
        #check if it's a subsequent query
        if of_those:
            target_names = set(of_those) & set(target_names)
        
        #check if it requires returning a type of genes
        if target_type[0]:
            target_names = target_type[0] & set(target_names)
            
        #provenance support
        self.send_background_support_db(tf_names, target_names, dbname, find_target=True, target_type=target_type[1])
        
        if out_gene:
            return target_names
            
        if target_names:
            reply = KQMLList('SUCCESS')
            target_json = self._get_genes_json(target_names)
            tmess = KQMLList()
            tmess.set('target-db', target_json)
            reply.set('targets', tmess)
        else:
            reply = KQMLList.from_string('(SUCCESS :targets NIL)')
        return reply
    
    def find_regulation_source(self, content):
        """
        Respond to find-regulation
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
        of_those_names,nouse = self._get_targets(content, 'of-those')
            
        if source_name == 'literature':
            #literature result
            try:
                stmt_types = stmt_type_map[keyword_name]
            except KeyError as e:
                reply = make_failure('INVALID_KEYWORD')
                return reply
            
            tf_json, gene_json, mirna_json, oth_json = [],[],[],[]
            if of_those_names:
                gene_json = self.get_regulator_indra_those(target_names, stmt_types, keyword_name,
                                                        of_those=of_those_names)
            else:
                tf_json, gene_json, mirna_json, oth_json = self.get_regulator_indra(target_names, stmt_types, keyword_name)
            
            descr_list = ['gene-literature','tf-literature','mirna-literature','other-literature']
            json_list = [gene_json, tf_json, mirna_json, oth_json]
            res = self._combine_json_list(descr_list, json_list)
            if res:
                reply = KQMLList('SUCCESS')
                reply.set('regulators', res)
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        elif source_name == 'geo rnai':
            #kinase regualtion
            kin_json = self.get_kinase_regulation(target_names, keyword_name, of_those=of_those_names)
            
            if len(kin_json):
                mes = KQMLList()
                mes.set('kinase-db', kin_json)
                reply = KQMLList('SUCCESS')
                reply.set('regulators', mes)
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        else:
            reply = make_failure('INVALID_SOURCE')
        return reply
        
    def find_regulation_agent(self, content):
        """
        Response to find-regulation request
        For example: Which kinases regulate the cfos gene?
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
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
        of_those_names,nouse = self._get_targets(content, descr='of-those')
            
        if agent_name == 'kinase':
            #kinase regualtion
            kin_json = self.get_kinase_regulation(target_names, keyword_name, of_those=of_those_names)
            if len(kin_json):
                mes = KQMLList()
                mes.set('kinase-db', kin_json)
                reply = KQMLList('SUCCESS')
                reply.set('regulators', mes)
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        elif agent_name in ['transcription factor', 'tf']:
            temp_reply = self.respond_find_target_tfs(content)
            tf_str = temp_reply.get('tfs')
            if tf_str != 'NIL':
                mes = KQMLList()
                mes.set('tf-db', tf_str)
                reply = KQMLList('SUCCESS')
                reply.set('regulators', mes)
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
        else:
            reply = make_failure('INVALID_REGULATOR')
        return reply
        
    def find_regulation_all(self, content):
        """
        Response to find-regulation request
        For example: what regulate MYC?
        """
        target_names,term_id = self._get_targets(content, descr='target')
        if not target_names:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
        if term_id:
            reply = self.wrap_family_message(term_id, 'FAMILY_NAME')
            return reply
        target_names = list(set(target_names))
        
        keyword_name = _get_keyword_name(content, descr='keyword')
        if not keyword_name:
            reply = make_failure('NO_KEYWORD')
            return reply
        
        #consider an optional parameter for subsequent query
        of_those_names,nouse = self._get_targets(content, 'of-those')
        
        #literature result
        try:
            stmt_types = stmt_type_map[keyword_name.lower()]
        except KeyError as e:
            reply = make_failure('INVALID_KEYWORD')
            return reply
        
        tf_json, gene_json, mirna_json, oth_json = [],[],[],[]
        if of_those_names:
            gene_json = self.get_regulator_indra_those(target_names, stmt_types, keyword_name,
                                                          of_those=set(of_those_names))
        else:
            tf_json, gene_json, mirna_json, oth_json = self.get_regulator_indra(target_names, stmt_types, keyword_name)
        
        #kinase regulation
        kin_json = self.get_kinase_regulation(target_names, keyword_name, of_those=of_those_names)
        
        #db result
        #only considering the regulation case
        if keyword_name == 'regulate':
            try:
                tf_names,dbname = self.tfta.find_tfs(target_names)
                if of_those_names:
                    tf_names = list(set(tf_names).intersection(set(of_those_names)))
                #provenance support
                self.send_background_support_db(tf_names, target_names, dbname, find_tf=True)
            except TargetNotFoundException:
                #provenance support
                self.send_background_support_db([], target_names, '')
                
                descr_list = ['kinase-db','gene-literature','tf-literature','mirna-literature','other-literature']
                json_list = [kin_json, gene_json, tf_json, mirna_json, oth_json]
                res = self._combine_json_list(descr_list, json_list)
                if res:
                    reply = KQMLList('SUCCESS')
                    reply.set('regulators', res)
                else:
                    reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
                return reply
        else:
            tf_names = []
        if len(tf_names):
            tf_db_json = self._get_genes_json(tf_names)
            descr_list = ['tf-db','kinase-db','gene-literature','tf-literature','mirna-literature','other-literature']
            json_list = [tf_db_json, kin_json, gene_json, tf_json, mirna_json, oth_json]
            res = self._combine_json_list(descr_list, json_list)
            reply = KQMLList('SUCCESS')
            reply.set('regulators', res)
        else:
            descr_list = ['kinase-db','gene-literature','tf-literature','mirna-literature','other-literature']
            json_list = [kin_json, gene_json, tf_json, mirna_json, oth_json]
            res = self._combine_json_list(descr_list, json_list)
            if res:
                reply = KQMLList('SUCCESS')
                reply.set('regulators', res)
            else:
                reply = KQMLList.from_string('(SUCCESS :regulators NIL)')
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
        regulator_name,term_id = self._get_targets(content, descr='regulator')
        if not regulator_name:
            reply = self.wrap_family_message(term_id, 'NO_REGULATOR_NAME')
            return reply
        
        target_name,term_id = self._get_targets(content, descr='target')
        if not target_name:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
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
        term_tuple = (regulator_name[0], target_name[0], keyword_name)
        if term_tuple not in self.stmts_indra:
            stmts,success = self.tfta.find_statement_indraDB(subj=regulator_name[0], obj=target_name[0], stmt_types=stmt_types)
            if success:
                self.stmts_indra[term_tuple] = stmts
        else:
            stmts = self.stmts_indra[term_tuple]
        self.send_background_support(stmts, regulator_name[0], target_name[0], keyword_name)
        if len(stmts):
            evidences = self.tfta.find_evidence_indra(stmts)
            if len(evidences):
                evi_message = _wrap_evidence_message(':literature', evidences)
        
        if keyword_name == 'regulate':
            db_names = self.tfta.find_evidence_dbname(regulator_name[0], target_name[0])
            if len(db_names):
                db_message = _wrap_dbname_message(':tf-db', db_names)
                
        message = _combine_messages([db_message, evi_message])
        if message:
            reply = KQMLList('SUCCESS')
            reply.set('evidence', '(' + message + ')')
        else:
            reply = KQMLList.from_string('(SUCCESS :evidence NIL)')
        return reply
    
    def respond_find_evidence_tfdb(self, content):
        """
        Respond to find-evidence request from tfdb
        """
        regulator_name,term_id = self._get_targets(content, descr='regulator')
        if not regulator_name:
            reply = self.wrap_family_message(term_id, 'NO_REGULATOR_NAME')
            return reply
        
        target_name,term_id = self._get_targets(content, descr='target')
        if not target_name:
            reply = self.wrap_family_message(term_id, 'NO_TARGET_NAME')
            return reply
            
        db_names = self.tfta.find_evidence_dbname(regulator_name[0], target_name[0])
        if len(db_names):
            db_message = _wrap_dbname_message(':tf-db', db_names)
            reply = KQMLList('SUCCESS')
            reply.set('evidence', '( ' + db_message + ')')
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
        if keyword_name:
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
        ogenes,nouse = self._get_targets(content, descr='gene')
        if ogenes:
            results = list(set(gene_list) & set(ogenes))
        else:
            results = gene_list
        if len(results):
            gene_json = self._get_genes_json(results)
            reply = KQMLList('SUCCESS')
            reply.set('genes', gene_json)
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
            
        gene_name,term_id = self._get_targets(content, descr='gene')
        if not gene_name:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
            
        #optional keyword
        keyword_name = _get_keyword_name(content, descr='keyword')
        if keyword_name == 'exclusive':
            try:
                is_express = self.tfta.is_tissue_gene_exclusive(tissue_name, gene_name[0])
            except Exception:
                reply = KQMLList.from_string('(SUCCESS :result FALSE)')
                return reply
        else:
            try:
                is_express = self.tfta.is_tissue_gene(tissue_name, gene_name[0])
            except Exception:
                reply = KQMLList.from_string('(SUCCESS :result FALSE)')
                return reply
        reply = KQMLList('SUCCESS')
        is_express_str = 'TRUE' if is_express else 'FALSE'
        reply.set('result', is_express_str)
        return reply
        
    def respond_go_enrichment(self, content):
        """
        Respond to GO-ENRICHMENT
        """
        gene_names,term_id = self._get_targets(content, descr='gene')
        if not gene_names:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
            
        results = self.go.go_enrichment_analysis(gene_names)
        #return GO_id, GO_name, p_bonferroni, study_items
        mes_json = []
        if results:
            for res in results:
                mes = KQMLList()
                mes.sets('GO-term', res.goterm.id)
                mes.sets('GO-name', res.name)
                mes.sets('p-bonferroni', str(res.p_bonferroni))
                r1 = res.ratio_in_study
                mes.sets('ratio_in_study', str(r1[0]) + '/' + str(r1[1]))
                r1 = res.ratio_in_pop
                mes.sets('ratio_in_pop', str(r1[0]) + '/' + str(r1[1]))
                gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in res.study_items]
                gene_json = self.make_cljson(gene_agent)
                mes.set('genes', gene_json)
                mes_json.append(mes.to_string())
            reply=KQMLList('SUCCESS')
            res_str = '(' + ' '.join(mes_json) + ')'
            reply.set('results', res_str)
            return reply
        else:
            reply = KQMLList.from_string('(SUCCESS :results NIL)')
            return reply
    
    def respond_go_annotation(self, content):
        """
        Respond to GO-ANNOTATION request
        """
        keyword_name = _get_keyword_name(content, hyphen=True)
        if not keyword_name:
            reply = make_failure('NO_KEYWORD_NAME')
            return reply
        
        results = self.go.get_annotations_keyword(keyword_name)
        if results:
            #return GO_ID, gene_name
            agent = []
            for key, value in results.items():
                go_terms = ','.join([entry['GO_ID'] for entry in value])
                agent.append(Agent(key, db_refs={'GO':go_terms}))
            a_json = self.make_cljson(agent)
            reply = KQMLList('SUCCESS')
            reply.set('genes', a_json)
            return reply
        else:
            reply = KQMLList.from_string('(SUCCESS :genes NIL)')
            return reply
            
    def respond_pathway_enrichment(self, content):
        """
        Respond to pathway-enrichment request
        """
        gene_names,term_id = self._get_targets(content, descr='gene')
        if not gene_names:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
        db_name = _get_keyword_name(content, descr='database', low_case=True)
        if not db_name:
            db_name = 'kegg'
            
        results = self.pw.get_ora_pathway(gene_names, db_name)
        mes_json = []
        if results:
            for res in results.keys():
                mes = KQMLList()
                mes.sets('name', res)
                mes.sets('dblink', results[res]['dblink'])
                mes.sets('p-bonferroni', str(results[res]['p-bonferroni']))
                gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in results[res]['gene']]
                gene_json = self.make_cljson(gene_agent)
                mes.set('genes', gene_json)
                mes_json.append(mes.to_string())
            reply=KQMLList('SUCCESS')
            res_str = '(' + ' '.join(mes_json) + ')'
            reply.set('results', res_str)
            return reply
        else:
            reply = KQMLList.from_string('(SUCCESS :results NIL)')
            return reply
            
    def respond_disease_enrichment(self, content):
        """
        Respond to disease-enrichment request
        """
        gene_names,term_id = self._get_targets(content, descr='gene')
        if not gene_names:
            reply = self.wrap_family_message(term_id, 'NO_GENE_NAME')
            return reply
        db_name = _get_keyword_name(content, descr='database', low_case=True)
        if not db_name:
            db_name = 'ctd'
            
        results = self.pw.get_ora_disease(gene_names, db_name)
        mes_json = []
        if results:
            for res in results.keys():
                mes = KQMLList()
                mes.sets('name', res)
                mes.sets('dblink', results[res]['dblink'])
                mes.sets('p-bonferroni', str(results[res]['p-bonferroni']))
                gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in results[res]['gene']]
                gene_json = self.make_cljson(gene_agent)
                mes.set('genes', gene_json)
                mes_json.append(mes.to_string())
            reply=KQMLList('SUCCESS')
            res_str = '(' + ' '.join(mes_json) + ')'
            reply.set('results', res_str)
            return reply
        else:
            reply = KQMLList.from_string('(SUCCESS :results NIL)')
            return reply
            
    def respond_mirna_disease_enrichment(self, content):
        """
        Respond to mirna-disease-enrichment request
        """
        miRNA_names = self._get_mirnas(content)
        if not miRNA_names:
            reply = make_failure('NO_MIRNA_NAME')
            return reply
        db_name = 'hmdd'
        results = self.pw.get_ora_disease(miRNA_names, db_str=db_name)
        mes_json = []
        if results:
            for res in results.keys():
                mes = KQMLList()
                mes.sets('name', res)
                mes.sets('dblink', results[res]['dblink'])
                mes.sets('p-bonferroni', str(results[res]['p-bonferroni']))
                gene_agent = [Agent(g, db_refs={'TYPE':'MIRNA'}) for g in results[res]['gene']]
                gene_json = self.make_cljson(gene_agent)
                mes.set('mirnas', gene_json)
                mes_json.append(mes.to_string())
            reply=KQMLList('SUCCESS')
            res_str = '(' + ' '.join(mes_json) + ')'
            reply.set('results', res_str)
            return reply
        else:
            reply = KQMLList.from_string('(SUCCESS :results NIL)')
            return reply
    
    def respond_make_heatmap(self, content):
        """
        Respond to make-heatmap request
        """
        path = _get_keyword_name(content, descr='filepath')
        if not path:
            reply = make_failure('NO_FILE_NAME')
            return reply
            
        heatmap_file, row_index, col_index = generate_heatmap(path.lower())
        if heatmap_file:
            reply = KQMLList('SUCCESS')
            reply.sets('heatmap-file', heatmap_file)
            reply.sets('row-reordered-ind', ','.join([str(r) for r in row_index]))
            reply.sets('col-reordered-ind', ','.join([str(r) for r in col_index]))
        else:
            reply = make_failure('INVALID_DATA_FORMAT')
        return reply
        
    task_func = {'IS-REGULATION':respond_is_regulation, 'FIND-TF':respond_find_tf,
                 'FIND-PATHWAY':respond_find_pathway, 'FIND-TARGET':respond_find_target,
                 'FIND-COMMON-PATHWAY-GENES':respond_find_common_pathway_genes,
                 'FIND-MIRNA-COUNT-GENE':respond_find_miRNA_count_gene,
                 'FIND-GENE-COUNT-MIRNA':respond_find_gene_count_miRNA,
                 'IS-MIRNA-TARGET':respond_is_miRNA_target,
                 'FIND-COMMON-TF-GENES':respond_find_common_tf_genes,
                 'IS-PATHWAY-GENE':respond_is_pathway_gene,
                 'FIND-TARGET-MIRNA':respond_find_target_miRNA,
                 'FIND-MIRNA-TARGET':respond_find_miRNA_target,
                 'IS-GENE-ONTO':respond_is_gene_onto,
                 'FIND-GENE-ONTO':respond_find_gene_onto,
                 'IS-TF-TARGET':respond_is_tf_target,
                 'FIND-TF-TARGET':respond_find_tf_targets,
                 'FIND-TARGET-TF':respond_find_target_tfs,
                 'FIND-TF-PATHWAY':respond_find_tf_pathway,
                 'FIND-GENE-PATHWAY':respond_find_gene_pathway,
                 'FIND-PATHWAY-KEYWORD':respond_find_pathway_keyword,
                 'FIND-TF-KEYWORD':respond_find_tf_keyword,
                 'FIND-GENE-GO-TF':respond_find_genes_go_tf2,
                 'IS-TF-TARGET-TISSUE':respond_is_tf_target_tissue,
                 'FIND-TARGET-TF-TISSUE':respond_find_target_tfs_tissue,
                 'FIND-EVIDENCE-MIRNA-TARGET':respond_find_evidence_miRNA_target,
                 'FIND-EVIDENCE-MIRNA-EXP':respond_find_evidence_miRNA_exp,
                 'FIND-PATHWAY-DB-KEYWORD':respond_find_pathway_db_keyword,
                 'FIND-TISSUE':respond_find_tissue_gene,
                 'FIND-KINASE-REGULATION':respond_find_kinase_regulation,
                 'FIND-TF-MIRNA':respond_find_tf_miRNA,
                 'FIND-REGULATION':respond_find_regulation,
                 'FIND-EVIDENCE':respond_find_evidence,
                 'FIND-GENE-TISSUE':respond_find_gene_tissue,
                 'IS-GENE-TISSUE':respond_is_gene_tissue,
                 'FIND-KINASE-PATHWAY':respond_find_kinase_pathway,
                 'GO-ENRICHMENT':respond_go_enrichment, 'GO-ANNOTATION':respond_go_annotation,
                 'IS-MIRNA-DISEASE':respond_is_mirna_disease,
                 'FIND-MIRNA-DISEASE':respond_find_mirna_disease,
                 'FIND-DISEASE-MIRNA':respond_find_disease_mirna,
                 'MAKE-HEATMAP':respond_make_heatmap, 'PATHWAY-ENRICHMENT':respond_pathway_enrichment,
                 'DISEASE-ENRICHMENT': respond_disease_enrichment,
                 'MIRNA-DISEASE-ENRICHMENT': respond_mirna_disease_enrichment}
    
    def receive_request(self, msg, content):
        """If a "request" message is received, decode the task and
        the content and call the appropriate function to prepare the
        response. A reply message is then sent back.
        """
        task_str = content.head().upper()
        if task_str not in self.task_func:
            logger.info('In receive_request: TFTA received the unkonwn task: {}.'.format(task_str))
            reply_content = make_failure2('NO-CAPABILITY', task_str)
        else:
            reply_content = self.task_func[task_str](self, content)
            logger.info('In receive_request: TFTA received the task: {}.'.format(task_str))
        
        reply_msg = KQMLPerformative('reply')
        reply_msg.set('content', reply_content)
        
        logger.info('In receive_request: {} is sending message...'.format(task_str))
        
        self.reply(msg, reply_msg)
        
    def get_regulator_indra(self, target_names, stmt_types, keyword_name):
        """
        wrap message for multiple targets case
        target_names: list
        stmt_types: indra statement type
        """
        #lit_messages = ''
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
                return None,None,None,None
        tfs, genes, mirnas, others, stmt_f = self.tfta.find_regulators_indra(stmts_d)
        tf_json, gene_json, mirna_json, oth_json = [],[],[],[]
        if len(tfs):
            tf_json = self._get_genes_json(tfs)
            #lit_messages.set('tf-literature', tf_json)
        if len(genes):
            gene_json = self._get_genes_json(genes)
            #lit_messages.set('gene-literature', gene_json)
        if len(mirnas):
            mirna_agent = [Agent(mir, db_refs={'TYPE':'MIRNA'}) for mir in mirnas]
            mirna_json = self.make_cljson(mirna_agent)
            #lit_messages.set('miRNA-literature', mirna_json)
        if len(others):
            oth_agent = [Agent(oth) for oth in others]
            oth_json = self.make_cljson(oth_agent)
            #lit_messages.set('other-literature', oth_json)
            
        #provenance support
        self.send_background_support(stmt_f, 'what', target_str, keyword_name)
        
        return tf_json, gene_json, mirna_json, oth_json
        
    def get_regulator_indra_those(self, target_names, stmt_types, keyword_name, of_those=None, target_type=None):
        """
        wrap message for multiple targets case, use json format
        target_names: list
        stmt_types: indra statement type
        of_those: set
        """
        lit_json = ''
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
                return lit_json
        genes, stmt_f = self.tfta.find_regulator_indra_targets(stmts_d, of_those=of_those, target_type=target_type)
        
        if len(genes):
            lit_json = self._get_genes_json(genes)
            #lit_messages = KQMLList()
            #lit_messages.set('gene-literature', gene_json)
        #provenance support
        self.send_background_support(stmt_f, 'what', target_str, keyword_name)
        
        return lit_json
        
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
            target_list = [Agent(target, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for target in genes]
            lit_messages = self.make_cljson(target_list)
             
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
            tf_list = [Agent(tf, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for tf in tfs]
            lit_messages = self.make_cljson(tf_list)
            
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
        row_list = ['<th class = table-borders>miRNA</th><th class = table-borders>Target</th> \
                   <th class = table-borders>Experiment</th><th class = table-borders>Support Type</th> \
                   <th class = table-borders>PMID</th>']
        for mirna,target,expe,st,pd in zip(mirna_name, target_name, experiment, support_type, pmid):
            pd_list = pd.split(',')
            pd_str = ''
            for p in pd_list:
                pd_str += '<a href=' + publink + p + ' target="_blank">' + p + '</a>;'
            row_list.append('<td class = table-borders>%s</td><td class = table-borders>%s</td> \
                             <td class = table-borders>%s</td><td class = table-borders>%s</td> \
                             <td class = table-borders>%s</td>' % (mirna, target, expe, st, pd_str[:-1]))
        html_str += '\n'.join(['  <tr>%s</tr>\n' % row_str
                               for row_str in row_list])
        html_str += '</table> <hr>'
        content = KQMLList('add-provenance')
        content.sets('html', html_str)
        return self.tell(content)
        
    def send_background_support_mirna(self, mirna_name, target_name, experiment, support_type, pmid, strength=None, find_mirna=False, find_target=False):
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
                if strength:
                    nl = 'does ' + mirna_name + 'have ' + strength + 'evidence for targeting ' + target_name + '?'
                else:
                    nl = 'does ' + mirna_name + ' regulate ' + target_name + '?'
                exp_str = ';\n'.join(experiment[target_name])
                suptype_str = ';\n'.join(support_type[target_name])
                pmid_str = ','.join(pmid[target_name])
                self.send_table_to_provenance_mirna([mirna_name], [target_name], [exp_str], [suptype_str], [pmid_str], nl)
            elif find_mirna:
                if len(target_name) > 1:
                    target_str = ', '.join(target_name[:-1]) + ' and ' + target_name[-1]
                else:
                    target_str = target_name[0]
                if strength:
                    nl = 'what microRNAs have strong evidence for targeting ' + target_str + '?'
                else:
                    nl = 'what microRNAs regulate ' + target_str + '?'
                target_list = []
                exp_list = []
                sup_list = []
                pmid_list = []
                mirna_list = []
                for mir in mirna_name:
                    for target in target_name:
                        mirna_list.append(mir)
                        target_list.append(target)
                        exp_list.append(';\n'.join(experiment[(mir,target)]))
                        sup_list.append(';\n'.join(support_type[(mir,target)]))
                        pmid_list.append(','.join(pmid[(mir,target)]))
                self.send_table_to_provenance_mirna(mirna_list, target_list, exp_list, sup_list, pmid_list, nl)
            elif find_target:
                if len(mirna_name) > 1:
                    mirna_str = ', '.join(mirna_name[:-1]) + ' and ' + mirna_name[-1]
                else:
                    mirna_str = mirna_name[0]
                if strength:
                    nl = 'what genes are regulated by ' + mirna_str + ' with strong evidence?'
                else:
                    nl = 'what genes are regulated by ' + mirna_str + '?'
                mirna_list = []
                target_list = []
                exp_list = []
                sup_list = []
                pmid_list = []
                for target in target_name:
                    for mir in mirna_name:
                        mirna_list.append(mir)
                        target_list.append(target)
                        exp_list.append(';\n'.join(experiment[(mir,target)]))
                        sup_list.append(';\n'.join(support_type[(mir,target)]))
                        pmid_list.append(','.join(pmid[(mir,target)]))
                self.send_table_to_provenance_mirna(mirna_list, target_list, exp_list, sup_list, pmid_list, nl)
            else:
                return
        else:
            for_what = 'your query'
            cause_txt = 'miRNA-target db'
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
        
        #logger.info('html_str={}'.format(html_str))
        
        return self.tell(content)
        
    def send_background_support_db(self, tf_name, target_name, dbname, tissue=None, find_tf=False, find_target=False, target_type=None):
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
                    
                if target_type:
                    t_str = str(target_type) + 's'
                else:
                    t_str = 'genes'
                if tissue:
                    nl = 'what ' + t_str + ' are regulated by ' + tf_str + ' in ' + tissue + '?'
                else:
                    nl = 'what ' + t_str + ' are regulated by ' + tf_str + '?'
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
        kin_json = ''
        kinase_names = []
        if keyword_name == 'regulate':
            try:
                kinase_names = self.tfta.find_kinase_target(target_names)
            except KinaseNotFoundException:
                return kin_json
        else:
            try:
                kinase_names = self.tfta.find_kinase_target_keyword(target_names, keyword_name)
            except KinaseNotFoundException:
                return kin_json
        
        if of_those:
            kinase_names = list(set(kinase_names).intersection(set(of_those)))
        if len(kinase_names):
            kin_json = self._get_genes_json(kinase_names)
        return kin_json
        
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
        
    def wrap_message(self, descr, gene_names):
        #use json format
        gene_list = []
        for gene in gene_names:
            hgnc_id = hgnc_client.get_hgnc_id(gene)
            if hgnc_id:
                gene_list.append(Agent(gene, db_refs={'HGNC': hgnc_id, 'TYPE':'ONT::GENE-PROTEIN'}))
            else:
                gene_list.append(Agent(gene, db_refs={'TYPE':'ONT::GENE-PROTEIN'}))
        gene_json = self.make_cljson(gene_list)
        reply = KQMLList('SUCCESS')
        reply.set(descr, gene_json)
        return reply
        
    def _get_genes_json(self, gene_names):
        #use json format
        gene_list = []
        for gene in gene_names:
            hgnc_id = hgnc_client.get_hgnc_id(gene)
            if hgnc_id:
                gene_list.append(Agent(gene, db_refs={'HGNC': hgnc_id, 'TYPE':'ONT::GENE-PROTEIN'}))
            else:
                gene_list.append(Agent(gene, db_refs={'TYPE':'ONT::GENE-PROTEIN'}))
        gene_json = self.make_cljson(gene_list)
        return gene_json
    
    def get_target_type_set(self, content):
        target_type = _get_keyword_name(content, descr='target-type')
        target_type_set = set()
        if target_type:
            try:
                target_type_set = self.tfta.get_onto_set(target_type)
            except GONotFoundException:
                return target_type_set, target_type
        return target_type_set,target_type
    
    def wrap_family_message2(self, term_id, msg):
        """
        parameter
        -------------
        term_id: dict, key is trips id or family name, value is Agent
        msg: str
        """
        #use json format
        #only consider one family for clarification for now
        members = []
        if term_id:
            term = list(term_id.keys())[0]
            members = self.tfta.find_member(term_id[term])
        else:
            reply = make_failure(msg)
        if members:
            #add +type+ to db_refs{}
            for a in members:
                a.db_refs.update({'TYPE':'ONT::GENE-PROTEIN'})
            mbj = self.make_cljson(members)
            res_str = KQMLList('resolve')
            res_str.set('term', self.make_cljson(Agent(term, db_refs={'TYPE':'ONT::PROTEIN-FAMILY'})))
            res_str.set('as', mbj)
            reply = make_failure_clarification('FAMILY_NAME', res_str)
        else:
            reply = make_failure(msg)
        return reply
    
    def wrap_family_message(self, term_id, msg):
        if term_id:
            term = list(term_id.keys())[0]
            reply = self.make_resolve_family_failure(term_id[term])
        else:
            reply = make_failure(msg)
        return reply
        
        
    def wrap_family_message_pathway2(self, term_id, descr='pathways', msg="PATHWAY_NOT_FOUND"):
        #term_id = _get_term_id(target_arg)
        if term_id:
            term = list(term_id.keys())[0]
            members = self.tfta.find_member(term_id[term])
            if members:
                #id = list(members.keys())[0]
                #add +type+ to db_refs{}
                for a in members:
                    a.db_refs.update({'TYPE':'ONT::GENE-PROTEIN'})
                mbj = self.make_cljson(members)
                res_str = KQMLList('resolve')
                res_str.set('family', self.make_cljson(Agent(term, db_refs={'TYPE':'ONT::PROTEIN-FAMILY'})))
                res_str.set('as', mbj)
                reply = make_failure_clarification('FAMILY_NAME', res_str)
            else:
                reply = make_failure(msg)
        else:
            reply = KQMLList('SUCCESS')
            reply.set(descr, 'NIL')
        return reply
    
    def wrap_family_message_pathway(self, term_id, descr='pathways', msg="PATHWAY_NOT_FOUND"):
        if term_id:
            try:
                term = list(term_id.keys())[0]
                reply = self.make_resolve_family_failure(term_id[term])
            except Exception:
                reply = KQMLList('SUCCESS')
                reply.set(descr, 'NIL')
        else:
            reply = KQMLList('SUCCESS')
            reply.set(descr, 'NIL')
        return reply
    
    def _wrap_pathway_genelist_message(self, pathwayName, dblink, genelist, pathway_names=None,
                                   gene_descr='tfs', of_gene_names=None, of_those=None):
        """
        parameters
        -------------
        pathwayName: dict[pthid]
        dblink: dict[pthid]
        genelist: dict[pthid]
        pathway_names: list or None
        gene_descr: str
        of_gene_names: list
        of_those: list or None
        """
        limit = 30
        num = 1
        pathway_list_json = []
        keys = list(genelist.keys())
        if keys:
            if pathway_names:
                for key in keys:
                    if of_those:
                        if pathwayName[key].lower() in of_those:
                            if _filter_subword(pathwayName[key], pathway_names):
                                if of_gene_names:
                                    genes = set(of_gene_names) & set(genelist[key])
                                else:
                                    genes = genelist[key]
                                if genes:
                                    gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in genes]
                                    gene_json = self.make_cljson(gene_agent)
                                    mes = KQMLList()
                                    mes.sets('name', pathwayName[key])
                                    mes.sets('dblink', dblink[key])
                                    mes.set(gene_descr, gene_json)
                                    pathway_list_json.append(mes.to_string())
                                    #check the limit
                                    num += 1
                                    if num > limit:
                                        break
                    else:
                        if _filter_subword(pathwayName[key], pathway_names):
                            if of_gene_names:
                                genes = set(of_gene_names) & set(genelist[key])
                            else:
                                genes = genelist[key]
                            if genes:
                                #for debug
                                #logger.info('wrap_pathway_genelist_message: The {}th pathwayname={}.'.format(num+1, pathwayName[key]))
                                gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in genes]
                                gene_json = self.make_cljson(gene_agent)
                                mes = KQMLList()
                                mes.sets('name', pathwayName[key])
                                mes.sets('dblink', dblink[key])
                                mes.set(gene_descr, gene_json)
                                pathway_list_json.append(mes.to_string())
                                #check the limit
                                num += 1
                                if num > limit:
                                    break
            else:
                for key in keys:
                    if of_those:
                        if pathwayName[key].lower() in of_those:
                            if of_gene_names:
                                genes = set(of_gene_names) & set(genelist[key])
                            else:
                                genes = genelist[key]
                            if genes:
                                gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in genes]
                                gene_json = self.make_cljson(gene_agent)
                                mes = KQMLList()
                                mes.sets('name', pathwayName[key])
                                mes.sets('dblink', dblink[key])
                                mes.set(gene_descr, gene_json)
                                pathway_list_json.append(mes.to_string())
                                #check the limit
                                num += 1
                                if num > limit:
                                    break
                    else:
                        if of_gene_names:
                            genes = set(of_gene_names) & set(genelist[key])
                        else:
                            genes = genelist[key]
                        if genes:
                            #for debug
                            #logger.info('wrap_pathway_genelist_message: The {}th pathwayname={}.'.format(num+1, pathwayName[key]))
                            gene_agent = [Agent(g, db_refs={'TYPE':'ONT::GENE-PROTEIN'}) for g in genes]
                            gene_json = self.make_cljson(gene_agent)
                            mes = KQMLList()
                            mes.sets('name', pathwayName[key])
                            mes.sets('dblink', dblink[key])
                            mes.set(gene_descr, gene_json)
                            pathway_list_json.append(mes.to_string())
                            #check the limit
                            num += 1
                            if num > limit:
                                break
            if pathway_list_json:
                reply = KQMLList('SUCCESS')
                pathway_list_str = '(' + ' '.join(pathway_list_json) + ')'
                reply.set('pathways', pathway_list_str)
            else:
                reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
        else:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
        return reply
    
    def _get_mirna_clarification(self, miRNA_mis):
        """
        miRNA_mis: dict, key is mirna name, value is trips term id or mirna name
        """
        #only consider the single miRNA case
        mirna = list(miRNA_mis.keys())[0]
        clari_mirna = self.tfta.get_similar_miRNAs(mirna)
        if clari_mirna:
            mir_agent = [Agent(mir, db_refs={'TYPE':'MIRNA'}) for mir in clari_mirna]
            mir_json = self.make_cljson(mir_agent)
            res = KQMLList('resolve')
            res.set('term', miRNA_mis[mirna])
            res.set('as', mir_json)
            reply = make_failure_clarification('MIRNA_NOT_FOUND', res)
        else:
            reply = make_failure('NO_SIMILAR_MIRNA')
        return reply
    
    def _get_targets(self, content, descr='target'):
        #parse json message format
        proteins = []
        family = dict()
        
        target_arg = content.get(descr)
        if not target_arg:
            return None,None
        try:
            agents = self.get_agent(target_arg)
        except Exception:
            return None,None
        if isinstance(agents, list):
            proteins = [a.name.upper() for a in agents if a is not None and ('UP' in a.db_refs or 'HGNC' in a.db_refs or len(a.db_refs)==0)]
            #family = {a.db_refs['TRIPS']:a.name for a in agents if a is not None and 'FPLX' in a.db_refs and a.name not in proteins}
            #consider +trips+ as an optional id
            for a in agents:
                if a is not None and 'FPLX' in a.db_refs and a.name.upper() not in proteins:
                    if 'TRIPS' in a.db_refs:
                        family[a.db_refs['TRIPS']] = a
                    else:
                        family[a.name.upper()] = a
        elif isinstance(agents, Agent):
            if 'UP' in agents.db_refs or 'HGNC' in agents.db_refs or len(agents.db_refs)==0:
                proteins = [agents.name.upper()]
            if not proteins and 'FPLX' in agents.db_refs:
                #family = {agents.db_refs['TRIPS']:agents.name}
                if 'TRIPS' in agents.db_refs:
                    family[agents.db_refs['TRIPS']] = agents
                else:
                    family[agents.name.upper()] = agents
        return proteins,family
        
    def _get_targets2(self, content, descr='target'):
        #parse json message format
        #expand members for family name
        proteins = []
        family = []
        fmembers = dict()
        
        target_arg = content.get(descr)
        if not target_arg:
            return None,None
        try:
            agents = self.get_agent(target_arg)
        except Exception:
            return None,None
        if isinstance(agents, list):
            proteins = [a.name.upper() for a in agents if a is not None and ('UP' in a.db_refs or 'HGNC' in a.db_refs or len(a.db_refs)==0)]
            #expand family name to members
            for a in agents:
                if a is not None and 'FPLX' in a.db_refs and a.name.upper() not in proteins:
                    family.append(a)
            fmembers = self.tfta.find_members(family)
        elif isinstance(agents, Agent):
            if 'UP' in agents.db_refs or 'HGNC' in agents.db_refs or len(agents.db_refs)==0:
                proteins = [agents.name.upper()]
            if not proteins and 'FPLX' in agents.db_refs:
                fmembers[agents.name.upper()] = self.tfta.find_member(agents)
        return proteins,fmembers
        
    def _get_of_those_pathway(self, content, descr='of-those'):
        #parse pathway names from JSON format
        pathways = []
        path = content.get(descr)
        try:
            agents = self.get_agent(path)
            if isinstance(agents, list):
                pathways = [a.name.lower() for a in agents if a is not None]
            elif isinstance(agents, Agent):
                pathways = [agents.name]
            return pathways
        except Exception:
            return None
        
    def _combine_json_list(self, descr_list, json_list):
        if any(json_list):
            res = KQMLList()
            for i in range(len(descr_list)):
                if json_list[i]:
                    res.set(descr_list[i], json_list[i])
            return res
        else:
            return None
    
    def _get_mirnas(self, content, descr='miRNA'):
        #JSON format
        #consider +trips+ if it exists
        mirna = {}
        try:
            mir_arg = content.get(descr)
        except Exception:
            return []
        try:
            agents = self.get_agent(mir_arg)
        except Exception:
            return []
        if isinstance(agents, list):
            for a in agents:
                if a is not None:
                    mname = _get_mirna_name(a.name.upper())
                    if 'TRIPS' in a.db_refs:
                        mirna[mname] = a.db_refs['TRIPS']
                    else:
                        mirna[mname] = mname
        elif isinstance(agents, Agent):
            if agents is not None:
                mname = _get_mirna_name(agents.name.upper())
                if 'TRIPS' in agents.db_refs:
                    mirna[mname] = agents.db_refs['TRIPS']
                else:
                    mirna[mname] = mname
        return mirna
        
    def _get_pathway_name(self, content, descr='pathway'):
        #JSON format
        pathways = set()
        keyword = ['signaling pathway', 'pathway']
        try:
            pth_arg = content.get(descr)
        except Exception:
            return []
        try:
            agents = self.get_agent(pth_arg)
        except Exception:
            return []
        if isinstance(agents, list):
            for a in agents:
                if a is not None:
                    p = _filter_pathway_name(a.name, keyword)
                    if p:
                        pathways.add(p.replace('-', ' '))
                    if a.db_refs.get('TEXT'):
                        p = _filter_pathway_name(a.db_refs.get('TEXT'), keyword)
                        if p:
                            pathways.add(p)
        elif isinstance(agents, Agent):
            if agents is not None:
                p = _filter_pathway_name(agents.name, keyword)
                if p:
                    pathways.add(p.replace('-', ' '))
                if agents.db_refs.get('TEXT'):
                    p = _filter_pathway_name(agents.db_refs['TEXT'], keyword)
                    if p:
                        pathways.add(p)
        return list(pathways)
        
    def _to_json_mirna(self, mirnas):
        mir_agent = [Agent(mir, db_refs={'type':'MIRNA'}) for mir in mirnas]
        mir_json = self.make_cljson(mir_agent)
        return mir_json
        
#------------------------------------------------------------------------#######

def _filter_pathway_name(path_name, keyword):
    p = ''
    for k in keyword:
        if path_name:
            p = path_name.lower().replace('signaling pathway', '').strip()
            p = p.replace('pathway', '').strip()
    return p

def make_failure(reason):
    msg = KQMLList('FAILURE')
    msg.set('reason', reason)
    return msg
    
def make_failure2(reason, task_str):
    msg1 = KQMLList(reason)
    msg1.set('name', task_str)
    msg = KQMLList('FAILURE')
    msg.set('reason', msg1)
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
    if descr:
        for d in descr:
            if d[-len(word):] == word:
                if len(d[:-len(word)-1]):
                    ds.append(d[:-len(word)-1])
            else:
                ds.append(d)
    return ds

def _get_mirna_name(str1):
    #handle two forms of input, like MIR-PUNC-MINUS-20-B-PUNC-MINUS-5-P and MIR-20-B-5-P
    if 'PUNC-MINUS' in str1:
        str2 = str1.replace('-PUNC-MINUS-','_')
        str2 = str2.replace('-','')
        str2 = str2.replace('_', '-')
        return str2.upper()
    else:
        plist = re.findall('([0-9]+-[a-zA-Z])', str1)
        s = str1
        for p in plist:
            p1 = p.replace('-','')
            s = s.replace(p, p1)
        return s.upper()

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
    
def _make_evidence_json(evids, limit = 10):
    ind = 0
    ev_json = []
    for ev in evids:
        if all(ev) and ind < limit:
            ev_json.append(Agent(ev[2], db_refs={'source_api':ev[0], 'pmid':ev[1]}))
            ind += 1
        else:
            break
    return ev_json

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
    
def _wrap_dbname_message(descr, db_names):
    db_str = ''
    for db in db_names:
        try:
            pmid = dbname_pmid_map[db]
        except KeyError as e:
            pmid = ''
        db_str += '(:name %s :pmid %s) ' % (db, pmid)
    db_str = descr + '( ' + db_str + ') '
    return db_str
    
def _combine_messages(mess_list):
    messages = ''
    for mess in mess_list:
        if mess:
            messages += mess
    return messages

def _wrap_pathway_message(pathwayName, dblink, keyword=None, of_those=None):
    """
    pathwayName: list
    dblink: list
    keyword: list or None
    of_those: list or None
    """
    pathway_list_str = ''
    #limit the number of pathways to return
    limit = 30
    num = 1
    if keyword:
        if type(keyword).__name__ == 'str':
            keyword = [keyword]
        if of_those:
            for pn, dbl in zip(pathwayName, dblink):
                if pn.lower() in of_those:
                    if _filter_subword(pn, keyword):
                        pnslash = '"' + pn +'"'
                        dbl = '"' + dbl +'"'
                        pathway_list_str += '(:name %s :dblink %s) ' % (pnslash, dbl)
                        num += 1
                        if num > limit:
                            break
        else:
            for pn, dbl in zip(pathwayName, dblink):
                if _filter_subword(pn, keyword):
                    pnslash = '"' + pn +'"'
                    dbl = '"' + dbl +'"'
                    pathway_list_str += '(:name %s :dblink %s) ' % (pnslash, dbl)
                    num += 1
                    if num > limit:
                        break
        if pathway_list_str:
            pathway_list_str = '(' + pathway_list_str + ')'
            reply = KQMLList('SUCCESS')
            reply.set('pathways', pathway_list_str)
        else:
            reply = KQMLList.from_string('(SUCCESS :pathways NIL)')
    else:
        if of_those:
            for pn, dbl in zip(pathwayName, dblink):
                if pn.lower() in of_those:
                    pnslash = '"' + pn + '"'
                    dbl = '"' + dbl + '"'
                    pathway_list_str += '(:name %s :dblink %s) ' % (pnslash, dbl)
                    num += 1
                    if num > limit:
                        break
        else:
            for pn, dbl in zip(pathwayName, dblink):
                pnslash = '"' + pn + '"'
                dbl = '"' + dbl + '"'
                pathway_list_str += '(:name %s :dblink %s) ' % (pnslash, dbl)
                num += 1
                if num > limit:
                    break
        pathway_list_str = '(' + pathway_list_str + ')'
        reply = KQMLList('SUCCESS')
        reply.set('pathways', pathway_list_str)
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
