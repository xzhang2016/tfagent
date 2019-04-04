from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor

#####################################
#Testing the following TFTA functions
#FIND-TARGET-MIRNA
#FIND-GENE-COUNT-MIRNA
#FIND-GENE-PATHWAY
#FIND-PATHWAY
#find-common-pathway-genes
#IS-PATHWAY-GENE
#FIND-PATHWAY-DB-KEYWORD
#FIND-MIRNA
#FIND-MIRNA-COUNT-GENE
#FIND-TISSUE
######################################

#define some functions
def get_gene_symbol(target_arg):
        agent = []
        targets = []
        ont1 = ['ONT::PROTEIN', 'ONT::GENE-PROTEIN', 'ONT::GENE']
        tp = TripsProcessor(target_arg)
        for term in tp.tree.findall('TERM'):
            if term.find('type').text in ont1:
                term_id = term.attrib['id']
                agent.append(tp._get_agent_by_id(term_id, None))
        for ag in agent:
            targets.append(ag.name)
        print('targets=' + ','.join(targets))

def get_family_name(target_arg):
        agent = []
        targets = []
        ont1 = ['ONT::PROTEIN-FAMILY', 'ONT::GENE-FAMILY']
        tp = TripsProcessor(target_arg)
        for term in tp.tree.findall('TERM'):
            if term.find('type').text in ont1:
                term_id = term.attrib['id']
                agent.append(tp._get_agent_by_id(term_id, None))
        for ag in agent:
            targets.append(ag.name)
        print('family_name = ' + ','.join(targets))

def get_family_id(target_arg):
        term_id = []
        ont1 = ['ONT::PROTEIN-FAMILY', 'ONT::GENE-FAMILY']
        tp = TripsProcessor(target_arg)
        for term in tp.tree.findall('TERM'):
            if term.find('type').text in ont1:
                term_id.append(term.attrib['id'])
        print('term_id = ' + ','.join(term_id))
        
        
#####################################################################
#TEST FIND-TARGET-MIRNA
#What genes does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-20b-5p')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 917, output
        
#What genes are regulated by miR-297? (subtask: find-target-mirna)
class TestFindTargetMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-297')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 190, output
        
#What genes are regulated by miR-20b-5p and MIR-29B-1-5P? (subtask: find-target-mirna)
class TestFindTargetMirna21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #mirna = ekb_kstring_from_text('miR-20b-5p, miR-29B-1-5P')
        mirna = ekb_kstring_from_text('miR-20b-5p and MIR-29B-1-5P')
        print("mirna=", mirna)
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 12, output
        
        
#What are the genes that have strong evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('mir-122-5p')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'strong')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 69, output
        
#What are the genes that have weak evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('mir-122-5p')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'weak')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 556, output

#clarification test        
#What are the genes that have weak evidence of being regulated by mir-128? (subtask: find-target-mirna)
class TestFindTargetMirna5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('mir-128')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'weak')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("output.get('reason') = ", output.get('reason'))
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        assert len(output.get('clarification').get('as')) == 3, output
        
#What genes are regulated by miR-200C? (subtask: find-target-mirna)
class TestFindTargetMirna6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-200C')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("output.get('reason') = ", output.get('targets'))
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        assert len(output.get('clarification').get('as')) == 2, output
        
#subsequent query test
#Which of those genes does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-20b-5p')
        of_those = 'STAT3,SRF,HRAS,KRAS,ELK1,JAK1,JAK2,FOS'
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 2, output
        
#Which of those genes does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna71(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna71, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-20b-5p')
        of_those = ekb_kstring_from_text('STAT3,SRF,HRAS,KRAS,ELK1,JAK1,JAK2,FOS')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 2, output
        
#What are the genes that have weak evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna8(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('mir-122-5p')
        of_those = ekb_kstring_from_text('STAT3,SRF,CDK4, CDK19,CSRP1,DZIP1L,HSPA4L')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'weak')
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 5, output
        
#What are the genes that have weak evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna81(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna81, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('mir-122-5p')
        of_those = 'STAT3,SRF,CDK4, CDK19,CSRP1,DZIP1L,HSPA4L'
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'weak')
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 5, output
        
#test target-type
#What kinases does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna9(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-20b-5p')
        target_type = 'kinase'
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('target-type', target_type)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 40, output
        
#What transcription facotrs does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna91(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna91, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-20b-5p')
        target_type = 'transcription factor'
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('target-type', target_type)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 130, output
        
#####################################################################################
#What genes are most frequently regulated by miR-335-5p, miR-155-5p, miR-145-5p, and miR-20a-5p?
#(subtask: FIND-GENE-COUNT-MIRNA)
class TestFindGeneCountMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneCountMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-335-5p, miR-155-5p, miR-145-5p, miR-20a-5p')
        count = 'count'
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 192, output
        
#What genes are most frequently regulated by miR-335-5p, miR-155-5p and miR-145-5p?
#(subtask: FIND-GENE-COUNT-MIRNA)
class TestFindGeneCountMirna11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneCountMirna11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-335-5p, miR-155-5p, miR-145-5p')
        count = 'count'
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 33, output
        
#clarification
#What genes are most frequently regulated by miR-128, miR-200c, and miR-20a-5p?
class TestFindGeneCountMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneCountMirna2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-122, miR-200c, and miR-20a-5p')
        print('mirna = ', str(mirna))
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        print("len(output.get('clarification').get('as'))=",len(output.get('clarification').get('as')))
        assert len(output.get('clarification').get('as')) == 2, output
        
#subsequent query
#Which of those genes are most frequently regulated by miR-335-5p, miR-155-5p, miR-145-5p, and miR-20a-5p?
#(subtask: FIND-GENE-COUNT-MIRNA)
class TestFindGeneCountMirna3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneCountMirna3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-335-5p, miR-155-5p, miR-145-5p, miR-20a-5p')
        of_those = 'stat3,srf,hras,CDK19, HSPA4L,FOXRED2,ZBTB25,cd28,sp4,TNKS2'
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 6, output

#Which of those genes are most frequently regulated by miR-335-5p, miR-155-5p, miR-145-5p, and miR-20a-5p?
#(subtask: FIND-GENE-COUNT-MIRNA)
class TestFindGeneCountMirna31(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneCountMirna31, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-335-5p, miR-155-5p, miR-145-5p, miR-20a-5p')
        of_those = ekb_kstring_from_text('stat3,srf,hras,CDK19, HSPA4L,FOXRED2,ZBTB25,cd28,sp4,TNKS2')
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 6, output


#####################################################################################
#TEST FIND-GENE-PATHWAY
#What genes are in the MAPK signaling pathway? (subtask: find-gene-pathway)
class TestFindGenePathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 21, output
        
#What genes are involved in the il-12 pathway?
class TestFindGenePathway11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('il-12')
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        content.set('subtype', 'gene')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 3, output
        
#What genes are in the immune system pathway?
class TestFindGenePathway12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        content.set('subtype', 'gene')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output

#subsequent query
#Which of these genes are in the immune system pathway?
class TestFindGenePathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        of_genes = 'STAT3,SRF,KRAS,HRAS,FZD8,JAK1,JAK2,FOS'
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_genes)
        print('content=',str(content))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output
        
#Which of STAT3,SRF,KRAS,HRAS,FZD8,JAK1,JAK2,FOS are in the immune system pathway?
class TestFindGenePathway21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        of_genes = ekb_kstring_from_text('STAT3,SRF,KRAS,HRAS,FZD8,JAK1,JAK2,FOS')
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_genes)
        print('content=',str(content))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output
        
#Which of STAT3,SRF, FZD8,JAK1,JAK2,FOS are in the immune system pathway?
class TestFindGenePathway22(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        of_genes = ekb_kstring_from_text('STAT3,SRF, FZD8,JAK1,JAK2,FOS')
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_genes)
        #print('content=',str(content))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output
###################################################################################
#TEST FIND-PATHWAY
#What pathways involve SRF? (subtask: find-pathway-gene)
class TestFindPathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SRF')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 27, output
        
#What pathways involve kras and elk1? (subtask: find-pathway-gene)
class TestFindPathway11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('KRAS, ELK1')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 23, output
        
##What pathways involve MEK? (subtask: find-pathway-gene)
class TestFindPathway12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('MEK')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output

##What pathways involve AKT? (subtask: find-pathway-gene)        
class TestFindPathway13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('AKT')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('pathways') == 'NIL', output
        
#Which Reactome pathways utilize SRF? (subtask: find-pathway-db-gene)
class TestFindPathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SRF')
        db = 'Reactome'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#Which KEGG pathways utilize SRF? (subtask: find-pathway-db-gene)
class TestFindPathway21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SRF')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
##Which KEGG pathways utilize MEK? (subtask: find-pathway-db-gene)
class TestFindPathway22(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('MEK')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#Which KEGG pathways utilize AKT? (subtask: find-pathway-db-gene)
class TestFindPathway23(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway23, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('AKT')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('pathways') == 'NIL', output
        
#What pathways involve calcium? (subtask: find-pathway-keyword)
class TestFindPathwayKeyword3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        content = KQMLList('FIND-PATHWAY-KEYWORD')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 9, output
        
#What pathways involve calcium? (subtask: find-pathway-keyword covered by find-pathway)
class TestFindPathway31(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway31, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 9, output

#What pathways involve immune system?        
class TestFindPathwayKeyword32(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword32, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'immune-system'
        content = KQMLList('FIND-PATHWAY-KEYWORD')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 4, output

#What pathways involve immune system? 
class TestFindPathway33(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway33, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'immune-system'
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 4, output
        
#What signaling pathways are shared by STAT3 and SRF? (subtask: find-common-pathway-genes)
class TestFindPathway4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3, SRF')
        count = 'count'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 6, output
        
#What immune pathways involve kras and elk1? (subtask: find-pathway-gene-keyword)
class TestFindPathway5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('kras, elk1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#What immune pathways involve tap1 and jak1? (subtask: find-pathway-gene-keyword)
class TestFindPathway51(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway51, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('tap1, jak1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 1, output
        
##What immune pathways involve MEK and jak1? (subtask: find-pathway-gene-keyword)
#here it returns results by just ignoring MEK gene
class TestFindPathway52(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway52, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('MEK, jak1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 3, output
        
#What immune pathways involve AKT and jak1? (subtask: find-pathway-gene-keyword)
class TestFindPathway53(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway53, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('AKT, jak1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('pathways') == 'NIL', output
        
##What immune pathways involve MEK? (subtask: find-pathway-gene-keyword)
class TestFindPathway54(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway54, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('MEK')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#Which immune system pathways involve stat3?
#(subtask: find-pathway-gene-keyword)
class TestFindPathway55(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway55, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        keyword = 'immune-system'
        #count = 'count'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output

#test the limit
class TestFindPathway6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'signaling-pathway'
        #count = 'count'
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_PATHWAY_NAME', output

#Which pathways utilize Stat3? (subtask: find-pathway-gene)
class TestFindPathway21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways')) =", len(output.get('pathways')))
        assert len(output.get('pathways')) == 50, output
#################################################################################
#find-common-pathway-genes
#Which pathways are shared by STAT3, SOCS3, IFNG, FOXO3, and CREB5 genes? 
#(subtask: find-common-pathway-genes)
class TestFindCommonPathwayGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, IFNG, FOXO3, CREB5')
        #count = 'count'
        content = KQMLList('find-common-pathway-genes')
        content.set('target', target)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 8, output
        
#What signaling pathways are shared by STAT3, AKT and SRF? (subtask: find-common-pathway-genes)
class TestFindCommonPathwayGene11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, AKT, SRF')
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 6, output
        
#What signaling pathways are shared by STAT3, MEK and SRF? (subtask: find-common-pathway-genes)
#MEK is ignored
class TestFindCommonPathwayGene12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, MEK, SRF')
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 6, output
        
#which of those are in the immune pathways? 
class TestFindCommonPathwayGene13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = 'STAT3, MEK1, SRF, HRAS, KRAS, JAK2, JAK1'
        keyword = 'immune'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways')) = ", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output
        
#which of those have common pathways? 
class TestFindCommonPathwayGene14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene14, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = 'STAT3, MEK1, SRF, HRAS, KRAS, JAK2, JAK1'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways')) = ", len(output.get('pathways')))
        assert len(output.get('pathways')) == 50, output
        
#Which immune pathways are shared by STAT3, SOCS3, and CREB5 genes?
#(subtask: find-common-pathway-genes-keyword)
class TestFindCommonPathwayGene15(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene15, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, CREB5')
        keyword = 'immune'
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
##Which immune pathways are shared by AKT, STAT3, SOCS3, and CREB5 genes?
#(subtask: find-common-pathway-genes-keyword)
class TestFindCommonPathwayGene16(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene16, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('AKT, STAT3, SOCS3, CREB5')
        keyword = 'immune'
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#Which immune pathways are shared by MEK, STAT3, SOCS3, and CREB5 genes?
#(subtask: find-common-pathway-genes-keyword)
class TestFindCommonPathwayGene17(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene17, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK, STAT3, SOCS3, CREB5')
        keyword = 'immune'
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#What KEGG pathways involve ERBB2, JUN, and MAPK8?
#sub-task:find-common-pathway-genes-db
class TestFindCommonPathwayGene2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('ERBB2, JUN, MAPK8')
        db = 'KEGG'
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('database', db)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#####################################################################################
#TEST IS-PATHWAY-GENE
#Does the mTor pathway utilize SGK1? (subtask: is-pathway-gene)
class TestIsPathwayGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SGK1')
        pathway = ekb_kstring_from_text('mTor pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 1, output
        
##Does the mTor pathway utilize MEK? (subtask: is-pathway-gene)
class TestIsPathwayGene11(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('MEK')
        pathway = ekb_kstring_from_text('mTor pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#Does the mTor pathway utilize AKT? (subtask: is-pathway-gene)
class TestIsPathwayGene12(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_from_text('AKT')
        get_gene_symbol(gene)
        get_family_name(gene)
        pathway = ekb_kstring_from_text('mTor pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', KQMLString(gene))
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert output.get('pathways') == 'NIL', output

#################################################################################
#FIND-PATHWAY-DB-KEYWORD
#What KEGG pathways involve immune signaling? (subtask: find-pathway-db-keyword)
class TestFindPathwayDbKeyword1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayDbKeyword1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        database = 'KEGG'
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY-DB-KEYWORD')
        content.set('database', database)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 1, output
        
#What KEGG pathways involve immune system? (subtask: find-pathway-db-keyword)
class TestFindPathwayDbKeyword11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayDbKeyword11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        database = 'KEGG'
        keyword = 'immune system'
        content = KQMLList('FIND-PATHWAY-DB-KEYWORD')
        content.set('database', database)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('pathways') == 'NIL', output
        
#What reactome pathways involve immune system? (subtask: find-pathway-db-keyword)
class TestFindPathwayDbKeyword12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayDbKeyword12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        database = 'reactome'
        keyword = 'immune system'
        content = KQMLList('FIND-PATHWAY-DB-KEYWORD')
        content.set('database', database)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output

###########################################################################
#TEST FIND-MIRNA
#What microRNAs target STAT3? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3')
        content = KQMLList('FIND-MIRNA')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('miRNAs')) == 80, output
        
##What microRNAs target MEK? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        content = KQMLList('FIND-MIRNA')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#What microRNAs target AKT? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('AKT')
        content = KQMLList('FIND-MIRNA')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('miRNAs') == 'NIL', output

######################################################################
#FIND-MIRNA-COUNT-GENE
#What miRNAs most frequently regulate EGFR, SRF, STAT3, JAK2, and SMAD3?
#(subtask: FIND-MIRNA-COUNT-GENE)
class TestFindMirnaCountGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirnaCountGene1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('EGFR, SRF, STAT3, JAK2, SMAD3')
        count = 'count'
        content = KQMLList('FIND-MIRNA-COUNT-GENE')
        content.set('target', target)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('miRNAs')) == 23, output

############################################################################
#find-tissue
#What tissues is STAT3 expressed in?
class TestFindTissue1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tissue'))=", str(len(output.get('tissue'))))
        assert len(output.get('tissue')) == 8, output

#What tissues is MEK expressed in?
class TestFindTissue11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('MEK')
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#What tissues is AKT expressed in?
class TestFindTissue12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('AKT')
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tissue') == 'NIL', output
        
#what tissues can I ask 
class TestFindTissue13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = ekb_kstring_from_text('AKT')
        content = KQMLList('FIND-TISSUE')
        #content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tissue'))=", len(output.get('tissue')))
        assert len(output.get('tissue')) == 30, output
        
#What tissues is frizzled8 expressed in?
class TestFindTissue14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue14, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('frizzled8')
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tissue'))=", str(len(output.get('tissue'))))
        assert len(output.get('tissue')) == 7, output
