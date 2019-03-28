from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor

#####################################
#Testing the following TFTA functions
#FIND-KINASE-PATHWAY
#FIND-EVIDENCE-MIRNA-TARGET
#TEST-FUNCTION
#UNIT TEST
#
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



##############################################################################
#FIND-KINASE-PATHWAY
#What kinases are in the MAPK signaling pathway?
class TestFindKinasePathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 21, output
        
#What kinases are in the immune system pathway?
class TestFindKinasePathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output
        
#subsequent query
#Which of these kinases are also in the MAPK signaling pathway?
class TestFindKinasePathway3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        of_those = 'AKT1, FGFR1,CAMK2B,MAP2K1,RAF1,MAPK1,JAK1,MAP3K1'
        #print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 18, output
        
#Which of these kinases are also in the MAPK signaling pathway?
class TestFindKinasePathway31(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway31, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        of_those = ekb_kstring_from_text('AKT1, FGFR1,CAMK2B,MAP2K1,RAF1,MAPK1,JAK1,MAP3K1')
        #print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 18, output
        
#What of these kinases are in the immune system pathway?
class TestFindKinasePathway4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        of_those = ekb_kstring_from_text('AKT1, FGFR1,CAMK2B,MAP2K1,RAF1,MAPK1,JAK1,MAP3K1')
        #print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output
        
#What of these kinases are in the immune system pathway?
class TestFindKinasePathway41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway41, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        of_those = 'AKT1, FGFR1,CAMK2B,MAP2K1,RAF1,MAPK1,JAK1,MAP3K1'
        #print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output

#########################################################################
#TEST FIND-EVIDENCE-MIRNA-TARGET
#show me evidence that miR-148a-3p targets DNMT1?
class TestFindEvidenceMirnaTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidenceMirnaTarget1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-148a-3p')
        target = ekb_kstring_from_text('DNMT1')
        content = KQMLList('FIND-EVIDENCE-MIRNA-TARGET')
        content.set('miRNA', mirna)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", len(output.get('evidence')))
        assert len(output.get('evidence')) == 7, output

#clarification
#show me evidence that miR-148a targets DNMT1?
class TestFindEvidenceMirnaTarget2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidenceMirnaTarget2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-148a')
        target = ekb_kstring_from_text('DNMT1')
        content = KQMLList('FIND-EVIDENCE-MIRNA-TARGET')
        content.set('miRNA', mirna)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        assert len(output.get('clarification').get('as')) == 2, output
        
#show me evidence that miR-148 targets DNMT1?
class TestFindEvidenceMirnaTarget3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidenceMirnaTarget3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-148')
        target = ekb_kstring_from_text('DNMT1')
        content = KQMLList('FIND-EVIDENCE-MIRNA-TARGET')
        content.set('miRNA', mirna)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_SIMILAR_MIRNA', output


############################################
#TEST-FUNCTION
class TestGetFamilyName(_IntegrationTest):
    def __init__(self, *args):
        super(TestGetFamilyName, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        family_arg = ekb_from_text('SMURF')
        print('family_arg=', family_arg)
        family_name = get_family_name(family_arg)
        content = KQMLList('TEST-FUNCTION')
        content.set('family', KQMLString(family_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("len(output.get('reason'))=", len(output.get('reason')))
        #It seems KQML calculate the number by space split
        assert len(output.get('reason')) == 4, output
    
#########################################
##UNIT TEST##############################
def test_find_tf_indra():
    tfta = TFTA()
    stmts,ind = tfta.find_statement_indraDB(obj='MYC', stmt_types=['IncreaseAmount'])
    print('len(stmts)=', len(stmts))
    print('ind=', ind)
    tfs, genes, mirnas, other = tfta.find_regulator_indra(stmts)
    print('len(tfs)=', (len(tfs)))
    print('len(mirnas)=', len(mirnas))
    print('len(genes)=', len(genes))
    print('len(other)=', len(other))
    assert(len(tfs)>0)
    assert(len(mirnas)>0)
    assert(len(genes)>0)
    assert(len(other)>0)
    
def test_map_exclusive_tissue_gene():
    tfta = TFTA()
    gene_exp_exclusive = tfta.map_exclusive_tissue_gene()
    print('len(gene_exp_exclusive)=', len(gene_exp_exclusive))
    assert(len(gene_exp_exclusive) == 26)


if __name__ == '__main__':
    TestIsRegulation1().run_test()