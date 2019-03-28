from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor

#####################################
#Testing the following TFTA functions
#IS-GENE-ONTO
#FIND-GENE-ONTO
#FIND-KINASE-REGULATION
#FIND-TF-MIRNA
#FIND-REGULATION
#FIND-EVIDENCE
#FIND-GENE-TISSUE
#IS-TISSUE-GENE
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


####################################################################################
#IS-GENE-ONTO
#Is stat3 a kinase?    
class TestIsGeneOnto1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        keyword = 'kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output

#Is stat3 a transcription factor?        
class TestIsGeneOnto2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        keyword = 'transcription factor'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
#Is stat3 a protein kinase?
class TestIsGeneOnto3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        keyword = 'protein kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
#Is jak1 a protein kinase?
class TestIsGeneOnto4(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('JAK1')
        keyword = 'protein kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

#Is PBRM1 a transcription factor?
class TestIsGeneOnto5(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('PBRM1')
        keyword = 'transcription factor'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

#TEST FAMILY NAME
#Is SMURF a transcription factor?
class TestIsGeneOnto6(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SMURF')
        keyword = 'transcription factor'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_GENE_NAME', output.get('reason')
        print("len(output.get('clarification').get('agents'))=",len(output.get('clarification').get('agents')))
        assert len(output.get('clarification').get('agents')) == 1, output
        
#################################################################################
#TEST FIND-GENE-ONTO
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are protein kinases?
class TestFindGeneOnto1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = "STAT3, JAK1, JAK2, ELK1, FOS"
        keyword = 'protein kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are protein kinases?
class TestFindGeneOnto11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text("STAT3, JAK1, JAK2, ELK1, FOS")
        keyword = 'protein kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are kinases?
class TestFindGeneOnto2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = "STAT3, JAK1, JAK2, ELK1, FOS"
        keyword = 'kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are kinases?
class TestFindGeneOnto21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text("STAT3, JAK1, JAK2, ELK1, FOS")
        keyword = 'kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are histone demethylase?
class TestFindGeneOnto3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = "STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B"
        keyword = 'histone demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output
        
#Among STAT3, JAK1, JAK2, ELK1, FOS and KDM4B, which are histone demethylase?
class TestFindGeneOnto31(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto31, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text("STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B")
        keyword = 'histone demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are demethylase?
#test use gene list in class variable
class TestFindGeneOnto4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = "STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B"
        keyword = 'demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are demethylase?
class TestFindGeneOnto41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto41, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text("STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B")
        keyword = 'demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output
        
#Among PBRM1, SMAD2, TBL1XR1, AKT1, CDK19, CDK8, CDK9, DDR1, GSK3A, GSK3B, MET,TRIM28,COL2A1,
# JAK1, PRMT1, RB1, SMURF2, TRAF4, and USP15, which are transcription factors?
class TestFindGeneOnto5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = "PBRM1, SMAD2, TBL1XR1, AKT1, CDK19, CDK8, CDK9, DDR1, \
                GSK3A, GSK3B, MET,TRIM28,COL2A1,JAK1, PRMT1, RB1, SMURF2, TRAF4, USP15"
        keyword = 'transcription factor'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print("len(output.get('genes'))=" + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 3, output
        
#complex query: find-target and find-gene-onto
#What genes regulated by FOS are kinases?
class TestFindGeneOnto6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = ekb_kstring_from_text('cfos')
        keyword = 'kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('regulator', regulator)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 5, output
        
#What genes regulated by FOS are kinases?
class TestFindGeneOnto61(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto61, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = 'FOS'
        keyword = 'kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('regulator', regulator)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 5, output

#######################################################################################
#FIND-KINASE-REGULATION
#Which kinases regulate the cfos gene?
class TestFindKinaseReg1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cfos')
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert len(output.get('kinase')) == 5, output
        
#test gene family
#Which kinases regulate the MEK gene?
class TestFindKinaseReg11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#Which kinases regulate the AKT gene?
class TestFindKinaseReg12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg12, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('akt')
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert output.get('kinase') == 'NIL', output
    
#Which kinases negatively regulate the cfos gene?
class TestFindKinaseReg2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cfos')
        keyword = 'decrease'
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert len(output.get('kinase')) == 3, output
        
#Which kinases positively regulate the cfos gene?
class TestFindKinaseReg21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg21, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cfos')
        keyword = 'increase'
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert len(output.get('kinase')) == 2, output
        
#test gene family
#Which kinases positively regulate the MEK gene?
class TestFindKinaseReg22(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg22, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        keyword = 'increase'
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output

#Which kinases positively regulate the AKT gene?
class TestFindKinaseReg23(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg23, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('AKT')
        keyword = 'increase'
        print(target)
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert output.get('kinase') == 'NIL', output

##############################################################################
#TEST FIND-TF-MIRNA
##what transcription factors does miR-124-3p regulate? 
class TestFindTfMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfMirna1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna_arg = ekb_from_text('miR-124-3p')
        print(mirna_arg, '\n')
        content = KQMLList('FIND-TF-MIRNA')
        content.set('miRNA', KQMLString(mirna_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 156, output

##what transcription factors does miR-200c regulate? 
class TestFindTfMirna11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfMirna11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna_arg = ekb_from_text('miR-200c')
        print(mirna_arg, '\n')
        content = KQMLList('FIND-TF-MIRNA')
        content.set('miRNA', KQMLString(mirna_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        print("len(output.get('clarification'))=", str(len(output.get('clarification'))))
        assert len(output.get('clarification')) == 5, output
        assert len(output.get('clarification').get('as')) == 2, output
        
##what transcription factors does miR-200c-3p regulate? 
class TestFindTfMirna12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfMirna12, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna_arg = ekb_from_text('miR-200c-3p')
        print(mirna_arg, '\n')
        content = KQMLList('FIND-TF-MIRNA')
        content.set('miRNA', KQMLString(mirna_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 39, output
        
#subsequent query
##which of these transcription factors does miR-200c-3p regulate? 
class TestFindTfMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfMirna2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna_arg = ekb_from_text('miR-200c-3p')
        #print(mirna_arg, '\n')
        of_those = ekb_kstring_from_text('ATRX, DNMT3B, MBD5,stat3,ZMAT3')
        content = KQMLList('FIND-TF-MIRNA')
        content.set('miRNA', KQMLString(mirna_arg))
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 3, output
        
##which of these transcription factors does miR-200c-3p regulate? 
class TestFindTfMirna21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfMirna21, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna_arg = ekb_from_text('miR-200c-3p')
        #print(mirna_arg, '\n')
        of_those = 'ATRX, DNMT3B, MBD5,stat3,ZMAT3'
        content = KQMLList('FIND-TF-MIRNA')
        content.set('miRNA', KQMLString(mirna_arg))
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 4, output


################################################################################
#FIND-REGULATION
###what regulate myc? 
class TestFindRegulation1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('myc')
        print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-db') is not None:
            print("len(output.get('regulators').get('tf-db'))=", str(len(output.get('regulators').get('tf-db'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 12, output

###what increase the amount of myc? 
class TestFindRegulation11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('myc')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-db') is not None:
            print("len(output.get('regulators').get('tf-db'))=", str(len(output.get('regulators').get('tf-db'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 10, output
    
###what regulate bcl2? 
class TestFindRegulation12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation12, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('bcl2')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        assert len(output.get('regulators')) == 12, output
    
##what regulate frizzled8? 
class TestFindRegulation2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('frizzled8')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        print("len(output.get('regulators').get('tf-db'))=", str(len(output.get('regulators').get('tf-db'))))
        if output.get('regulators').get('tf-literature'):
            print("len(output.get('regulators').get('tf-literature'))=", 
                str(len(output.get('regulators').get('tf-literature'))))
        if output.get('regulators').get('gene-literature'):
            print("len(output.get('regulators').get('gene-literature'))=", 
                str(len(output.get('regulators').get('gene-literature'))))
        if output.get('regulators').get('mirna-literature'):
            print("len(output.get('regulators').get('mirna-literature'))=", 
                str(len(output.get('regulators').get('miRNA-literature'))))
        if output.get('regulators').get('other-literature'):
            print("len(output.get('regulators').get('other-literature'))=", 
                str(len(output.get('regulators').get('other-literature'))))
        print("type(output.get('regulators'))=",type(output.get('regulators')))
        assert len(output.get('regulators')) == 10, output

###what regulate GLUL? 
class TestFindRegulation3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('GLUL')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-db') is not None:
            print("len(output.get('regulators').get('tf-db'))=", str(len(output.get('regulators').get('tf-db'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 10, output
        
###what regulate GFAP? 
class TestFindRegulation4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('GFAP')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-db') is not None:
            print("len(output.get('regulators').get('tf-db'))=", str(len(output.get('regulators').get('tf-db'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 10, output

#What regulate cfos from literature?
class TestFindRegulation5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('FOS')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        content.set('source', 'literature')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 8, output

###What regulates FOS from the GEO RNAi database?
class TestFindRegulation6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation6, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('FOS')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        content.set('source', 'geo rnai')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 2, output
        
#Which kinases regulate the cfos gene?
class TestFindRegulation61(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation61, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('FOS')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        content.set('regulator', 'kinase')
        print("content=", str(content))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        assert len(output.get('regulators')) == 2, output
        
#Which transcription factors regulate frizzled8?
class TestFindRegulation62(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation62, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('frizzled8')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        content.set('regulator', 'transcription factor')
        print("content=", str(content))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        #print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        assert len(output.get('regulators')) == 2, output
        
#Which transcription factors regulate frizzled8 in liver? (subtask: find-target-tf-tissue)
class TestFindRegulation63(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation63, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('frizzled8')
        tissue = 'liver'
        content = KQMLList('FIND-REGULATION')
        content.set('target', target)
        content.set('tissue', tissue)
        content.set('regulator', 'transcription factor')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('regulators') == 'NIL', output
        
###What are the regulators of MAPK14 in bladder?
class TestFindRegulation7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation7, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('MAPK14')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        #content.set('keyword', 'regulate')
        content.set('tissue', 'bladder')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        #print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        print(output.get('regulators').get('tf-db'))
        print("len(output.get('regulators'))=", len(output.get('regulators')))
        assert len(output.get('regulators')) == 2, output
        
###What are the regulators of SMURF2 in liver?
class TestFindRegulation8(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation8, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('SMURF2')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        #content.set('keyword', 'regulate')
        content.set('tissue', 'liver')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        #print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        print(output.get('regulators'))
        assert output.get('regulators') == 'NIL', output
        
#what increase the amount of mzd8? 
class TestFindRegulation81(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation81, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('fzd8')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", str(len(output.get('regulators'))))
        if output.get('regulators').get('tf-literature') is not None:
            print("len(output.get('regulators').get('tf-literature'))=", str(len(output.get('regulators').get('tf-literature'))))
        assert len(output.get('regulators')) == 8, output
        
#Which transcription factors regulate frizzled8 in liver? (subtask: find-target-tf-tissue)
class TestFindRegulation9(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('frizzled8')
        tissue = 'liver'
        content = KQMLList('FIND-REGULATION')
        content.set('target', target)
        content.set('tissue', tissue)
        content.set('regulator', 'transcription factor')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('regulators') == 'NIL', output
        
#test of-those
#Which of those regulate myc? 
class TestFindRegulation100(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation100, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('myc')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        of_those = 'AATF,ATF5,E2F4,MXD1,HIF1A,STAT3, SRF,EGFR,TRIM28'
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", len(output.get('regulators')))
        assert len(output.get('regulators')) == 6, output

###which of those increase the amount of myc? 
class TestFindRegulation101(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation101, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('myc')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        of_those = 'EGFR, GSK3A, FOXA1, GATA3, HMGA1,AKAP4, AKT3, MAP4K4, STAT3, SRF,E2F4'
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'increase')
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", len(output.get('regulators')))
        assert len(output.get('regulators')) == 4, output

#Which of those regulate cfos from literature?
class TestFindRegulation102(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindRegulation102, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('FOS')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        of_those = 'AATF,ATF5,E2F4,MXD1,HIF1A,STAT3, SRF,EGFR,TRIM28'
        content = KQMLList('FIND-REGULATION')
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        content.set('source', 'literature')
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('regulators'))=", len(output.get('regulators')))
        assert len(output.get('regulators')) == 2, output

#################################################################################
#FIND-EVIDENCE
##Show me evidence that kras regulate frizzled8? 
class TestFindEvidence1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('kras')
        target_arg = ekb_from_text('fzd8')
        print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print(output.get('evidence'))
        #print("len(output.get('evidence').get('source_api'))=", str(len(output.get('evidence').get('source_api'))))
        assert len(output.get('evidence')) == 2, output
        
##show me evidence that kras increase frizzled8? 
class TestFindEvidence2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('kras')
        target_arg = ekb_from_text('fzd8')
        print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        assert len(output.get('evidence')) == 2, output
        
##show me evidence that kras decrease frizzled8? 
class TestFindEvidence3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('kras')
        target_arg = ekb_from_text('fzd8')
        print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'decrease')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("type(output.get('evidence'))=",type(output.get('evidence')))
        print(output.get('evidence'))
        #print("len(output.get('evidence').get('source_api'))=", str(len(output.get('evidence').get('source_api'))))
        assert len(output.get('evidence')) == 2, output
        
##Show me the evidence that IL6 increase the amount of SOCS1.
class TestFindEvidence4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('il6')
        target_arg = ekb_from_text('socs1')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        assert len(output.get('evidence')) == 2, output
        
##Show me the evidence that SRF binds to the FOS gene.
class TestFindEvidence5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('SRF')
        target_arg = ekb_from_text('cfos')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'bind')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        assert len(output.get('evidence')) == 2, output

##Show me the evidence that SRF regulate FOS gene.
class TestFindEvidence6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence6, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('SRF')
        target_arg = ekb_from_text('cfos')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        assert len(output.get('evidence')) == 4, output
        
#IncreaseAmount(miR_491(), GFAP())
class TestFindEvidence7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence7, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator_arg = ekb_from_text('miR_491')
        target_arg = ekb_from_text('GFAP')
        target = get_gene_symbol(target_arg)
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', KQMLString(regulator_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', 'increase')
        print(content, '\n')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_REGULATOR_NAME', output
        
########################################################################
#FIND-GENE-TISSUE        
#what genes are expressed in liver? 
class TestFindGeneTissue1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = ekb_kstring_from_text('AKT')
        tissue = 'liver'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 1929, output
        
#among stat3,srf, kras, and hras, what genes are expressed in liver? 
class TestFindGeneTissue11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = 'stat3, srf, kras, hras'
        tissue = 'liver'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 1, output
        
#what genes are exclusively expressed in liver? 
class TestFindGeneTissue12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = 'stat3, srf, kras, hras'
        tissue = 'liver'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        #content.set('gene', gene)
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 31, output
        
#what genes are exclusively expressed in brain? 
class TestFindGeneTissue13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = 'stat3, srf, kras, hras'
        tissue = 'brain'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        #content.set('gene', gene)
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 44, output

########################################################################
#IS-TISSUE-GENE
###Is stat3 expressed in liver? 
class TestIsTissueGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('stat3')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'liver')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
###Is kras expressed in liver? 
class TestIsTissueGene2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('kras')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'liver')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
###Is kras expressed in brain? 
class TestIsTissueGene3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('kras')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'brain')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
###Is frizzled8 expressed in liver? 
class TestIsTissueGene4(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('frizzled8')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'liver')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
###Is stat3 exclusively expressed in liver? 
class TestIsTissueGene5(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('stat3')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'liver')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
###Is GYS2 exclusively expressed in liver? 
class TestIsTissueGene6(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene6, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('GYS2')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'liver')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
###Is NEUROD2 exclusively expressed in brain? 
class TestIsTissueGene7(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene7, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('NEUROD2')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'brain')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
###Is frizzled8 expressed in BRAIN? 
class TestIsTissueGene8(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene8, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('frizzled8')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'brain')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
###Is GAST expressed in stomach? 
class TestIsTissueGene9(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene9, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target_arg = ekb_from_text('GAST')
        #print(target_arg, '\n')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', KQMLString(target_arg))
        content.set('tissue', 'stomach')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
