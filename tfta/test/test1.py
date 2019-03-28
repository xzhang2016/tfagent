from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor

#####################################
#Testing the following TFTA functions
#IS-REGULATION
#IS-MIRNA-TARGET
#FIND-TF
#FIND-TF-PATHWAY
#FIND-COMMON-TF-GENES
#FIND-TARGET
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


#############################################################################
#TEST IS-REGULATION
##Does STAT3 regulate the c-fos gene? (subtask: is-tf-target)
class TestIsRegulation1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#protein family
#Does stat3 regulate the SMURF gene?
class TestIsRegulation11(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target_arg = ekb_from_text('SMURF')
        print('target=', target_arg)
        get_gene_symbol(target_arg)
        get_family_name(target_arg)
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', KQMLString(target_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert len(output.get('reason')) == 1, output
        print(type(output.get('reason')))

#Does stat3 regulate the MEK gene?
class TestIsRegulation12(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation12, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target_arg = ekb_from_text('MEK')
        print('target=', target_arg)
        get_gene_symbol(target_arg)
        get_family_name(target_arg)
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', KQMLString(target_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert len(output.get('reason')) == 1, output
        print(type(output.get('reason')))
        assert output.get('reason').to_string() == \
        "((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))", output.get('reason')
        
#Does stat3 regulate the AKT gene?
#AKT is in ['ONT::PROTEIN', 'ONT::GENE', 'ONT::PROTEIN-FAMILY']
class TestIsRegulation13(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation13, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target_arg = ekb_from_text('AKT')
        print('target=', target_arg)
        get_gene_symbol(target_arg)
        get_family_name(target_arg)
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', KQMLString(target_arg))
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output

##Does MEK regulate the c-fos gene? (subtask: is-tf-target)
class TestIsRegulation14(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation14, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('MEK')
        target = ekb_kstring_from_text('c-fos')
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
##Does AKT regulate the MEK gene? (subtask: is-tf-target)
class TestIsRegulation15(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation15, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('AKT')
        target = ekb_kstring_from_text('MEK')
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
##does srf regulate acta1?
class TestIsRegulation16(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation16, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf_arg = ekb_from_text('srf')
        target_arg = ekb_from_text('acta1')
        target = get_gene_symbol(target_arg)
        content = KQMLList('IS-REGULATION')
        content.set('tf', KQMLString(tf_arg))
        content.set('target', KQMLString(target_arg))
        #print(content, '\n')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('result'))=", str(len(output.get('result'))))
        print('result=', output.get('result'))
        print('db=', output.get('db'))
        print('literature=', output.get('literature'))
        assert output.get('result') == 'TRUE', output
        assert output.get('db') == 'TRUE', output
        
#Does STAT3 regulate the JUN gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #since trips parser cannot recognize JUN as a gene, so instead try
        #c-JUN or JUN gene
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-JUN')
        tissue = 'lung'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#Does STAT3 regulate the c-fos gene in the liver? (subtask: is-tf-target-tissue)
class TestIsRegulation21(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        tissue = 'liver'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'FALSE', output
        
#Does STAT3 regulate the cfos gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation22(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        tissue = 'lung'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#TEST failures
#Does STAT3 regulate the MEK gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation23(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation23, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('MEK')
        tissue = 'lung'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#Does AKT regulate the c-fos gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation24(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation24, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('AKT')
        target = ekb_kstring_from_text('c-fos')
        tissue = 'lung'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'FALSE', output
        
#Does STAT3 regulate the c-fos gene in the breast cancer? (subtask: is-tf-target-tissue)
class TestIsRegulation25(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation25, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        tissue = 'breast cancer'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.gets('reason') == 'INVALID_TISSUE', output

#keyword
#Does STAT3 increase transcription of the c-fos gene?
class TestIsRegulation3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf_arg = ekb_from_text('stat3')
        target_arg = ekb_from_text('c-fos')
        target = get_gene_symbol(target_arg)
        keyword = 'increase'
        content = KQMLList('IS-REGULATION')
        content.set('tf', KQMLString(tf_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', keyword)
        #print(content, '\n')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('result'))=", str(len(output.get('result'))))
        print('result=', output.get('result'))
        print('db=', output.get('db'))
        print('literature=', output.get('literature'))
        assert output.get('result') == 'TRUE', output
        
#Does STAT3 upregulate MEK gene?
class TestIsRegulation31(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation31, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf_arg = ekb_from_text('stat3')
        target_arg = ekb_from_text('MEK')
        print(target_arg)
        target = get_gene_symbol(target_arg)
        keyword = 'increase'
        content = KQMLList('IS-REGULATION')
        content.set('tf', KQMLString(tf_arg))
        content.set('target', KQMLString(target_arg))
        content.set('keyword', keyword)
        #print(content, '\n')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("len(output.get('reason'))=", len(output.get('reason')))
        assert len(output.get('reason'))==1, output

######################################################################################
#IS-MIRNA-TARGET
#Does miR-20b-5p target STAT3? (subtask: is-mirna-target)
class TestIsMirnaTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaTarget1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3')
        mirna = ekb_kstring_from_text('miR-20b-5p')
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('is-miRNA-target') == 'TRUE', output

#protein family
#Does miR-20b-5p target MEK? (subtask: is-mirna-target)
class TestIsMirnaTarget11(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaTarget11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        mirna = ekb_kstring_from_text('miR-20b-5p')
        print('mirna=', str(mirna))
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#clarification test
#Does miR-200c target STAT3? (subtask: is-mirna-target)
class TestIsMirnaTarget2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaTarget2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3')
        mirna = ekb_kstring_from_text('miR-200c')
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        assert len(output.get('clarification').get('as')) == 2, output

####################################################################################
#TEST FIND-TF
#Which transcription factors regulate frizzled8? (subtask: find-target-tf)
class TestFindTf1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('frizzled8')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
    
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 2, output
        
##What transcription factors regulate SMURF2?
#What are the regulators of SMURF2? (subtask: find-target-tf)
class TestFindTf11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SMURF2')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output
        
#test failures, such as using gene family as input
#What transcription factors regulate MEK? (subtask: find-target-tf)
class TestFindTf12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf12, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
##What transcription factors regulate AKT? (subtask: find-target-tf)
class TestFindTf13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf13, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('AKT')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tfs') == 'NIL', output
        
##what transcription factors regulate SOS?
class TestFindTf14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf14, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SOS')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        #assert output.get('reason').to_string() == \
        #'((:for V106700 :error FAMILY_NAME_NOT_ALLOWED))', output

#check keyword parameter
#Which transcription factors upregulate cfos?
class TestFindTf15(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf15, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cfos')
        print(target)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 42, output
        
#Which transcription factors downregulate cfos?
class TestFindTf16(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf16, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cfos')
        print(target)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('keyword', 'decrease')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 24, output

#Which transcription factors downregulate MEK?
class TestFindTf17(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf17, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        print(target)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('keyword', 'decrease')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output

#tissues
#Which transcription factors regulate frizzled8 in liver? (subtask: find-target-tf-tissue)
class TestFindTf2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('frizzled8')
        tissue = 'liver'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tfs') == 'NIL', output
        
#Which transcription factors regulate mapk14 in bladder? (subtask: find-target-tf-tissue)
class TestFindTf21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('mapk14')
        tissue = 'bladder'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 4, output
        
#Which transcription factors regulate ELK1 in brain? (subtask: find-target-tf-tissue)
class TestFindTf22(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('ELK1')
        tissue = 'brain'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 55, output
        
#test failures
#Which transcription factors regulate AKT in liver?
class TestFindTf23(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf23, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('AKT')
        print(target)
        tissue = 'liver'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tfs') == 'NIL', output

#Which transcription factors regulate MEK in liver?
class TestFindTf24(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf24, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('MEK')
        print(target)
        tissue = 'liver'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#Which transcription factors regulate mapk14 in breast cancer?
class TestFindTf25(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf25, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('mapk14')
        #print(target)
        tissue = 'breast'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'INVALID_TISSUE', output
        
##What transcription factors regulate insulin?
class TestFindTf3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('insulin')
        print(target, '\n')
        target_arg = ekb_from_text('insulin')
        get_gene_symbol(target_arg)
        content = KQMLList('FIND-TF')
        content.set('target', target)
       
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 31, output
        
##What transcription factors regulate cofilin? (which was taken as a protein, then get COF1)
class TestFindTf31(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf31, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cofilin')
        print(target, '\n')
        target_arg = ekb_from_text('cofilin')
        get_gene_symbol(target_arg)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tfs') == 'NIL', output
        
##What transcription factors regulate cofilin gene? (which was taken as a protein, then get COF1)
class TestFindTf32(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf32, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cofilin gene')
        print(target, '\n')
        target_arg = ekb_from_text('cofilin gene')
        get_gene_symbol(target_arg)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tfs') == 'NIL', output
        
##What transcription factors regulate cofilin 1? (CFL1)
class TestFindTf33(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf33, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cofilin 1')
        print(target, '\n')
        target_arg = ekb_from_text('cofilin 1')
        get_gene_symbol(target_arg)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 17, output

#subsequent query
##What of those regulate cofilin 1? (CFL1)
class TestFindTf4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cofilin 1')
        of_those = "RELA, TGFB1, SMARCA2, JAK1, FOS, HRAS,SHC1"
        print(target, '\n')
        target_arg = ekb_from_text('cofilin 1')
        get_gene_symbol(target_arg)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", len(output.get('tfs')))
        assert len(output.get('tfs')) == 4, output

##What of RELA, TGFB1, SMARCA2, JAK1, FOS, AND HRAS regulate cofilin 1? (CFL1)
class TestFindTf41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf41, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('cofilin 1')
        of_those = ekb_kstring_from_text('RELA, TGFB1, SMARCA2, JAK1, FOS, HRAS,SHC1')
        print(target, '\n')
        target_arg = ekb_from_text('cofilin 1')
        get_gene_symbol(target_arg)
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", len(output.get('tfs')))
        assert len(output.get('tfs')) == 4, output

###################################################################################
#FIND-TF-PATHWAY
#What transcription factors are in the calcium regulated pathways? (subtask: find-tf-keyword)
class TestFindTfPathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = ekb_kstring_from_text('calcium regulated pathways')
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 3, output
        
#Which transcription factors are in the MAPK signaling pathway? (subtask: find-tf-pathway)
class TestFindTfPathway11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway11, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 14, output
        
#subsequent query
#which of these transcription factors are in the calcium regulated pathways? (subtask: find-tf-keyword)
class TestFindTfPathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = ekb_kstring_from_text('calcium regulated pathways')
        of_those = 'STAT3,SRF,NFAT5,JAK1'
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', keyword)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert output.get('pathways') == 'NIL', output
        
#Which of those transcription factors are in the MAPK signaling pathway? (subtask: find-tf-pathway)
class TestFindTfPathway21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway21, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        of_those = ekb_kstring_from_text('STAT3,SRF,NFAT5,JAK1')
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 3, output


###############################################################################
#FIND-COMMON-TF-GENES
#What transcription factors are shared by the SRF, HRAS, and elk1 genes? (subtask: find-common-tf-genes)
class TestFindCommonTfGenes1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SRF, HRAS, elk1')
        #count = 'count'
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output
        
#What transcription factors are in common to the STAT3, SOCS3, and CREB5 genes?
class TestFindCommonTfGenes11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, CREB5')
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output

#What transcription factors are in common to the STAT3, SOCS3, IFNG, FOXO3, and CREB5 genes?        
class TestFindCommonTfGenes12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, IFNG, FOXO3, SOCS3, CREB5')
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('tfs'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 8, output
        
#test gene family
#What transcription factors are in common to the STAT3, SOCS3, and MEK genes?
#MEK will be ignored in this case
class TestFindCommonTfGenes13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, MEK')
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 1, output
        
#What transcription factors are in common to the STAT3, SOCS3, and AKT genes?
class TestFindCommonTfGenes14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes14, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, AKT')
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 1, output

#subsequent query
#Which of these transcription factors are shared by the SRF, HRAS, FOS, and elk1 genes? (subtask: find-common-tf-genes)
class TestFindCommonTfGenes2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SRF, HRAS, cfos, elk1')
        of_those = 'stat3,ELK1,TFAP2A,CREB1,TP53'
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", len(output.get('tfs')))
        assert len(output.get('tfs')) == 4, output
        
#Which of these transcription factors are shared by the SRF, HRAS, FOS, and elk1 genes? (subtask: find-common-tf-genes)
class TestFindCommonTfGenes21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SRF, HRAS, cfos, elk1')
        of_those = ekb_kstring_from_text('stat3,ELK1,TFAP2A,CREB1,TP53')
        print('of_those=', of_those)
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", len(output.get('tfs')))
        assert len(output.get('tfs')) == 3, output


#####################################################################################
#TEST FIND-TARGET
#what genes are regulated by smad2? (subtask: find-tf-target)
class TestFindTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('smad2')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 117, output
        
#what genes are regulated by elk1 and smad2? (subtask: find-tf-target)
class TestFindTarget11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('ELK1, SMAD2')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 28
        
#what genes are regulated by MEK (subtask: find-tf-target)
class TestFindTarget13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('MEK')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#what genes are regulated by AKT (subtask: find-tf-target)
class TestFindTarget14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget14, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('akt')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('targets') == 'NIL', output
        
#What genes does stat3 regulate in lung? (subtask: find-tf-target-tissue)
class TestFindTarget2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        tissue = 'lung'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 57, output
        
#what genes does stat3 regulate in liver? (subtask: find-tf-target-tissue)
class TestFindTarget21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        tissue = 'liver'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('targets') == 'NIL', output
        
##what genes does MEK regulate in liver? (subtask: find-tf-target-tissue)
class TestFindTarget22(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('MEK')
        tissue = 'liver'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#what genes does AKT regulate in liver? (subtask: find-tf-target-tissue)
class TestFindTarget23(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget23, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('AKT')
        tissue = 'liver'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('targets') == 'NIL', output
        
#What genes does smad2 upregulate?
class TestFindTarget3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('smad2')
        keyword = 'increase'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 77, output
        
#What genes does smad2 downregulate?
class TestFindTarget31(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget31, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('smad2')
        keyword = 'decrease'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 48, output
        
#What genes does MEK downregulate?
class TestFindTarget32(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget32, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('MEK')
        keyword = 'decrease'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').to_string() == \
        '((:for V14939346 :error FAMILY_NAME_NOT_ALLOWED))', output
        
#test sequencing query
#Which of these genes are regulated by STAT3? (subtask: find-tf-target)
class TestFindTarget4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        of_those = "A2M, ABCA2, AKAP12, AKT1, PBRM1, SMAD2,CEBPA"
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('of-those', of_those)
        print("content=", content)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 5, output
        
#Which of A2M, ABCA2, AKAP12, AKT1, PBRM1, SMAD2, and CEBPA genes are regulated by STAT3? (subtask: find-tf-target)
#A2M is under ONT::MUTATION
class TestFindTarget41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget41, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        of_those = ekb_kstring_from_text('A2M, ABCA2, AKAP12, AKT1, PBRM1, SMAD2,CEBPA')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('of-those', of_those)
        print("content=", content)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 4, output

#target-type test
#What transcription factors are regulated by stat3?
class TestFindTarget5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = 'STAT3'
        target_type = 'tf'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('target-type', target_type)
        #print("content=", content)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 135, output
        
#What genes are regulated by stat3?
class TestFindTarget51(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget51, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = 'STAT3'
        #target_type = 'tf'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        #content.set('target-type', target_type)
        #print("content=", content)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 423, output

#What TFs does stat3 regulate in lung? (subtask: find-tf-target-tissue)
class TestFindTarget6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = 'STAT3'
        tissue = 'lung'
        target_type = 'tf'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        content.set('target-type', 'tf')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 55, output
        
#What TFs does smad2 upregulate?
class TestFindTarget7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('smad2')
        keyword = 'increase'
        target_type = 'tf'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('keyword', keyword)
        content.set('target-type', 'tf')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 22, output

#What kinases are regulated by stat3?
class TestFindTarget51(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget51, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = 'STAT3'
        target_type = 'kinase'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('target-type', target_type)
        #print("content=", content)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 12, output
        

#What kinases does stat3 regulate in lung? (subtask: find-tf-target-tissue)
class TestFindTarget61(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget61, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = 'STAT3'
        tissue = 'lung'
        target_type = 'kinase'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        content.set('target-type', target_type)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert output.get('targets') == 'NIL', output
        
#What kinases does smad2 upregulate?
class TestFindTarget71(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget71, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('smad2')
        keyword = 'increase'
        target_type = 'kinase'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('keyword', keyword)
        content.set('target-type', target_type)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", str(len(output.get('targets'))))
        assert len(output.get('targets')) == 6, output
        
