from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor

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

#################################################################################
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
class TestFindPathway3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway3, self).__init__(TFTA_Module)
        
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
class TestFindPathway32(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway32, self).__init__(TFTA_Module)
        
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
        assert len(output.get('pathways')) == 74, output
        
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
        assert len(output.get('reason')) == 1, output.get('reason')
        #assert output.get('reason').to_string() == \
        #'((:for V35246 :error FAMILY_NAME_NOT_ALLOWED))', output
        
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

#########################################################################
#functions
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

def get_term_id(target_arg):
    term_id = []
    ont1 = ['ONT::PROTEIN-FAMILY', 'ONT::GENE-FAMILY']
    tp = TripsProcessor(target_arg)
    for term in tp.tree.findall('TERM'):
        if term.find('type').text in ont1:
            term_id.append(term.attrib['id'])
    print('term_id = ' + ','.join(term_id))
    return ','.join(term_id)
    
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

