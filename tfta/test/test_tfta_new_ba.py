from kqml import KQMLList
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest

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
        
#Does STAT3 regulate the JUN gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('JUN')
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
        
#Does miR-20b-5p target STAT3? (subtask: is-mirna-target)
class TestIsRegulation3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3')
        mirna = ekb_kstring_from_text('miR-20b-5p')
        content = KQMLList('IS-REGULATION')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('is-miRNA-target') == 'TRUE', output
        
#
    
