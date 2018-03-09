from kqml import KQMLList
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest


#FIND-TF-TARGET
class _TestFindTfTarget(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindTfTarget, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text(self.tf)
        content = KQMLList('FIND-TF-TARGET')
        content.set('tf', tf)
        return get_request(content), content

    def check_response_to_message(self, output):
        # Here we check the details of the response and make sure everything
        # we are looking for is in the output
        # First, we check that the response is a success
        assert output.head() == 'SUCCESS', output
        # Then we check that we got the expected number of target
        assert len(output.get('targets')) == 115
        # We could do further checks here to see if a given expected target
        # shows up in the response, etc.


#what genes are regulated by smad2?
class TestFindTfTarget1(_TestFindTfTarget):
    tf = 'SMAD2'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 115, output

#What genes does stat3 regulate?


#what genes are regulated by elk1 and smad2?
class TestFindTfTarget2(_TestFindTfTarget):
    tf = 'ELK1, SMAD2'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 28

#test IS-TF-TARGET
class _TestIsTfTarget(_IntegrationTest):
    def __init__(self, *args):
        super(_TestIsTfTarget, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text(self.tf)
        target = ekb_kstring_from_text(self.target)
        content = KQMLList('IS-TF-TARGET')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content

#Does STAT3 regulate the c-fos gene?
class TestIsTfTarget1(_TestIsTfTarget):
    tf = 'STAT3'
    target = 'c-fos'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#
