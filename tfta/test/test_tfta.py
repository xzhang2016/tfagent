from kqml import KQMLList
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest


class TestFindTfTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfTarget1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('SMAD2')
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


class TestFindTfTarget2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfTarget2, self).__init__(TFTA_Module)

    def create_message(self):
        tf = ekb_kstring_from_text('ELK1, SMAD2')
        content = KQMLList('FIND-TF-TARGET')
        content.set('tf', tf)
        return get_request(content), content

    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 28
