from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request, agent_clj_from_text
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor
from indra.statements import Agent
from bioagents import Bioagent

#####################################
#Testing the following TFTA functions
#IS-REGULATION
#
#
#
#
#
######################################

def _get_targets(target_arg):
        proteins = None
        family = None
        agents = Bioagent.get_agent(target_arg)
        if isinstance(agents, list):
            proteins = [a.name for a in agents if a is not None and ('UP' in a.db_refs or 'HGNC' in a.db_refs)]
            family = [a.name for a in agents if a is not None and 'FPLX' in a.db_refs and a.name not in proteins]
        elif isinstance(agents, Agent):
            if 'UP' in agents.db_refs or 'HGNC' in agents.db_refs:
                proteins = [agents.name]
            if not proteins and 'FPLX' in agents.db_refs:
                family = [agents.name]
        if proteins:
            print('genes=', ','.join(proteins))
        else:
            print('Genes = None\n')
        if family:
            print('family=', ','.join(family))
        else:
            print('family = None\n')
        return proteins,family
        
#############################################################################
#TEST IS-REGULATION
##Does STAT3 regulate the c-fos gene? (subtask: is-tf-target)
class TestIsRegulation1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        target = agent_clj_from_text('c-fos')
        _get_targets(target)
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
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        target_arg = agent_clj_from_text('SMURF')
        _get_targets(target_arg)
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target_arg)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert output.get('reason') == 'FAMILY_NAME', output
        assert len(output.get('clarification')) == 5, output

