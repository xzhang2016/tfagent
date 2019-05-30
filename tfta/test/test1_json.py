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
# Testing the following TFTA capabilities
# IS-REGULATION
# FIND-TF
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
class TestIsRegulation2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation2, self).__init__(TFTA_Module)

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
        
#Does stat3 regulate the MEK gene?
class TestIsRegulation3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        print('tf=', str(tf))
        target_arg = agent_clj_from_text('MEK')
        _get_targets(target_arg)
        print('target=', str(target_arg))
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
        
#Does stat3 regulate the AKT gene?
class TestIsRegulation4(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        print('tf=', str(tf))
        target_arg = agent_clj_from_text('AKT')
        _get_targets(target_arg)
        print('target=', str(target_arg))
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

##Does AKT regulate the MEK gene? (subtask: is-tf-target)
class TestIsRegulation5(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('AKT')
        _get_targets(tf)
        print('tf=', str(tf))
        target = agent_clj_from_text('MEK')
        _get_targets(target)
        print('target=', str(target))
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert output.get('reason') == 'FAMILY_NAME', output
        assert len(output.get('clarification')) == 5, output

##does srf regulate acta1?
class TestIsRegulation6(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation6, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('srf')
        _get_targets(tf)
        print('tf=', str(tf))
        target = agent_clj_from_text('acta1')
        _get_targets(target)
        print('target=', str(target))
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
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
class TestIsRegulation7(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #since trips parser cannot recognize JUN as a gene, so instead try
        #c-JUN or JUN gene
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        print('tf=', str(tf))
        target = agent_clj_from_text('c-JUN')
        _get_targets(target)
        print('target=', str(target))
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
class TestIsRegulation8(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        print('tf=', str(tf))
        target = agent_clj_from_text('c-fos')
        _get_targets(target)
        print('target=', str(target))
        tissue = 'liver'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'FALSE', output

#Does AKT regulate the c-fos gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation9(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('AKT')
        _get_targets(tf)
        print('tf=', str(tf))
        
        target = agent_clj_from_text('c-fos')
        _get_targets(target)
        print('target=', str(target))
        
        tissue = 'lung'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#Does STAT3 regulate the c-fos gene in the breast cancer? (subtask: is-tf-target-tissue)
class TestIsRegulation10(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation10, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('STAT3')
        _get_targets(tf)
        print('tf=', str(tf))
        
        target = agent_clj_from_text('c-fos')
        _get_targets(target)
        print('target=', str(target))
        
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
class TestIsRegulation11(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = agent_clj_from_text('stat3')
        _get_targets(tf)
        print('tf=', str(tf))
        
        target = agent_clj_from_text('c-fos')
        _get_targets(target)
        print('target=', str(target))
        
        keyword = 'increase'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('keyword', keyword)
        print(content, '\n')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('result'))=", str(len(output.get('result'))))
        print('result=', output.get('result'))
        print('db=', output.get('db'))
        print('literature=', output.get('literature'))
        assert output.get('result') == 'TRUE', output
        assert output.get('literature') == 'TRUE', output
        assert output.get('db') == 'FALSE', output

##############################################################################
#

