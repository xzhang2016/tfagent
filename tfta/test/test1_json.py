from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request, agent_clj_from_text
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor
from indra.statements import Agent
from bioagents import Bioagent
from indra.sources import trips


#####################################
# Testing the following TFTA capabilities
# IS-REGULATION
# FIND-TF
# FIND-TARGET
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

def agents_clj_from_text(text):
    ekb_xml = ekb_from_text(text)
    tp = trips.process_xml(ekb_xml)
    agents = tp.get_agents()
    clj = Bioagent.make_cljson(agents)
    return clj

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
# TEST FIND-TF
#Which transcription factors regulate frizzled8? (subtask: find-target-tf)
class TestFindTf1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('frizzled8')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
    
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 2, output
        
##What transcription factors regulate SMURF2?
#What are the regulators of SMURF2? (subtask: find-target-tf)
class TestFindTf2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('SMURF2')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output

#test failures, such as using gene family as input
#What transcription factors regulate MEK? (subtask: find-target-tf)
class TestFindTf3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('MEK')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

##what transcription factors regulate SOS?
class TestFindTf4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('SOS')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output
        
#check keyword parameter
#Which transcription factors upregulate cfos?
class TestFindTf5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('cfos')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 42, output

#Which transcription factors downregulate MEK?
class TestFindTf6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('MEK')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('keyword', 'decrease')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#tissues
#Which transcription factors regulate frizzled8 in liver? (subtask: find-target-tf-tissue)
class TestFindTf7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('frizzled8')
        _get_targets(target)
        print('target=', str(target))
        
        tissue = 'liver'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('tfs') == 'NIL', output

#Which transcription factors regulate mapk14 in bladder? (subtask: find-target-tf-tissue)
class TestFindTf8(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('mapk14')
        _get_targets(target)
        print('target=', str(target))
        
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
class TestFindTf9(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('ELK1')
        _get_targets(target)
        print('target=', str(target))
        
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
class TestFindTf10(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf10, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('AKT')
        _get_targets(target)
        print('target=', str(target))
        
        tissue = 'liver'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#Which transcription factors regulate mapk14 in breast cancer?
class TestFindTf11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('mapk14')
        _get_targets(target)
        print('target=', str(target))
        
        tissue = 'breast'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'INVALID_TISSUE', output

##What transcription factors regulate insulin?
class TestFindTf12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf12, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('insulin')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", str(len(output.get('tfs'))))
        assert len(output.get('tfs')) == 31, output

##What transcription factors regulate cofilin? (which was taken as a protein, then get COF1)
class TestFindTf13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf13, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('cofilin')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#subsequent query
##What of those regulate cofilin 1? (CFL1)
class TestFindTf14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf14, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('cofilin 1')
        _get_targets(target)
        print('target=', str(target))
        
        of_those = agents_clj_from_text("RELA, TGFB1, SMARCA2, JAK1, FOS, HRAS, SHC1")
        _get_targets(of_those)
        print('of_those=', str(of_those))
        
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", len(output.get('tfs')))
        assert len(output.get('tfs')) == 4, output
        
#####################################################################################
#TEST FIND-TARGET


