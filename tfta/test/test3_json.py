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
import re

#####################################
# Testing the following TFTA capabilities
# IS-MIRNA-TARGET
# FIND-TARGET-MIRNA
# FIND-MIRNA-TARGET
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

def agents_clj_from_text(text):
    ekb_xml = ekb_from_text(text)
    tp = trips.process_xml(ekb_xml)
    agents = tp.get_agents()
    clj = Bioagent.make_cljson(agents)
    return clj

def get_mirnas(mir_arg):
    mirna = []
    try:
        agents = Bioagent.get_agent(mir_arg)
    except Exception:
        return []
    if isinstance(agents, list):
        mirna = [_get_mirna_name(a.name) for a in agents if a is not None]
    elif isinstance(agents, Agent):
        mirna = [_get_mirna_name(agents.name)]
    if mirna:
        print('mirna=', ','.join(mirna))
    else:
        print('mirna = None')
    return mirna
    
def _get_mirna_name(str1):
    plist = re.findall('([0-9]+-[a-zA-Z])', str1)
    s = str1
    for p in plist:
        p1 = p.replace('-','')
        s = s.replace(p, p1)
    return s

#In order to byparse trips due to the wrong interpretation of some microRNA names
def make_mirna_cljson(mir_str):
    mir_str = mir_str.replace(' ', '')
    mir_list = mir_str.split(',')
    mir_agent = [Agent(mir) for mir in mir_list]
    mir_json = Bioagent.make_cljson(mir_agent)
    return mir_json
    
#############################################################################
# IS-MIRNA-TARGET
#Does miR-20b-5p target STAT3? (subtask: is-mirna-target)
class TestIsMirnaTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaTarget1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('STAT3')
        mirna = agent_clj_from_text('miR-20b-5p')
        agents = Bioagent.get_agent(mirna)
        print('name=', agents.name)
        print('db_refs=', agents.db_refs)
        
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('is-miRNA-target') == 'TRUE', output

#protein family
#Does miR-20b-5p target MEK? (subtask: is-mirna-target)
class TestIsMirnaTarget2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaTarget2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('MEK')
        mirna = agent_clj_from_text('miR-20b-5p')
        agents = Bioagent.get_agent(mirna)
        print('name=', agents.name)
        print('db_refs=', agents.db_refs)
        
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output
        
#clarification test
#Does miR-200c target STAT3? (subtask: is-mirna-target)
class TestIsMirnaTarget3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaTarget3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('STAT3')
        mirna = agent_clj_from_text('miR-200c')
        agents = Bioagent.get_agent(mirna)
        print('name=', agents.name)
        print('db_refs=', agents.db_refs)
        
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        assert len(output.get('clarification').get('as')) == 2, output

####################################################################################
#FIND-TARGET-MIRNA
#What genes does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('miR-20b-5p')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 917, output

#What genes are regulated by miR-297? (subtask: find-target-mirna)
class TestFindTargetMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('miR-297')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 190, output
        
#What genes are regulated by miR-20b-5p and MIR-29B-1-5P? (subtask: find-target-mirna)
class TestFindTargetMirna3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agents_clj_from_text('miR-20b-5p, MIR-29B-1-5P')
        print("mirna=", mirna)
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 12, output

#What are the genes that have strong evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('mir-122-5p')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'strong')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 69, output
        
#What are the genes that have weak evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('mir-122-5p')
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
class TestFindTargetMirna6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('mir-128')
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
class TestFindTargetMirna7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('miR-200C')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        print("output.get('reason') = ", output.get('targets'))
        assert output.get('reason') == 'MIRNA_NOT_FOUND', output
        assert len(output.get('clarification').get('as')) == 2, output

#Which of those genes does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna8(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('miR-20b-5p')
        of_those = agents_clj_from_text("STAT3, SRF, HRAS, KRAS, ELK1, JAK1, JAK2, FOS")
        _get_targets(of_those)
        print('target=', str(of_those))
        
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 2, output
        
#What are the genes that have weak evidence of being regulated by mir-122-5p.? (subtask: find-target-mirna)
class TestFindTargetMirna9(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('mir-122-5p')
        of_those = agents_clj_from_text("STAT3, SRF, CDK4, CDK19, CSRP1, DZIP1L, HSPA4L")
        _get_targets(of_those)
        print('target=', str(of_those))
        
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('strength', 'weak')
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets')) = ", len(output.get('targets')))
        assert len(output.get('targets')) == 5, output

#test target-type
#What kinases does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna10(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna10, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('miR-20b-5p')
        target_type = 'kinase'
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('target-type', target_type)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 40, output
        
#What transcription facotrs does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTargetMirna11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTargetMirna11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = agent_clj_from_text('miR-20b-5p')
        target_type = 'transcription factor'
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        content.set('target-type', target_type)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('targets'))=", len(output.get('targets')))
        assert len(output.get('targets')) == 130, output

######################################################################################
#FIND-MIRNA-TARGET
#What microRNAs target STAT3? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('STAT3')
        content = KQMLList('FIND-MIRNA-TARGET')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('miRNAs'))=", len(output.get('miRNAs')))
        assert len(output.get('miRNAs')) == 80, output
        
#What microRNAs target il2? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('il2')
        content = KQMLList('FIND-MIRNA-TARGET')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('miRNAs'))", len(output.get('miRNAs')))
        assert len(output.get('miRNAs')) == 5, output

##What microRNAs target MEK? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('MEK')
        content = KQMLList('FIND-MIRNA-TARGET')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#test of-those
#Which of those mirnas target il2?
class TestFindMirna4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('il2')
        of_those = agents_clj_from_text('miR-20b-5p, miR-125b-5p, miR-337-3p, miR-155-5p, miR-877-3p, miR-181c-5p, miR-503-3p')
        get_mirnas(of_those)
        print(of_those)
        
        content = KQMLList('FIND-MIRNA-TARGET')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('miRNAs'))=", len(output.get('miRNAs')))
        assert len(output.get('miRNAs')) == 2, output
        
#Which of those mirnas also target stat3? (use output from what micrornas target il2?)
class TestFindMirna5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('stat3')
        of_those = agents_clj_from_text('miR-181c-5p, miR-484, MIR-155-5p, let-7i-5p, miR-503-3p')
        get_mirnas(of_those)
        print(of_those)
        
        content = KQMLList('FIND-MIRNA-TARGET')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('miRNAs'))=", len(output.get('miRNAs')))
        assert len(output.get('miRNAs')) == 1, output

#####################################################################################
#


