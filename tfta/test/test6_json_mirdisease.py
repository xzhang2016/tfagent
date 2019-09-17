from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
from tfta.mirDisease import mirDisease
from tfta.tfta_module import TFTA_Module
from enrichment.GO import GOEnrich
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request, agent_clj_from_text
from bioagents.tests.integration import _IntegrationTest, _FailureTest
from indra.sources.trips.processor import TripsProcessor
from indra.statements import Agent
from bioagents import Bioagent
from indra.sources import trips
import re

#######################################################################################
# IS-MIRNA-DISEASE
#
#
########################################################################################


def _get_targets(target):
        proteins = None
        family = None
        agents = Bioagent.get_agent(target)
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
            print('Genes = None')
        if family:
            print('family=', ','.join(family))
        else:
            print('family = None')
        return proteins,family
        
def get_genes(target):
    agents = Bioagent.get_agent(target)
    genes = [agent.name for agent in agents]
    print('genes=', ','.join(genes))
    
def make_json(text):
    genes = text.split(',')
    agents = [Agent(g, db_refs={'HGNC':13}) for g in genes]
    clj = Bioagent.make_cljson(agents)
    return clj

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

##################################################################################
# IS-MIRNA-DISEASE
# Is mir-218-1 associated with urinary bladder cancer?
class TestIsMirnaDisease1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaDisease1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mir = make_mirna_cljson('mir-218-1')
        disease = 'Urinary Bladder Cancer'
        content = KQMLList('IS-MIRNA-DISEASE')
        content.set('mirna', mir)
        content.set('disease', disease)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

# Is mir-218-5p associated with urinary bladder cancer?
class TestIsMirnaDisease2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsMirnaDisease2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mir = make_mirna_cljson('mir-218-5p')
        disease = 'Urinary Bladder Cancer'
        content = KQMLList('IS-MIRNA-DISEASE')
        content.set('mirna', mir)
        content.set('disease', disease)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

##########################################################################
# FIND-MIRNA-DISEASE
# Which mirnas are associated with urinary bladder cancer?
class TestFindMirnaDisease1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirnaDisease1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        disease = 'Urinary Bladder Cancer'
        content = KQMLList('FIND-MIRNA-DISEASE')
        content.set('disease', disease)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('miRNAs'))=", len(output.get('miRNAs')))
        assert len(output.get('miRNAs')) == 132, output
        
    
# Which of those mirnas are associated with breast cancer?
class TestFindMirnaDisease2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirnaDisease2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        disease = 'breast Neoplasms'
        of_those = make_mirna_cljson('mir-200c, mir-195,mir-188,mir-181a')
        content = KQMLList('FIND-MIRNA-DISEASE')
        content.set('disease', disease)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('miRNAs'))=", len(output.get('miRNAs')))
        assert len(output.get('miRNAs')) == 3, output
        
##################################################################################
# FIND-DISEASE-MIRNA
#what diseases are associated with mir-200c-5p?
class TestFindDiseaseMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindDiseaseMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = make_mirna_cljson('mir-200c-5p')
        content = KQMLList('FIND-DISEASE-MIRNA')
        content.set('mirna', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('disease'))=", len(output.get('disease')))
        assert len(output.get('disease')) == 3, output

#what diseases are associated with mir-200c and mir-195?
class TestFindDiseaseMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindDiseaseMirna2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = make_mirna_cljson('mir-200c, mir-195')
        content = KQMLList('FIND-DISEASE-MIRNA')
        content.set('mirna', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('disease'))=", len(output.get('disease')))
        assert len(output.get('disease')) == 39, output

#what diseases can I ask about for microRNA?
class TestFindDiseaseMirna3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindDiseaseMirna3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #mirna = make_mirna_cljson('mir-200c, mir-195')
        content = KQMLList('FIND-DISEASE-MIRNA')
        #content.set('mirna', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('disease'))=", len(output.get('disease')))
        assert len(output.get('disease')) == 894, output


