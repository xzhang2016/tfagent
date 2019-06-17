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
# FIND-PATHWAY
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


###########################################################
# FIND-PATHWAY
#What pathways involve SRF? (subtask: find-pathway-gene)
class TestFindPathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('SRF')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 27, output

#What pathways involve kras and elk1? (subtask: find-pathway-gene)
class TestFindPathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('KRAS, ELK1')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 23, output
        
##What pathways involve MEK? (subtask: find-pathway-gene)
class TestFindPathway3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('MEK')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        print("len(output.get('clarification').get('as')=", len(output.get('clarification').get('as')))
        assert len(output.get('clarification')) == 5, output
        assert len(output.get('clarification').get('as')) == 2, output

#What pathways involve stat and il2? (subtask: find-pathway-gene)
class TestFindPathway4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('stat, il2')
      
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        print("len(output.get('clarification').get('as')=", len(output.get('clarification').get('as')))
        assert len(output.get('clarification')) == 5, output
        assert len(output.get('clarification').get('as')) == 8, output

#test regulator
#What pathways involve genes regulated by smad2? (subtask: find-pathway-gene)
class TestFindPathway5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('smad2')
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 13, output

#Which Reactome pathways utilize SRF? (subtask: find-pathway-db-gene)
class TestFindPathway6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('SRF')
        db = 'Reactome'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#Which KEGG pathways utilize SRF? (subtask: find-pathway-db-gene)
class TestFindPathway7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('SRF')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
##Which KEGG pathways utilize MEK? (subtask: find-pathway-db-gene)
class TestFindPathway8(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('MEK')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        print("len(output.get('clarification').get('as')=", len(output.get('clarification').get('as')))
        assert len(output.get('clarification')) == 5, output
        assert len(output.get('clarification').get('as')) == 2, output

#test regulator for find-pathway-db-gene
#which kegg pathways utilize genes regulated by smad2?
class TestFindPathway9(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('smad2')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output
        
#What pathways involve calcium? (subtask: find-pathway-keyword)
class TestFindPathwayKeyword10(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword10, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 9, output

#What pathways involve immune system?        
class TestFindPathwayKeyword11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword11, self).__init__(TFTA_Module)
        
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
        
#What immune pathways involve kras and elk1? (subtask: find-pathway-gene-keyword)
class TestFindPathway12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('kras, elk1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#What immune pathways involve tap1 and jak1? (subtask: find-pathway-gene-keyword)
class TestFindPathway13(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway13, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('tap1, jak1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 1, output
        
##What immune pathways involve AKT and jak1? (subtask: find-pathway-gene-keyword)
class TestFindPathway14(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway14, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('AKT, jak1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        print("len(output.get('clarification').get('as')=", len(output.get('clarification').get('as')))
        assert len(output.get('clarification')) == 5, output
        assert len(output.get('clarification').get('as')) == 3, output

#test regulator for find-pathway-gene-keyword
#Which immune system pathways involve genes regulated by stat3?
class TestFindPathway15(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway15, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('STAT3')
        keyword = 'immune-system'
        
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output

#test the limit (find-pathway-keyword)
class TestFindPathway16(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway16, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'signaling-pathway'
        
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_PATHWAY_NAME', output

###########################################################################
# 


