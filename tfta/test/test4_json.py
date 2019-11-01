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
# IS-PATHWAY-GENE
# FIND-GENE-PATHWAY
# FIND-KINASE-PATHWAY
# FIND-TF-PATHWAY
# find-common-pathway-genes
# FIND-PATHWAY-DB-KEYWORD
# FIND-PATHWAY-KEYWORD
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
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 27, output
        
class TestFindPathway1_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway1_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        of_those = make_pathway_json(['Cell Differentiation - Index', 'HTLV-I infection - Homo sapiens (human)', 
                                      'Mesodermal Commitment Pathway', 'il12 pathway'])
        gene = agent_clj_from_text('SRF')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output

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
        print("len(output.get('pathways'))=", len(output.get('pathways')))
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
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 50, output
        
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
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 26, output

#test regulator
#What pathways involve genes regulated by smad2? 
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
        assert len(output.get('pathways')) == 30, output
        
class TestFindPathway5_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway5_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('smad2')
        of_those = make_pathway_json(['Gene Expression', 'Generic Transcription Pathway', 
                   'Pathways in cancer - Homo sapiens (human)', 'Cell Differentiation - Index'])
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output
        
#What pathways involve genes regulated by FOS (protein family)? 
class TestFindPathway51(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway51, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('FOS')
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 30, output

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
        print("len(output.get('pathways'))=" + str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 4, output
        
class TestFindPathway6_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway6_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('SRF')
        db = 'Reactome'
        of_those = make_pathway_json(['Cell Differentiation - Index', 'RHO GTPase Effectors',
                   'Mesodermal Commitment Pathway', 'Signaling by Rho GTPases'])
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=" + str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 2, output
        
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
        print("len(output.get('pathways'))=" + str(len(output.get('pathways'))))
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
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=" + str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 50, output

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
        assert len(output.get('pathways')) == 30, output
        
class TestFindPathway9_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway9_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('smad2')
        db = 'KEGG'
        of_those = make_pathway_json(['Cell Differentiation - Index', 'Pathways in cancer - Homo sapiens (human)',
                   'Transcriptional misregulation in cancer - Homo sapiens (human)', 'HTLV-I infection - Homo sapiens (human)'])
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('database', db)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output
        
#which kegg pathways utilize genes regulated by FOS?
class TestFindPathway9_2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway9_2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('FOS')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 30, output
        
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
        
class TestFindPathwayKeyword10_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword10_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        of_those = make_pathway_json(['Calcium Regulation in the Cardiac Cell', 
                   'Highly calcium permeable nicotinic acetylcholine receptors'])
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 2, output

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
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 2, output
        
class TestFindPathway12_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway12_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('kras, elk1')
        keyword = 'immune'
        of_those = make_pathway_json(['Cell Differentiation - Index', 'immune pathway',
                   'Immune System'])
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 1, output
        
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
        print("len(output.get('pathways'))=", len(output.get('pathways')))
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
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 2, output
        
##What immune pathways involve AKT? (subtask: find-pathway-gene-keyword)
class TestFindPathway141(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway141, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('AKT')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 3, output

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
        
class TestFindPathway15_1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway15_1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('STAT3')
        keyword = 'immune-system'
        of_those = make_pathway_json(['Cell Differentiation - Index', 'immune pathway',
                   'Immune System', 'Innate Immune System'])
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('keyword', keyword)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output

#Which immune system pathways involve genes regulated by FOS?
class TestFindPathway16(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway16, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('FOS')
        keyword = 'immune-system'
        
        content = KQMLList('FIND-PATHWAY')
        content.set('regulator', regulator)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print("len(output.get('pathways'))=" + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output

#test the limit (find-pathway-keyword)
class TestFindPathway17(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway17, self).__init__(TFTA_Module)
        
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
# IS-PATHWAY-GENE
#Does the mTor pathway utilize SGK1? (subtask: is-pathway-gene)
class TestIsPathwayGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('SGK1')
        pathway = agent_clj_from_text('mTor pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 1, output
        
##Does the mTor pathway utilize MEK? (subtask: is-pathway-gene)
class TestIsPathwayGene2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('MEK')
        pathway = agent_clj_from_text('mTor pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        print("len(output.get('clarification').get('as'))=", len(output.get('clarification').get('as')))
        assert len(output.get('clarification')) == 5, output
        assert len(output.get('clarification').get('as')) == 2, output
        
#Is STAT3 in the EMT pathway?
class TestIsPathwayGene3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('stat3')
        pathway = agent_clj_from_text('EMT pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert output.get('pathways') == 'NIL', output
        
#Is TGFB1 in the EMT pathway?
class TestIsPathwayGene4(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsPathwayGene4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('TGFB1')
        pathway = agent_clj_from_text('EMT pathway')
        print(pathway)
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 1, output

####################################################################################
# FIND-GENE-PATHWAY
#What genes are in the MAPK signaling pathway? (subtask: find-gene-pathway)
class TestFindGenePathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('MAPK')
        print(pathway)
        
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 18, output

#What genes are involved in the il-12 pathway?
class TestFindGenePathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('il-12 pathway')
        print(pathway)
        
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 2, output
        
#What genes are in the immune system pathway?
class TestFindGenePathway3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('immune system')
        print(pathway)
        
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output

#subsequent query
#Which of these genes are in the immune system pathway?
class TestFindGenePathway4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGenePathway4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('immune system')
        print(pathway)
        of_genes = agents_clj_from_text('STAT3,SRF,KRAS,HRAS,FZD8,JAK1,JAK2,FOS')
        
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_genes)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output

###################################################################################
# FIND-KINASE-PATHWAY
#What kinases are in the MAPK signaling pathway?
class TestFindKinasePathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('MAPK')
        print(pathway)
        
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 18, output

#What kinases are in the immune system pathway?
class TestFindKinasePathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('immune system')
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
        pathway = agent_clj_from_text('MAPK')
        of_those = agents_clj_from_text('AKT1, FGFR1,CAMK2B,MAP2K1,RAF1,MAPK1,JAK1,MAP3K1')
        #print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 15, output
        
#What of these kinases are in the immune system pathway?
class TestFindKinasePathway4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinasePathway4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('immune system')
        of_those = agents_clj_from_text('AKT1, FGFR1,CAMK2B,MAP2K1,RAF1,MAPK1,JAK1,MAP3K1')
        #print(pathway)
        content = KQMLList('FIND-KINASE-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output

#################################################################################
# FIND-TF-PATHWAY
#What transcription factors are in the calcium regulated pathways? (subtask: find-tf-keyword)
class TestFindTfPathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('calcium regulated pathways')
        print(pathway)
        
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 3, output

#Which transcription factors are in the MAPK signaling pathway? (subtask: find-tf-pathway)
class TestFindTfPathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway2, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('MAPK')
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 12, output
        
#subsequent query
#Which of those transcription factors are in the MAPK signaling pathway? (subtask: find-tf-pathway)
class TestFindTfPathway3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTfPathway3, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('MAPK')
        of_those = agents_clj_from_text('STAT3,SRF,NFAT5,JAK1, smad2, jak2')
        
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', pathway)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 3, output

###############################################################################
# find-common-pathway-genes
#Which pathways are shared by STAT3, SOCS3, IFNG, FOXO3, and CREB5 genes? 
#(subtask: find-common-pathway-genes)
class TestFindCommonPathwayGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, SOCS3, IFNG, FOXO3, CREB5')
        print('target=', target)
       
        content = KQMLList('find-common-pathway-genes')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 30, output

#What signaling pathways are shared by STAT3, AKT and SRF? (subtask: find-common-pathway-genes)
class TestFindCommonPathwayGene2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('KRAS, MEK, FOS')
        
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 23, output
        
#which of those are in the immune pathways? 
class TestFindCommonPathwayGene3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, MEK1, SRF, HRAS, KRAS, JAK2, JAK1')
        keyword = 'immune'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways')) = ", len(output.get('pathways')))
        assert len(output.get('pathways')) == 4, output
        
#which of those have common pathways? 
class TestFindCommonPathwayGene4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, MEK1, SRF, HRAS, KRAS, JAK2, JAK1')
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways')) = ", len(output.get('pathways')))
        assert len(output.get('pathways')) == 30, output
        
#Which immune pathways are shared by STAT3, SOCS3, and CREB5 genes?
#(subtask: find-common-pathway-genes-keyword)
class TestFindCommonPathwayGene5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, SOCS3, CREB5')
        keyword = 'immune'
       
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 2, output
        
##Which immune pathways are shared by AKT, STAT3, SOCS3, and CREB5 genes?
#(subtask: find-common-pathway-genes-keyword)
class TestFindCommonPathwayGene6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('AKT, STAT3, SOCS3, CREB5')
        keyword = 'immune'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 2, output
        
#What KEGG pathways involve ERBB2, JUN, and MAPK8?
#sub-task:find-common-pathway-genes-db
class TestFindCommonPathwayGene7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonPathwayGene7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('ERBB2, JUN, MAPK8')
        db = 'KEGG'
        
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('pathways'))=", str(len(output.get('pathways'))))
        assert len(output.get('pathways')) == 4, output

###################################################################################
# FIND-PATHWAY-DB-KEYWORD
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
class TestFindPathwayDbKeyword2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayDbKeyword2, self).__init__(TFTA_Module)
        
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
class TestFindPathwayDbKeyword3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayDbKeyword3, self).__init__(TFTA_Module)
        
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

################################################################################
# FIND-PATHWAY-KEYWORD
#What pathways involve calcium? (subtask: find-pathway-keyword)
class TestFindPathwayKeyword1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        content = KQMLList('FIND-PATHWAY-KEYWORD')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 9, output
        
#What pathways involve immune system?        
class TestFindPathwayKeyword2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathwayKeyword2, self).__init__(TFTA_Module)
        
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

###################################################################################

#test no-capability message
class TestNoCapability(_IntegrationTest):
    def __init__(self, *args):
        super(TestNoCapability, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = agent_clj_from_text('calcium regulated pathways')
        
        content = KQMLList('FIND-T-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason').head() == 'NO-CAPABILITY', output
        assert output.get('reason').get('name') == 'FIND-T-PATHWAY', output



