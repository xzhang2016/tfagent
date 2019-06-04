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
# IS-GENE-ONTO
# FIND-GENE-ONTO
# FIND-KINASE-REGULATION
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
#IS-GENE-ONTO
#Is stat3 a kinase?    
class TestIsGeneOnto1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('STAT3')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output

#Is stat3 a transcription factor?        
class TestIsGeneOnto2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('STAT3')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'transcription factor'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
#Is stat3 a protein kinase?
class TestIsGeneOnto3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('STAT3')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'protein kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output

#Is jak1 a protein kinase?
class TestIsGeneOnto4(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('JAK1')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'protein kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

#Is PBRM1 a transcription factor?
class TestIsGeneOnto5(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('PBRM1')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'transcription factor'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

#TEST FAMILY NAME
#Is SMURF a transcription factor?
class TestIsGeneOnto6(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('SMURF')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'transcription factor'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output.get('reason')
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#is stat a kinase? (STAT is grounded as a gene, not a family)
class TestIsGeneOnto7(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('STAT')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'KINASE'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#test protein and gene category
#is map3k7 a protein?
class TestIsGeneOnto8(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('map3k7')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'protein'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

#is stat3 a gene?
class TestIsGeneOnto9(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto9, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('stat3')
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'gene'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

##################################################################################
##TEST FIND-GENE-ONTO
#Among STAT3, JAK1, JAK2, ELK1, ELK2, HRAS, and FOS, which are protein kinases?
class TestFindGeneOnto1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text("STAT3, JAK1, JAK2, ELK1, FOS, HRAS, ELK2")
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'protein kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print("len(output.get('genes'))=", str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output

#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are histone demethylase?
class TestFindGeneOnto2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text("STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B")
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'histone demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output

#Among PBRM1, SMAD2, TBL1XR1, AKT1, CDK19, CDK8, CDK9, DDR1, GSK3A, GSK3B, MET,TRIM28,COL2A1,
# JAK1, PRMT1, RB1, SMURF2, TRAF4, and USP15, which are transcription factors?
class TestFindGeneOnto3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text("PBRM1, SMAD2, TBL1XR1, AKT1, CDK19, CDK8, CDK9, DDR1, \
                GSK3A, GSK3B, MET,TRIM28,COL2A1,JAK1, PRMT1, RB1, SMURF2, TRAF4, USP15")
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'transcription factor'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print("len(output.get('genes'))=" + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 3, output

#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are demethylase?
class TestFindGeneOnto4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text("STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B")
        _get_targets(gene)
        print('target=', str(gene))
        
        keyword = 'demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output

#complex query: find-target and find-gene-onto
#What genes regulated by FOS are kinases?
class TestFindGeneOnto5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('cfos')
        _get_targets(regulator)
        print('target=', str(regulator))
        
        keyword = 'kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('regulator', regulator)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 5, output

###############################################################################
# FIND-KINASE-REGULATION


