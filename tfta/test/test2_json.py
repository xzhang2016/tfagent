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
# IS-GENE-TISSUE
# FIND-GENE-TISSUE
# FIND-TISSUE
# FIND-COMMON-TF-GENES
# FIND-EVIDENCE
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
#Which kinases regulate the cfos gene?
class TestFindKinaseReg1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('cfos')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('kinase')) == 5, output

#test gene family
#Which kinases regulate the MEK gene?
class TestFindKinaseReg2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('MEK')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#Which kinases negatively regulate the cfos gene?
class TestFindKinaseReg3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('cfos')
        _get_targets(target)
        print('target=', str(target))
        
        keyword = 'decrease'
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert len(output.get('kinase')) == 3, output
        
#Which kinases positively regulate the cfos gene?
class TestFindKinaseReg4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('cfos')
        _get_targets(target)
        print('target=', str(target))
        
        keyword = 'increase'
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'SUCCESS', output
        assert len(output.get('kinase')) == 2, output

#Which kinases positively regulate the AKT gene?
class TestFindKinaseReg5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindKinaseReg5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('AKT')
        _get_targets(target)
        print('target=', str(target))
        
        keyword = 'increase'
        content = KQMLList('FIND-KINASE-REGULATION')
        content.set('target', target)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output

#######################################################################################
#IS-GENE-TISSUE
###Is stat3 expressed in liver? 
class TestIsTissueGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('stat3')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', target)
        content.set('tissue', 'liver')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
###Is kras expressed in brain? 
class TestIsTissueGene2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('kras')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', target)
        content.set('tissue', 'brain')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output

###Is stat3 exclusively expressed in liver? 
class TestIsTissueGene3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('stat3')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', target)
        content.set('tissue', 'liver')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'FALSE', output
        
###Is GYS2 exclusively expressed in liver? 
class TestIsTissueGene4(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('GYS2')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', target)
        content.set('tissue', 'liver')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output
        
###Is NEUROD2 exclusively expressed in brain? 
class TestIsTissueGene5(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('NEUROD2')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', target)
        content.set('tissue', 'brain')
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

###Is GAST expressed in stomach? 
class TestIsTissueGene6(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsTissueGene6, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agent_clj_from_text('GAST')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('IS-GENE-TISSUE')
        content.set('gene', target)
        content.set('tissue', 'stomach')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE', output

######################################################################################
#FIND-GENE-TISSUE
#what genes are expressed in liver? 
class TestFindGeneTissue1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = ekb_kstring_from_text('AKT')
        tissue = 'liver'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 1929, output
        
#among stat3,srf, kras, and hras, what genes are expressed in liver? 
class TestFindGeneTissue2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agents_clj_from_text('stat3, srf, kras, hras')
        tissue = 'liver'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 1, output
        
#what genes are exclusively expressed in liver? 
class TestFindGeneTissue3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = 'stat3, srf, kras, hras'
        tissue = 'liver'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 31, output
        
#what genes are exclusively expressed in brain? 
class TestFindGeneTissue4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneTissue4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tissue = 'brain'
        content = KQMLList('FIND-GENE-TISSUE')
        content.set('tissue', tissue)
        content.set('keyword', 'exclusive')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        assert len(output.get('genes')) == 44, output

###############################################################################
# FIND-TISSUE
#What tissues is STAT3 expressed in?
class TestFindTissue1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('STAT3')
        _get_targets(gene)
        print('target=', str(gene))
        
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tissue'))=", str(len(output.get('tissue'))))
        assert len(output.get('tissue')) == 8, output

#What tissues is MEK expressed in?
class TestFindTissue2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('MEK')
        _get_targets(gene)
        print('target=', str(gene))
        
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        print("len(output.get('clarification'))=", len(output.get('clarification')))
        assert len(output.get('clarification')) == 5, output
        
#what tissues can I ask 
class TestFindTissue3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        content = KQMLList('FIND-TISSUE')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tissue'))=", len(output.get('tissue')))
        assert len(output.get('tissue')) == 30, output
        
#What tissues is frizzled8 expressed in?
class TestFindTissue4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = agent_clj_from_text('frizzled8')
        _get_targets(gene)
        print('target=', str(gene))
        
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tissue'))=", str(len(output.get('tissue'))))
        assert len(output.get('tissue')) == 7, output

####################################################################################
#FIND-COMMON-TF-GENES
#What transcription factors are shared by the SRF, HRAS, and elk1 genes? (subtask: find-common-tf-genes)
class TestFindCommonTfGenes1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('SRF, HRAS, elk1')
        _get_targets(target)
        print('target=', str(target))
      
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output

#What transcription factors are in common to the STAT3, SOCS3, IFNG, FOXO3, and CREB5 genes?
class TestFindCommonTfGenes2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, IFNG, FOXO3, SOCS3, CREB5')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('tfs'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 8, output

#test gene family
#What transcription factors are in common to the STAT3, SOCS3, and MEK genes?
#MEK will be ignored in this case
class TestFindCommonTfGenes3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, SOCS3, MEK')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 1, output

#What transcription factors are in common to the STAT3, SOCS3, and AKT genes?
class TestFindCommonTfGenes4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('STAT3, AKT, MEK')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'FAMILY_NAME', output
        assert len(output.get('clarification').get('as')) == 2, output

#Which of these transcription factors are shared by the SRF, HRAS, FOS, and elk1 genes? (subtask: find-common-tf-genes)
class TestFindCommonTfGenes5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindCommonTfGenes5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = agents_clj_from_text('SRF, HRAS, cfos, elk1')
        _get_targets(target)
        print('target=', str(target))
        
        of_those = agents_clj_from_text('stat3,ELK1,TFAP2A,CREB1,TP53')
        _get_targets(of_those)
        print('target=', str(of_those))
        
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        content.set('of-those', of_those)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('tfs'))=", len(output.get('tfs')))
        assert len(output.get('tfs')) == 3, output

######################################################################################
# FIND-EVIDENCE
##Show me evidence that kras regulate frizzled8? 
class TestFindEvidence1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('kras')
        target = agent_clj_from_text('fzd8')
        _get_targets(target)
        print('target=', str(target))
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', regulator)
        content.set('target', target)
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("len(output.get('evidence').get('literature'))=", str(len(output.get('evidence').get('literature'))))
        assert len(output.get('evidence')) == 2, output
        assert len(output.get('evidence').get('literature')) == 1, output
        
##show me evidence that kras increase frizzled8? 
class TestFindEvidence2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence2, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('kras')
        target = agent_clj_from_text('fzd8')
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', regulator)
        content.set('target', target)
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("len(output.get('evidence').get('literature'))=", str(len(output.get('evidence').get('literature'))))
        assert len(output.get('evidence')) == 2, output
        assert len(output.get('evidence').get('literature')) == 1, output
        
##show me evidence that kras decrease frizzled8? 
class TestFindEvidence3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('kras')
        target = agent_clj_from_text('fzd8')
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', regulator)
        content.set('target', target)
        content.set('keyword', 'decrease')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("type(output.get('evidence'))=",type(output.get('evidence')))
        print("len(output.get('evidence').get('literature'))=", str(len(output.get('evidence').get('literature'))))
        assert len(output.get('evidence')) == 2, output
        assert len(output.get('evidence').get('literature')) == 1, output
        
##Show me the evidence that IL6 increase the amount of SOCS1.
class TestFindEvidence4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence4, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('il6')
        target = agent_clj_from_text('socs1')
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator',regulator)
        content.set('target',target)
        content.set('keyword', 'increase')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("len(output.get('evidence').get('literature'))=", str(len(output.get('evidence').get('literature'))))
        assert len(output.get('evidence')) == 2, output
        assert len(output.get('evidence').get('literature')) == 9, output
        
##Show me the evidence that SRF binds to the FOS gene.
class TestFindEvidence5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence5, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('SRF')
        target = agent_clj_from_text('cfos')
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', regulator)
        content.set('target', target)
        content.set('keyword', 'bind')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("len(output.get('evidence').get('tf-db'))=", str(len(output.get('evidence').get('tf-db'))))
        assert len(output.get('evidence')) == 2, output
        assert len(output.get('evidence').get('tf-db')) == 2, output

##Show me the evidence that SRF regulate FOS gene.
class TestFindEvidence6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence6, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        regulator = agent_clj_from_text('SRF')
        target = agent_clj_from_text('cfos')
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', regulator)
        content.set('target', target)
        content.set('keyword', 'regulate')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('evidence'))=", str(len(output.get('evidence'))))
        print("len(output.get('evidence').get('literature'))=", str(len(output.get('evidence').get('literature'))))
        assert len(output.get('evidence')) == 4, output
        assert len(output.get('evidence').get('literature')) == 2, output

#IncreaseAmount(miR_491(), GFAP())
class TestFindEvidence7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindEvidence7, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        agents = Bioagent.get_agent(agent_clj_from_text('miR-491'))
        print(agents)
        print('name=', agents.name)
        print('db_refs=', agents.db_refs)
        regulator = agent_clj_from_text('miR-491')
        target = agent_clj_from_text('GFAP')
        
        content = KQMLList('FIND-EVIDENCE')
        content.set('regulator', regulator)
        content.set('target', target)
        content.set('keyword', 'increase')
        print(content, '\n')
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_REGULATOR_NAME', output





