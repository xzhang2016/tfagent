from kqml import KQMLList
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest


#FIND-TF-TARGET
class _TestFindTfTarget(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindTfTarget, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text(self.tf)
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


#what genes are regulated by smad2?
class TestFindTfTarget1(_TestFindTfTarget):
    tf = 'SMAD2'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 115, output

#What genes does stat3 regulate?


#what genes are regulated by elk1 and smad2?
class TestFindTfTarget2(_TestFindTfTarget):
    tf = 'ELK1, SMAD2'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 28

#test IS-TF-TARGET
class _TestIsTfTarget(_IntegrationTest):
    def __init__(self, *args):
        super(_TestIsTfTarget, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text(self.tf)
        target = ekb_kstring_from_text(self.target)
        content = KQMLList('IS-TF-TARGET')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content

#Does STAT3 regulate the c-fos gene?
class TestIsTfTarget1(_TestIsTfTarget):
    tf = 'STAT3'
    target = 'c-fos'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#TEST FIND-TARGET-TF
class _TestFindTargetTf(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindTargetTf, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text(self.target)
        content = KQMLList('FIND-TARGET-TF')
        content.set('target', target)
        return get_request(content), content

    def check_response_to_message(self, output):
        # Here we check the details of the response and make sure everything
        # we are looking for is in the output
        # First, we check that the response is a success
        assert output.head() == 'SUCCESS', output
        # Then we check that we got the expected number of target
        assert len(output.get('targets')) == 116
        # We could do further checks here to see if a given expected target
        # shows up in the response, etc.

#Which transcription factors regulate frizzled8?
class TestFindTargetTf1(_TestFindTargetTf):
    target = 'frizzled8'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 2, output

#What regulates SMURF2?
#What are the regulators of SMURF2?
class TestFindTfTarget2(_TestFindTfTarget):
    target = 'SMURF2'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3
    
#TEST FIND-PATHWAY-GENE
class _TestFindPathwayGene(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindPathwayGene, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text(self.gene)
        content = KQMLList('FIND-PATHWAY-GENE')
        content.set('gene', gene)
        return get_request(content), content
        
#What pathways involve SRF?
class TestFindPathwayGene1(_TestFindPathwayGene):
    gene = 'SRF'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 27, output
        
#What pathways involve kras and elk1?
class TestFindPathwayGene2(_TestFindPathwayGene):
    gene = 'KRAS, ELK1'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 23, output

#TEST FIND-PATHWAY-DB-GENE
class _TestFindPathwayDBGene(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindPathwayDBGene, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text(self.gene)
        db = self.database
        content = KQMLList('FIND-PATHWAY-DB-GENE')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
#Which Reactome pathways utilize SRF?
class TestFindPathwayDBGene1(_TestFindPathwayDBGene):
    gene = 'SRF'
    database = 'Reactome'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#Which KEGG pathways utilize SRF?
class TestFindPathwayDBGene1(_TestFindPathwayDBGene):
    gene = 'SRF'
    database = 'KEGG'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output

#TEST FIND-TF-PATHWAY
class _TestFindTfPathway(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindTfPathway, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text(self.pathway)
        content = KQMLList('FIND-TF-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
#Which transcription factors are in the MAPK signaling pathway?
class TestFindTfPathway1(_TestFindTfPathway):
    pathway = 'MAPK'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 17, output
        
#TEST FIND-GENE-PATHWAY
class _TestFindGenePathway(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindGenePathway, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text(self.pathway)
        content = KQMLList('FIND-GENE-PATHWAY')
        content.set('pathway', pathway)
        return get_request(content), content
        
#What genes are in the MAPK signaling pathway?
class TestFindGenePathway1(_TestFindGenePathway):
    pathway = 'MAPK'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 23, output
        
#What genes are involved in the il-12 pathway?
class TestFindGenePathway2(_TestFindGenePathway):
    pathway = 'il-12'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 3, output
        
#What genes are in the immune system pathway?
class TestFindGenePathway3(_TestFindGenePathway):
    pathway = 'immune system'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#TEST FIND-PATHWAY-KEYWORD
class _TestFindPathwayKeyword(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindPathwayKeyword, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = self.keyword
        content = KQMLList('FIND-PATHWAY-KEYWORD')
        content.set('keyword', keyword)
        return get_request(content), content

#What pathways involve calcium?
class TestFindPathwayKeyword1(_TestFindPathwayKeyword):
    keyword = 'calcium'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 10, output
        
#FIND-TF-KEYWORD
class _TestFindTfKeyword(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindTfKeyword, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = self.keyword
        content = KQMLList('FIND-TF-KEYWORD')
        content.set('keyword', keyword)
        return get_request(content), content

#What transcription factors are in the calcium regulated pathways?
class TestFindTfKeyword1(_TestFindTfKeyword):
    keyword = 'calcium'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#TEST FIND-COMMON-TF-GENES
class _TestFindCommonTfGene(_IntegrationTest):
    def __init__(self, *args):
        super(_TestFindCommonTfGene, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text(self.target)
        content = KQMLList('FIND-PATHWAY-GENE')
        content.set('target', target)
        return get_request(content), content
        
#What transcription factors are shared by the SRF, HRAS, and elk1 genes?
class TestFindCommmonTfGene1(_TestFindCommonTfGene):
    target = 'SRF, HRAS, elk1'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 54, output
        
#What transcription factors are in common to the STAT3, SOCS3, and CREB5 genes?
class TestFindCommonTfGene1(_TestFindCommonTfGene):
    target = 'STAT3, SOCS3, CREB5'
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 94, output
        
#
