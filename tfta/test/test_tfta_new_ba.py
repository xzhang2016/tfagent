from kqml import KQMLList
from tfta.tfta import TFTA
from tfta.tfta_module import TFTA_Module
from bioagents.tests.util import ekb_from_text, ekb_kstring_from_text, \
                                 get_request
from bioagents.tests.integration import _IntegrationTest, _FailureTest

#TEST IS-REGULATION
##Does STAT3 regulate the c-fos gene? (subtask: is-tf-target)
class TestIsRegulation1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#Does STAT3 regulate the JUN gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation2(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #since trips parser cannot recognize JUN as a gene, so instead try
        #c-JUN or JUN gene
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-JUN')
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
class TestIsRegulation21(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        tissue = 'liver'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'FALSE', output
        
#Does STAT3 regulate the cfos gene in the lung? (subtask: is-tf-target-tissue)
class TestIsRegulation22(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        target = ekb_kstring_from_text('c-fos')
        tissue = 'lung'
        content = KQMLList('IS-REGULATION')
        content.set('tf', tf)
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.gets('result') == 'TRUE', output
        
#Does miR-20b-5p target STAT3? (subtask: is-mirna-target)
class TestIsRegulation3(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsRegulation3, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3')
        mirna = ekb_kstring_from_text('miR-20b-5p')
        content = KQMLList('IS-MIRNA-TARGET')
        content.set('target', target)
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('is-miRNA-target') == 'TRUE', output
        
#TEST FIND-TF
#Which transcription factors regulate frizzled8? (subtask: find-target-tf)
class TestFindTf1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf1, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('frizzled8')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
    
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 2, output
        
##What regulates SMURF2?
#What are the regulators of SMURF2? (subtask: find-target-tf)
class TestFindTf11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf11, self).__init__(TFTA_Module)

    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SMURF2')
        content = KQMLList('FIND-TF')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3
        
#test sequencing query
#Which of these genes are regulated by STAT3? (subtask: find-tf-target)
class TestFindTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        of_targets = "A2M, ABCA2, AKAP12, AKT1, PBRM1, SMAD2,CEBPA"
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('of-targets', of_targets)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 5, output
        
#Which transcription factors regulate frizzled8 in liver? (subtask: find-target-tf-tissue)
class TestFindTf2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('frizzled8')
        tissue = 'liver'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'TARGET_NOT_FOUND', output
        
#Which transcription factors regulate mapk14 in bladder? (subtask: find-target-tf-tissue)
class TestFindTf21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('mapk14')
        tissue = 'bladder'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 5, output
        
#Which transcription factors regulate ELK1 in brain? (subtask: find-target-tf-tissue)
class TestFindTf22(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf22, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('ELK1')
        tissue = 'brain'
        content = KQMLList('FIND-TF')
        content.set('target', target)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 55, output
        
#What transcription factors are in the calcium regulated pathways? (subtask: find-tf-keyword)
class TestFindTf3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        content = KQMLList('FIND-TF')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 6, output
        
#What transcription factors are shared by the SRF, HRAS, and elk1 genes? (subtask: find-common-tf-genes)
class TestFindTf4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('SRF, HRAS, elk1')
        #count = 'count'
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output
        
#What transcription factors are in common to the STAT3, SOCS3, and CREB5 genes?
class TestFindTf41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf41, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, CREB5')
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 3, output

#What transcription factors are in common to the STAT3, SOCS3, IFNG, FOXO3, and CREB5 genes?        
class TestFindTf42(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf42, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, IFNG, FOXO3, SOCS3, CREB5')
        content = KQMLList('FIND-COMMON-TF-GENES')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('tfs'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tfs')) == 8, output
        
#Which transcription factors are in the MAPK signaling pathway? (subtask: find-tf-pathway)
class TestFindTf5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTf5, self).__init__(TFTA_Module)
    
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        content = KQMLList('FIND-TF')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 17, output
        
#TEST FIND-TARGET
#what genes are regulated by smad2? (subtask: find-tf-target)
class TestFindTarget1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('smad2')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 115, output
        
#what genes are regulated by elk1 and smad2? (subtask: find-tf-target)
class TestFindTarget11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('ELK1, SMAD2')
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 28
        
#test sequencing query
#Which of these genes are regulated by STAT3? (subtask: find-tf-target)
class TestFindTarget12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        of_targets = "A2M, ABCA2, AKAP12, AKT1, PBRM1, SMAD2,CEBPA"
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('of-targets', of_targets)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 5, output
        
#What genes does stat3 regulate in lung? (subtask: find-tf-target-tissue)
class TestFindTarget2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        tissue = 'lung'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 57, output
        
#hat genes does stat3 regulate in liver? (subtask: find-tf-target-tissue)
class TestFindTarget21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        tf = ekb_kstring_from_text('STAT3')
        tissue = 'liver'
        content = KQMLList('FIND-TARGET')
        content.set('tf', tf)
        content.set('tissue', tissue)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'NO_TARGET_FOUND', output
        
#What genes does miR-20b-5p target? (subtask: find-target-mirna)
class TestFindTarget3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-20b-5p')
        content = KQMLList('FIND-TARGET-MIRNA')
        content.set('miRNA', mirna)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 917, output
        
#What genes are most frequently regulated by miR-335-5p, miR-155-5p, miR-145-5p, and miR-20a-5p?
#(subtask: FIND-GENE-COUNT-MIRNA)
class TestFindTarget4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-335-5p, miR-155-5p, miR-145-5p, miR-20a-5p')
        count = 'count'
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 380, output
        
#What genes are most frequently regulated by miR-335-5p, miR-155-5p and miR-145-5p?
#(subtask: FIND-GENE-COUNT-MIRNA)
class TestFindTarget41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTarget41, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        mirna = ekb_kstring_from_text('miR-335-5p, miR-155-5p, miR-145-5p')
        count = 'count'
        content = KQMLList('FIND-GENE-COUNT-MIRNA')
        content.set('miRNA', mirna)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('targets')) == 151, output
        
#TEST FIND-GENE
#What genes are in the MAPK signaling pathway? (subtask: find-gene-pathway)
class TestFindGene1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGene1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('MAPK')
        content = KQMLList('FIND-GENE')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 23, output
        
#What genes are involved in the il-12 pathway?
class TestFindGene11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGene11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('il-12')
        content = KQMLList('FIND-GENE')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 3, output
        
#What genes are in the immune system pathway?
class TestFindGene12(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGene12, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        pathway = ekb_kstring_from_text('immune system')
        content = KQMLList('FIND-GENE')
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#not implemented for find-gene-go-tf

#TEST FIND-PATHWAY
#What pathways involve SRF? (subtask: find-pathway-gene)
class TestFindPathway1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SRF')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 27, output
        
#What pathways involve kras and elk1? (subtask: find-pathway-gene)
class TestFindPathway11(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway11, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('KRAS, ELK1')
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 23, output
        
#Which Reactome pathways utilize SRF? (subtask: find-pathway-db-gene)
class TestFindPathway2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SRF')
        db = 'Reactome'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#Which KEGG pathways utilize SRF? (subtask: find-pathway-db-gene)
class TestFindPathway21(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway21, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SRF')
        db = 'KEGG'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('database', db)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#What pathways involve calcium? (subtask: find-pathway-keyword)
class TestFindPathway3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'calcium'
        content = KQMLList('FIND-PATHWAY')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 10, output
        
#Which pathways are shared by STAT3, SOCS3, IFNG, FOXO3, and CREB5 genes? 
#(subtask: find-common-pathway-genes)
class TestFindPathway4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, IFNG, FOXO3, CREB5')
        #count = 'count'
        content = KQMLList('find-common-pathway-genes')
        content.set('target', target)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 8, output
        
#What signaling pathways are shared by STAT3 and SRF? (subtask: find-common-pathway-genes)
class TestFindPathway41(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway41, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3, SRF')
        count = 'count'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 6, output
        
#What immune pathways involve kras and elk1? (subtask: find-pathway-gene-keyword)
class TestFindPathway5(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway5, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('kras, elk1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#What immune pathways involve tap1 and jak1? (subtask: find-pathway-gene-keyword)
class TestFindPathway51(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway51, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('tap1, jak1')
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('gene', gene)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 1, output
        
#Which immune pathways are shared by STAT3, SOCS3, and CREB5 genes?
#(subtask: find-common-pathway-genes-keyword)
class TestFindPathway6(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway6, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3, SOCS3, CREB5')
        keyword = 'immune'
        #count = 'count'
        content = KQMLList('FIND-COMMON-PATHWAY-GENES')
        content.set('target', target)
        content.set('keyword', keyword)
        #content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('pathways'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#Does the mTor pathway utilize SGK1? (subtask: is-pathway-gene)
class TestFindPathway7(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway7, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('SGK1')
        pathway = ekb_kstring_from_text('mTor pathway')
        content = KQMLList('IS-PATHWAY-GENE')
        content.set('gene', gene)
        content.set('pathway', pathway)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 3, output
        
#What KEGG pathways involve immune signaling? (subtask: find-pathway-db-keyword)
class TestFindPathway8(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway8, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        database = 'KEGG'
        keyword = 'immune'
        content = KQMLList('FIND-PATHWAY')
        content.set('database', database)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 2, output
        
#What KEGG pathways involve immune system? (subtask: find-pathway-db-keyword)
class TestFindPathway81(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway81, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        database = 'KEGG'
        keyword = 'immune system'
        content = KQMLList('FIND-PATHWAY')
        content.set('database', database)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'FAILURE', output
        assert output.get('reason') == 'PathwayNotFoundException', output
        
#What reactome pathways involve immune system? (subtask: find-pathway-db-keyword)
class TestFindPathway82(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindPathway82, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        database = 'reactome'
        keyword = 'immune system'
        content = KQMLList('FIND-PATHWAY')
        content.set('database', database)
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('pathways')) == 4, output
        
#TEST FIND-MIRNA
#What microRNAs target STAT3? (subtask: FIND-MIRNA-TARGET)
class TestFindMirna1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('STAT3')
        content = KQMLList('FIND-MIRNA')
        content.set('target', target)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('miRNAs')) == 80, output
        
#What miRNAs most frequently regulate EGFR, SRF, STAT3, JAK2, and SMAD3?
#(subtask: FIND-MIRNA-COUNT-GENE)
class TestFindMirna2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindMirna2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        target = ekb_kstring_from_text('EGFR, SRF, STAT3, JAK2, SMAD3')
        count = 'count'
        content = KQMLList('FIND-MIRNA')
        content.set('target', target)
        content.set('count', count)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('miRNAs')) == 23, output
        
#find tissue
#What tissues is STAT3 expressed in?
class TestFindTissue1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindTissue1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
        content = KQMLList('FIND-TISSUE')
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert len(output.get('tissue')) == 29, output

#Is stat3 a kinase?    
class TestIsGeneOnto1(_IntegrationTest):
    def __init__(self, *args):
        super(TestIsGeneOnto1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3')
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
        gene = ekb_kstring_from_text('STAT3')
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
        gene = ekb_kstring_from_text('STAT3')
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
        gene = ekb_kstring_from_text('JAK1')
        keyword = 'protein kinase'
        content = KQMLList('is-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        assert output.get('result') == 'TRUE, output

#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are protein kinases?
class TestFindGeneOnto1(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3, JAK1, JAK2, ELK1, FOS')
        keyword = 'protein kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are kinases?
class TestFindGeneOnto2(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto2, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3, JAK1, JAK2, ELK1, FOS')
        keyword = 'kinase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 2, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are histone demethylase?
class TestFindGeneOnto3(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto3, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        gene = ekb_kstring_from_text('STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B')
        keyword = 'histone demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output
        
#Among STAT3, JAK1, JAK2, ELK1, and FOS, which are histone demethylase?
#test use gene list in class variable
class TestFindGeneOnto4(_IntegrationTest):
    def __init__(self, *args):
        super(TestFindGeneOnto4, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        #gene = ekb_kstring_from_text('STAT3, JAK1, JAK2, ELK1, FOS, SMAD2, KDM4B')
        keyword = 'histone demethylase'
        content = KQMLList('find-gene-onto')
        content.set('keyword', keyword)
        #content.set('gene', gene)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        print('len(output)=' + str(len(output.get('genes'))))
        assert output.head() == 'SUCCESS', output
        assert len(output.get('genes')) == 1, output

if __name__ == '__main__':
    TestIsRegulation1().run_test()

