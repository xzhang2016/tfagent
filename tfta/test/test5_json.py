from kqml import KQMLList, KQMLString
from tfta.tfta import TFTA
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


###################################################
# GO-ANNOTATION
# GO-ENRICHMENT
#
#
###################################################

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
    

###############################################################################
# GO-ANNOTATION
#Find the GO annotation for genes which have growth in their name.
class TestGOAnnotation1(_IntegrationTest):
    def __init__(self, *args):
        super(TestGOAnnotation1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        keyword = 'growth'
        content = KQMLList('GO-ANNOTATION')
        content.set('keyword', keyword)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("len(output.get('genes'))=", len(output.get('genes')))
        get_genes(output.get('genes'))
        assert len(output.get('genes')) == 171, output


#############################################################################
# GO-ENRICHMENT
#find enriched GO terms for genes which have growth in their name (returned from previous query).
class TestGOEnrichment1(_IntegrationTest):
    def __init__(self, *args):
        super(TestGOEnrichment1, self).__init__(TFTA_Module)
        
    def create_message(self):
        # Here we create a KQML request that the TFTA needs to respond to
        genes = 'FAM168B,MEGF11,TGFBR3L,TGFB1,IGF2BP3,MSTN,HGS,FGF10,TGFB1I1,FGF16,FIBP,FRS3,GAS2,VEGFD,FGF17,GPRIN2,GDF9,GAS7,MEGF6,GADD45B,FGF18,GADD45G,GDF11,FGF19,EGFR,PDGFB,EGF,TGFA,NGF,GH2,IGF2,PDGFA,NTRK1,IGF1,FGF1,IGF1R,MET,FGF4,IGFBP1,FGF2,CXCL1,GIG44,PDGFRB,TGFB3,KIT,FGF6,GHR,FGFR1,FGF3,FGF5,TDGF1,HGF,VEGFA,PDGFRA,IGFBP3,FLT1,IGFBP2,EGR1,FGF7,FGFR2,FGFR4,FGFR3,IGFBP4,GADD45A,IGFBP6,IGFBP5,MST1,GDF1,FGF9,IGFALS,FLT4,KDR,EPS15,GDF5,PGF,VEGFB,VEGFC,HDGF,TDGF1P3,GAS1,FGF8,GDF10,FGF12,TGFB2,GRB2,GHRHR,TGFBR3,HGFAC,EGR4,EGR3,EPS8,GRB10,LRRC32,GAS6,GRB14,GRB7,FGFBP1,LTBP1,LTBP2,PDGFRL,TGFBI,IGFBP7,NTRK3,NTRK2,MIA,TBRG1,GRTP1,OGFRL1,IGFL4,GDF6,IGFL1,IGFL2,IGFL3,GPRIN1,NEGR1,GDF7,HDGFL2,MEGF8,TMEM219,NRROS,EGFL6,LTBP4,FGFRL1,DNER,GADD45GIP1,FGFBP3,EPS8L3,EPS8L1,FRS2,TGFBRAP1,IGFBPL1,ING5,GHSR,FGF13,FGF11,FGF14,MYDGF,RERG,MEGF10,HBEGF,CGREF1,CGRRF1,EGFL8,GDF15,FGFBP2,PDGFD,FGF23,ING2,MEGF9,GHITM,RERGL,EPS8L2,FGF22,FGF20,GDF3,PDGFC,LTBP3,FGF21,LYAR,ING3,IGF2BP1,OGFR,EPS15L1,EGFL7,OSGIN1,GDF2,ING1,ING4,OSGIN2,HDGFL3,IGF2BP2'
        clj = make_json(genes)
        content = KQMLList('GO-ENRICHMENT')
        content.set('gene', clj)
        return get_request(content), content
        
    def check_response_to_message(self, output):
        assert output.head() == 'SUCCESS', output
        print("\nlen(output.get('results'))=", len(output.get('results')))
        assert len(output.get('results')) == 14, output





