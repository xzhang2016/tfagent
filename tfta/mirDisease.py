import os
import time
#import wget
#import pandas as pd
import logging
import sqlite3
import pickle
from collections import defaultdict


logging.basicConfig(format='%(levelname)s: %(name)s - %(message)s',
                    level=logging.INFO)
logger = logging.getLogger('TFTA-mirDisease')

_resource_dir = os.path.dirname(os.path.realpath(__file__)) + '/../resources/'

class mirDisease:
    def __init__(self):
        #load db file
        self.mirdb = _load_db()
        self.mirna_precursor = _load_mirna_precursor_mapping()
        
    def __del__(self):
        self.mirdb.close()
    
    def is_mirna_disease(self, mirna_name_dict, disease):
        """
        If the mirna associated with the disease, return ture; otherwise return false.
        
        parameter
        -----------
        mirna_name_dict: dict
        disease: str
        """
        mirna_name = list(mirna_name_dict.keys())[0]
        mapped_mirna = self._to_precursor([mirna_name])
        if self.mirdb is not None:
            dstr = '%' + disease + '%'
            for prec in mapped_mirna[mirna_name]:
                t = (prec, dstr)
                res = self.mirdb.execute("SELECT mirna FROM mir2disease WHERE mirna LIKE ? "
                                    "AND disease LIKE ?", t).fetchone()
                if res:
                    return True
        return False
        
    def find_disease_mirna(self, mirna_names_dict):
        """
        For a given list of mirnas, return the associated diseases.
        
        paramter
        -------------
        mirna_names_dict: dict
        """
        mirna_names = list(mirna_names_dict.keys())
        mapped_mirnas = self._to_precursor(mirna_names)
        #id = set()
        disease = defaultdict(set)
        fdisease = set()
        if self.mirdb is not None:
            for mir in mirna_names:
                associated = False
                for prec in mapped_mirnas[mir]: 
                    t = (prec,)
                    res = self.mirdb.execute("SELECT disease FROM mir2disease WHERE mirna LIKE ?", t).fetchall()
                    if res:
                        #id = set([r[0] for r in res])
                        disease[mir] = disease[mir].union(set([r[0] for r in res]))
                        associated = True
                if not associated:
                    return fdisease
            
            #get the intersection if each mirna has associated diseases
            fdisease = disease[mirna_names[0]]
            for m in mirna_names:
                fdisease = fdisease.intersection(disease[m])
        return fdisease
        
    def find_mirna_disease(self, disease):
        """
        Return the mirnas associated with the given disease.
        
        parameter
        -------------
        disease: str
        """
        mirnas = []
        if self.mirdb is not None:
            dstr = '%' + disease + '%'
            t = (dstr,)
            res = self.mirdb.execute("SELECT DISTINCT mirna FROM mir2disease WHERE disease LIKE ?", t).fetchall()
            if res:
                mirnas = [r[0] for r in res]
        return mirnas
        
    def find_disease(self):
        """
        Return all the diseases in the db file.
        """
        disease = set()
        if self.mirdb is not None:
            res = self.mirdb.execute("SELECT DISTINCT disease FROM mir2disease").fetchall()
            disease = set([r[0] for r in res])
        return disease
        
    def find_common_disease_mirna(self, mirna_dict):
        """
        Return shared diseases by some of the mirnas.
        
        parameter
        --------------
        mirna_dict: dict
        """
        mirna_names = list(mirna_dict.keys())
        disease = defaultdict(set)
        if self.mirdb is not None:
            pass
        
        
                
    def _to_precursor(self, mirna):
        """
        Mapping each matured mirna to its precursors.
    
        parameter
        ----------
        mirna: list
        """
        mir2pre = defaultdict(list)
        for mir in mirna:
            if mir.endswith('-3P') or mir.endswith('-5P'):
                try:
                    mir2pre[mir] = self.mirna_precursor[mir]
                except KeyError:
                    mir2pre[mir] = None
            else:
                mir2pre[mir] = [mir]
        return mir2pre

def _load_db():
    db_file = os.path.join(_resource_dir, 'mirnaDisease.db')
        
    #check file size to determine if it need regenerate
    if os.path.isfile(db_file):
        statinfo = os.stat(db_file)
    else:
        statinfo = None
        
    if statinfo and statinfo.st_size > 0:
        mirDisease = sqlite3.connect(db_file, check_same_thread=False)
        logger.info('TFTA loaded mirnaDisease database.')
    else:
        logger.info('Generating the mirnaDisease db file...')
        num = _generate_mirDisease_db(db_file)
        if os.path.isfile(db_file):
            mirDisease = sqlite3.connect(db_file, check_same_thread=False)
        else:
            logger.error('TFTA could not load mirnaDisease database.')
            mirDisease = None;
    return mirDisease
    
def _load_mirna_precursor_mapping():
    fn = os.path.join(_resource_dir, "mirna_precursor_dict.pickle")
    if os.path.isfile(fn):
        with open(fn, "rb") as pickle_in:
            mirna_precursor = pickle.load(pickle_in)
    else:
        mirna_precursor = _map_mirna_precursor(fn)
    logger.info('TFTA loaded mirna_precursor mapping file.')
    return mirna_precursor

def _generate_mirDisease_db(db_file, data_folder=_resource_dir):
    """
    Generate a sqlite db file which contains MicroRNAs and the associated diseases.
    """
    #read data
    with open(os.path.join(data_folder, 'mirna_disease.txt'), 'r') as fr:
        fr.readline()
        lines = fr.readlines()
    
    #generate table
    t0 = time.perf_counter()
        
    conn = sqlite3.connect(db_file)
    c = conn.cursor()
    c.execute('''CREATE TABLE mir2disease
        (id integer, mirna text, disease text, pmid text)''')
    
    records = []
    num = 1
    for line in lines:
        record = _parse_line(line, num)
        if record:
            records.append(record)
            num += 1

    c.executemany('INSERT INTO mir2disease VALUES (?,?,?,?)', records)
    conn.commit()
    time.sleep(0.1)
    c.close()
    conn.close()
        
    t1 = time.perf_counter()
    logger.info('Created {}. The mir2disease table has {} items.'.format(db_file.split('/')[-1], num-1))
    logger.info('Used {} seconds to generate {} db file.'.format(t1-t0, db_file.split('/')[-1]))
    return num-1
    
def _parse_line(line, ind):
    t = []
    s = line.strip().split('\t')
    t.append(ind)
    t.append(s[0])
        
    #process disease field
    disease = s[1]
    #if disease.startswith('"') and disease.endswith('"'):
    #    disease = disease[1:-1]
    if ',' in disease:
        disease = disease.split(', ')
        disease.reverse()
        disease = ' '.join(disease)
    t.append(disease)
    t.append(s[2])

    return tuple(t)

def _map_mirna_precursor(fn_pickle, data_folder=_resource_dir):
    """
    Generate the matured mirna to precursor mapping file since the disease db contains the precursors
    instead of matured mirnas.
    Use miRNA.dat at http://www.mirbase.org/ftp.shtml. (2019-09-04)
    
    parameter
    -----------
    fn_pickle: str
    data_folder: str
    """
    mirna_precursor = defaultdict(list)
    with open(os.path.join(data_folder, 'hsa_precursor_mature_mirna.txt'), 'r') as fr:
        #skip title line
        fr.readline()
        lines = fr.readlines()
        
    for line in lines:
        s = line.strip().split('\t')
        mirna_precursor[s[0].upper()] = s[1].upper().split(',')
        
    with open(fn_pickle, 'wb') as pickle_out:
        pickle.dump(mirna_precursor, pickle_out)
    
    return mirna_precursor
    

