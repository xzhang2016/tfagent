TF agent: TFTA
===========================
TFTA (Transcription Factor and Target Agent): The TFTA's task is to search for targets for known transcription factor, and TFs for known target; TF and target involved in specific pathways; TF and target involved in specific GO terms.

KQML messaging classes used by the agents are available at: https://github.com/bgyori/pykqml

INDRA is available at: https://github.com/sorgerlab/indra

Installing the TFTA
========================
Note that currently the TFTA has limited usage on its own. It's
meant to be launched in the context of a communication system. 

The bioagents depend on the following non-default python packages: objectpath,
rdflib, functools32, requests, lxml, pandas, suds

Please follow the more detailed instructions on the [INDRA page](https://github.com/sorgerlab/indra) 
to install it and its dependencies:

`pip install git+https://github.com/sorgerlab/indra.git`

INDRA depends on [PySB](http://pysb.org), which is best installed from Github:

`pip install git+https://github.com/pysb/pysb.git`

PySB depends on [BioNetGen](http://bionetgen.org/index.php/Download). Make sure
that BioNetGen is unzipped into /usr/local/share/BioNetGen, such that BNG2.pl is located at /usr/local/share/BioNetGen/BNG2.pl. Alternatively, set BNGPATH 
to the folder in which BNG2.pl is.


