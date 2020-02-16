# TFTA
===========================
TFTA (Transcription Factor and Target Agent): The TFTA's task is to search for targets for known transcription factors, and TFs for known targets; TF and target involved in specific pathways and tissues; TF and target involved in specific GO terms;microRNA and target relationship; find evidence for given regulation from both our database and literatures. I also provides enrichment analysis, such as GO enrichment, pathway enrichment, and disease enrichment.

TFTA is designed to be launched in a commucation system, which ability is limited when running by itself. Welcome to visit our [biodialog system](http://54.84.114.146/).

KQML messaging classes used by the agent are available [here](https://github.com/bgyori/pykqml).

INDRA is available [here](https://github.com/sorgerlab/indra).

# Installing TFTA
========================
Note that currently the TFTA has limited usage on its own. It's
meant to be launched in the context of a communication system. 

`pip install git+https://github.com/xzhang2016/tfagent.git`

# Note
=======================
The master was merged with dev branch and now it's using JSON as input and output format. (merged on Nov. 1, 2019)
