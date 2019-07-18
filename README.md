TF agent: TFTA
===========================
TFTA (Transcription Factor and Target Agent): The TFTA's task is to search for targets for known transcription factor, and TFs for known target; TF and target involved in specific pathways; TF and target involved in specific GO terms;microRNA and target relationship; find evidence for given regulation from both our database and literature.

KQML messaging classes used by the agent are available at: https://github.com/bgyori/pykqml

INDRA is available at: https://github.com/sorgerlab/indra

Installing the TFTA
========================
Note that currently the TFTA has limited usage on its own. It's
meant to be launched in the context of a communication system. 

`pip install git+https://github.com/xzhang2016/tfagent.git`

Note
=======================
The master branch is using EKB XML as input format, while the dev branch is using JSON as input and output format. Enrichment is only at dev branch.
