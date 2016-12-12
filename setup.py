from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup

def main():
    setup(name='tfta',
          version='0.0.1',
          description='TF Agent',
          long_description='TF Agent',
          author='Xue Zhang',
          author_email='xue.zhang@tufts.edu',
          url='http://github.com/sorgerlab/bioagents',
          packages=['tfta'],
          install_requires=['pysb', 'indra', 'pykqml', 'objectpath', 'rdflib',
                            'functools32', 'requests', 'lxml',
                            'pandas', 'suds'],
          include_package_data=True,
          keywords=['systems', 'biology', 'model', 'pathway', 'assembler',
                    'nlp', 'mechanism', 'biochemistry'],
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          )
if __name__ == '__main__':
    main()
