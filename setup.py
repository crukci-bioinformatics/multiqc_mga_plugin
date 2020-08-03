#!/usr/bin/env python
"""
Multi Genome Alignment plugin for MultiQC.

For more information about MultiQC, see http://multiqc.info
"""

from setuptools import setup, find_packages

version = '1.6'

setup(
    name = 'multiqc_mga_plugin',
    version = version,
    author = 'Richard Bowers',
    author_email = 'richard.bowers@cruk.cam.ac.uk',
    description = "Multi Genome Alignment MultiQC plugin",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/crukci-bioinformatics/multiqc_mga_plugin',
    download_url = 'https://github.com/crukci-bioinformatics/multiqc_mga_plugin/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires = [
        'multiqc', 'lxml', 'natsort'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'mga = multiqc_mga.modules.mga:MultiqcModule',
        ],
        'multiqc.hooks.v1': [
            'execution_start = multiqc_mga.custom_code:mga_plugin_execution_start'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)
