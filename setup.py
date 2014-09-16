#!/usr/bin/env python3
# encoding: utf-8

from setuptools import setup

config = {
    'name': 'pyRVT',
    'version': '0.1',
    'description': 'Seismologic Random Vibration Theory',
    'author': 'Albert Kottke',
    'author_email': 'albert.kottke@gmail.com',
    'url': 'http://github.com/arkottke/pyrvt',
    'entry_points': {
        'console_scripts': [
            'rvt_operator = pyrvt.runner:main',
            ],
        },
    'packages': ['pyrvt'],
    'package_data': {
        'pyrvt': ['data/*']
    },
    'requires': [
        'nose',
        'numpy',
        'scipy',
    ],
    'test_suite': 'nose.collector',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ],
}

setup(**config)
