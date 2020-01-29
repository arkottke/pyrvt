#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('README.rst') as fp:
    readme = fp.read()

with open('HISTORY.rst') as fp:
    history = fp.read()

setup(
    name='pyRVT',
    version='0.7.1',
    packages=find_packages(exclude=['.*tests.*']),
    entry_points={
        'console_scripts': [
            'pyrvt = pyrvt.runner:main',
        ],
    },
    install_requires=[
        'numpy',
        'numba',
        'setuptools',
        'scipy',
    ],
    extras_require={
        'xls': ['xlwt', 'xlrd'],
        'xlsx': ['openpyxl'],
        'all': ['openpyxl', 'xlwt', 'xlrd'],
    },
    package_data={
        'pyrvt': ['data/*'],
    },
    zip_safe=False,
    test_suite='tests',
    author='Albert Kottke',
    author_email='albert.kottke@gmail.com',
    description='Random vibration theory for earthquake ground motions.',
    license='MIT',
    long_description=readme + '\n\n' + history,
    url='http://github.com/arkottke/pyrvt',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research',
    ], )
