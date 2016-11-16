#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'numpy',
    'openpyxl',
    'pyprind',
    'setuptools',
    'scipy',
    'xlrd',
    'xlwt',
]

test_requirements = [
    'coveralls',
    'flake8',
    'pytest',
    'pytest-cov',
    'pytest-flake8',
    'pytest-runner'
]

setup(
    name='pyRVT',
    version='0.5.6',
    description='Ground motion models implemented in Python.',
    long_description=readme + '\n\n' + history,
    author='Albert Kottke',
    author_email='albert.kottke@gmail.com',
    url='http://github.com/arkottke/pyrvt',
    packages=['pyrvt'],
    entry_points={
        'console_scripts': [
            'rvt_operator = pyrvt.runner:main',
            ],
        },
    package_dir={'pyrvt': 'pyrvt'},
    include_package_data=True,
    install_requires=requirements,
    license='MIT',
    zip_safe=False,
    keywords='pyrvt',
    package_data={'pyrvt': ['data/*'],},
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
