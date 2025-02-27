"""
Python installation file for welly project.

:copyright: 2021 Agile Scientific
:license: Apache 2.0
"""
from setuptools import setup
import re

verstr = 'unknown'
VERSIONFILE = "welly/_version.py"
with open(VERSIONFILE, "r") as f:
    verstrline = f.read().strip()
    pattern = re.compile(r"__version__ = ['\"](.*)['\"]")
    mo = pattern.search(verstrline)
if mo:
    verstr = mo.group(1)
    print("Version "+verstr)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

REQUIREMENTS = ['numpy',
                'scipy',
                'matplotlib',
                'lasio',
                'striplog',
                'tqdm',
                ]

TEST_REQUIREMENTS = ['pytest',
                     'coveralls',
                     'pytest-cov',
                     'pytest-mpl',
                     ]

# Test command is:
# python run_tests.py

CLASSIFIERS = ['Development Status :: 4 - Beta',
               'Intended Audience :: Science/Research',
               'Natural Language :: English',
               'License :: OSI Approved :: Apache Software License',
               'Operating System :: OS Independent',
               'Programming Language :: Python',
               'Programming Language :: Python :: 3.6',
               'Programming Language :: Python :: 3.7',
               'Programming Language :: Python :: 3.8',
               'Programming Language :: Python :: 3.9',
               ]

setup(name='welly',
      version=verstr,
      description='Tools for making and managing well data.',
      url='http://github.com/agile-geoscience/welly',
      author='Agile Scientific',
      author_email='hello@agilescientific.com',
      license='Apache 2',
      packages=['welly'],
      tests_require=TEST_REQUIREMENTS,
      test_suite='run_tests',
      install_requires=REQUIREMENTS,
      classifiers=CLASSIFIERS,
      zip_safe=False,
      )
