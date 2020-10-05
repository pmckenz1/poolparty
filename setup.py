#!/usr/bin/env python

import os
import re
from setuptools import setup

# parse version from init.py
with open("poolparty/__init__.py") as init:
    CUR_VERSION = re.search(
        r"^__version__ = ['\"]([^'\"]*)['\"]",
        init.read(),
        re.M,
    ).group(1)


# setup installation
setup(
    name="poolparty",
    packages=["poolparty"],
    version=CUR_VERSION,
    author="Patrick McKenzie",
    author_email="p.mckenzie@columbia.edu",
    install_requires=["numpy"],
    license='MIT',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',                
    ],
)
