#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Some python tools for the InSight mars mission.

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2019
'''
from os.path import join as pjoin

from setuptools import setup

setup(
    name='detect_events',
    version='0.1',
    packages=['detect_events'],
    package_data={'detect_events': [pjoin('data', '*')]},
    url='https://github.com/sstaehler/detect_events',
    license='GPLv3',
    author='Simon Staehler',
    author_email='mail@simonstaehler.com',
    description='Detect synthetic Marsquakes in real seismic noise',
    install_requires=['numpy', 'scipy', 'matplotlib', 'obspy', 'instaseis'],
    scripts=['bin/detect_marsquakes'],
    )
