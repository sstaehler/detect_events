#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Some python tools for the InSight mars mission.

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2019
'''
from setuptools import setup, find_packages

setup(
    name='detect_events',
    version='0.1',
    packages=find_packages(),
    url='https://github.com/sstaehler/detect_events',
    license='GPLv3',
    author='Simon Staehler',
    author_email='mail@simonstaehler.com',
    description='Detect synthetic Marsquakes in real seismic noise'
    )
