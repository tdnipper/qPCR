#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:07:23 2023

@author: tnipper
"""
import pandas as pd

# Import files for each replicate from .xlsx and aggregate for statistics
rep1_input = input('Input filepath for replicate 1 input:')
rep1_input = pd.read_excel(rep1_input)
rep1_enrich = input('Rep1 enrich:')
rep1_enrich = pd.read_excel(rep1_enrich)
