"""
Title: write_ispec_to_moog.py
Author: Quin Aicken Davies
Date: 01/10/2024

Description: Uses the iSpec to convert the atomic linelist to a MOOG readable format.
"""

#%%
import numpy as np
import glob
import os
from astropy.io import fits
import pandas as pd
import time
import scipy as sp
import sys

#--- iSpec directory -------------------------------------------------------------
#ispec_dir = os.path.dirname(os.path.realpath(__file__)) + "/"
ispec_dir = '/home/users/qai11/iSpec_v20201001/'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

from ispec.lines import __moog_write_atomic_linelist
from ispec.lines import write_atomic_linelist
from ispec.lines import read_atomic_linelist

#%%

#--- Read atomic linelist -------------------------------------------------------

linelist_save_loc = '/home/users/qai11/Documents/Fixed_fits_files/quinlinelist.in'
linelist_filename = '/home/users/qai11/iSpec_v20201001/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv'
linelist = read_atomic_linelist(linelist_filename, wave_base=None, wave_top=None)

line_list = write_atomic_linelist(linelist, linelist_filename=linelist_save_loc, code='moog', tmp_dir=None)

# __moog_write_atomic_linelist(linelist, linelist_filename=None, tmp_dir=None)
# %%
