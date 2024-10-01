"""
Title: extract_lines.py
Author: Heather Sinclair-Wentworth, Quin Aicken Davies
Date: 30/09/24

Description: This script reads in the linelist from VisieR
and extracts the lines that are in the wavelength range
grouped by element.
"""
# %%
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.patches as pt
import pandas as pd
from astroquery.simbad import Simbad

#Taking GES linelist and splitting into a list per element

#%%

hdul = fits.open("/home/users/qai11/Documents/quin-masters-code/asu.fit")
data = hdul[1].data
data = data.byteswap().newbyteorder('=')
data = data.byteswap().newbyteorder('=')
data = data.byteswap().newbyteorder('=')
hdr1 = hdul[1].header
hdr0 = hdul[0].header


lines = pd.DataFrame(data)
lines = lines[(lines["gfflag"] != "-") & (lines["gfflag"] != "N") & (lines["synflag"] != "N") & (lines["synflag"] != "-")]
# lines = lines[(lines["lambda"].between(4800,6800))]
# lines = lines[~lines["lambda"].between(5776, 5836)]

lines.to_csv("/home/users/qai11/Documents/quin-masters-code/Linelists/ges_lines.tsv", index = False, sep = "\t")

#%%
# elements = ["Mg"]
elements = lines["Element"].unique()

all_list = pd.DataFrame()

for e in elements:
   el = lines[(lines["Element"] == e)].reset_index(drop=True)
   el_list = el[["lambda", "Element", "Ion"]]
   el_list["wave_peak"] = el_list["lambda"].round(4) / 10
   el_list["wave_base"] = (el_list["lambda"].round(4) / 10 - 0.05).round(4)
   el_list["wave_top"] = (el_list["lambda"].round(4) / 10 + 0.05).round(4)
   el_list["Ion"] = el_list["Ion"].astype(str)
   el_list["note"] = el_list["Element"] + " " + el_list["Ion"]
   el_list = el_list.drop(columns = ["Element", "Ion", "lambda"])
   
   el_list.to_csv(f"/home/users/qai11/Documents/quin-masters-code/Linelists/{e}_lines.csv", index = False, sep = "\t")
   all_list = pd.concat([all_list, el_list])
   all_list.to_csv(f"/home/users/qai11/Documents/quin-masters-code/Linelists/all_lines.csv", index = False, sep = "\t")
   

#%%
mg = lines[lines["Element"] == "Mg"].reset_index()
