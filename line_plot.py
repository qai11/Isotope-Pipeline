
from matplotlib import pyplot as plt
import pandas as pd
import os
import sys

#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

star_name = 'hd_45588'

synth = pd.read_csv(f"/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/synth_regions_Ba/synth_Ba_585.36663.txt",delimiter=' ')

star_spectrum = ispec.read_spectrum(f"/home/users/qai11/Documents/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")

plt.figure()
plt.plot(synth['waveobs'],synth['flux'])
plt.plot(star_spectrum['waveobs'],star_spectrum['flux'])
plt.show()