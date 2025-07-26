"""
Title: make_iSpec_linelist.py
Author: Quin Aicken Davies
Date: 13/09/24

Description: This script reads in the linelist made by linemake.f
and a iSpec linelist, Reformat the linemake linelist to be added
to the iSpec linelist. This is to be used in the iSpec analysis and
gets converted to the moog list in moog.py saved in documents.
"""
#%%
import pandas as pd
import numpy as np
import os
import sys
from astropy.io import fits
import matplotlib.pyplot as plt
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
import matplotlib.patches as pt
from astroquery.simbad import Simbad

#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

from ispec import read_chemical_elements

from ispec import read_molecular_symbols


#%%
molecules = read_molecular_symbols(ispec_dir + '/input/abundances/molecular_symbols.dat')
elements = read_chemical_elements(ispec_dir + '/input/abundances/chemical_elements_symbols.dat')
elements = elements.filled()
elements = np.array(elements, dtype=[('id', 'i4'), ('symbol', 'U2'), ('name', 'U20'), ('group', 'i4'), ('period', 'i4'), ('block', 'i4'), ('atomic_num', 'i4'), ('ion', 'i4')])

# Load both linelists
linelist_iSpec = pd.read_csv(ispec_dir + '/input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv', delimiter='\t', engine='python')
try:
    linelist_quin = pd.read_csv('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/quinlist.MgH', delimiter=r'\s+', engine='python')
    linelist_quin2 = pd.read_csv('/home/users/qai11/Documents/quin-masters-code/quinlist5500.MgH', delimiter=r'\s+', engine='python')
except:
    linelist_quin = pd.read_csv('/Users/quin/quin-masters-code/quinlist.MgH', delimiter=r'\s+', engine='python')

# Convert the molecule array to a dictionary for easy mapping
molecule_dict = {item['atomic_num']: item['symbol'] for item in molecules}

#Rename the columns for concatenating the 5100-5200 linelists
linelist_quin.columns = ['wave_A', 'element', 'lower_state_eV', 'loggf', 'waals', 'd0', 'equivalent_width', 'comment']

#Rename the columns for concatenating the 5500-5200 linelists
linelist_quin2.columns = ['wave_A', 'element', 'lower_state_eV', 'loggf', 'waals', 'd0', 'equivalent_width', 'comment']

# Add a new column to track the source of each DataFrame
linelist_quin['source'] = 'quin'
linelist_quin2['source'] = 'quin2'
linelist_iSpec['source'] = 'iSpec'

# Concatenate the DataFrames
merged_df = pd.concat([linelist_iSpec, linelist_quin, linelist_quin2], axis=0, ignore_index=True)

# Fill in missing values with 0 for the non-overlapping columns
merged_df = merged_df.fillna(0)

# Set 'theoretical_depth' to 0.01 for rows added from 'linelist_quin' and 'linelist_quin2'
merged_df.loc[merged_df['source'].isin(['quin', 'quin2']), 'theoretical_depth'] = 0.01

# Drop the 'source' column if it's no longer needed
merged_df = merged_df.drop(columns='source')

# Update the wave_nm column for the quin linelist section
merged_df.loc[merged_df['wave_nm'] == 0, 'wave_nm'] = merged_df['wave_A'] / 10

# Convert non-numeric values in 'waals_single_gamma_format' to NaN
merged_df['waals'] = pd.to_numeric(merged_df['waals'], errors='coerce')
# Replace NaN values with 0
merged_df['waals'].fillna(0, inplace=True)


# Update values in 'waals_single_gamma_format' with values from 'waals' where 'waals_single_gamma_format' is 0
merged_df.loc[merged_df['waals_single_gamma_format'] == 0, 'waals_single_gamma_format'] = merged_df['waals']

# Update spectrum_fudge_factor values of 0 to 1
merged_df.loc[merged_df['spectrum_fudge_factor'] == 0, 'spectrum_fudge_factor'] = 1
# Update lower_state_cm1 values of 0 to the calculated value for eV to cm^-1
merged_df.loc[merged_df['lower_state_cm1'] == 0, 'lower_state_cm1'] = merged_df['lower_state_eV'] * 8065.54429

# Function to check if the first three digits are in the molecule_dict
def check_molecule(element_value):
    try:
        # Try to extract the first three digits (convert to int)
        first_three_digits = int(str(element_value).split('.')[0])
        # Check if the first three digits are in the dictionary
        if first_three_digits in molecule_dict:
            return 'T'
        else:
            return 'F'
    except ValueError:
        # If the value cannot be converted to an integer (e.g., 'La 1'), return 'F'
        return 'F'

# Apply the function to create the 'molecule' column
merged_df['molecule'] = merged_df['element'].apply(check_molecule)

# Update the 'spectrum_moog_species' with the values in 'element' with moog format
merged_df.loc[merged_df['spectrum_moog_species'] == 0, 'spectrum_moog_species'] = merged_df['element']

# Function to check if the first three digits are in the molecule_dict
def check_ion(element_value):
    try:
        # Try to extract the first three digits (convert to int)
        first_digit_after_decimal = int(str(element_value).split('.')[1][0])
        # If the first digit after the decimal is 0-6, return the corresponding ion number
        if first_digit_after_decimal == 0:
            return 1
        elif first_digit_after_decimal == 1:
            return 2
        elif first_digit_after_decimal == 2:
            return 3
        elif first_digit_after_decimal == 3:
            return 4
        elif first_digit_after_decimal == 4:
            return 5
        elif first_digit_after_decimal == 5:
            return 6
    except ValueError:
        # If the value cannot be converted to an integer (e.g., 'La 1'), return 0
        return 0
# Apply the function to create the 'ion' column
merged_df['ion'] = merged_df['spectrum_moog_species'].apply(check_ion)

# Change the 'moog_support' column to 'T' where it is 0 from the moog linelist
merged_df.loc[merged_df['moog_support'] == 0, 'moog_support'] = "T"

# Function to map numbers to molecule names
def map_to_molecule(value):
    if isinstance(value, str):
        # Skip values that are already strings (like "La 1")
        return value
    else:
        # Extract the integer part before the decimal
        key = int(str(value).split('.')[0])
        # Return the corresponding molecule name from the dictionary or the original value if not found
        return molecule_dict.get(key, value)

# Apply the function to the 'elements' column
merged_df['element'] = merged_df['element'].apply(map_to_molecule)

# Create a dictionary from the ndarray for quick lookup
element_dict = {row['id']: (row['symbol'], row['ion']) for row in elements}

# Function to map numbers to element names with ion numbers
def map_to_element(value, ion):
    if isinstance(value, str):
        # Skip values that are already strings (like "La 1")
        return value
    else:
        # Extract the integer part and the first three digits before the decimal
        key = int(str(value).split('.')[0])
        # Get the element symbol and ion number
        symbol, ion_number = element_dict.get(key, (value, ion))
        # Return the formatted string
        return f"{symbol} {ion}".strip()

# Apply the function to the 'elements' column, passing the 'ion' column for context
merged_df['element'] = merged_df.apply(lambda row: map_to_element(row['element'], row['ion']), axis=1)

#Remove the last 3 columns
merged_df = merged_df.iloc[:, :-3]



# %%
"""ADD in the EuII isotopes list to the linelist from the GES linelist"""
# Load the EuII isotopes list
try:
    hdul = fits.open("/home/users/qai11/Documents/quin-masters-code/asu.fit")
except:
    hdul = fits.open("/Users/quin/quin-masters-code/asu.fit")
    
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

# lines.to_csv("/home/users/qai11/Documents/quin-masters-code/Linelists/ges_lines.tsv", index = False, sep = "\t")

# elements = ["Mg"]
elements = lines["Element"].unique()

#%%
# Now I want to reformat the lines to match the iSpec format

# Create a new DataFrame to hold the reformatted lines
el_list_Eu = pd.DataFrame()

'''Pull out the Eu lines I want to reformat them to match iSpec'''
el_Eu = lines[(lines["Element"] == "Eu")].reset_index(drop=True)

# Add the correct values to the new DataFrame for Europium
el_list_Eu['element'] = el_Eu['Element']+ " " + el_Eu['Ion'].astype(str)
el_list_Eu['wave_A'] = el_Eu['lambda']
el_list_Eu['wave_nm'] = el_Eu['lambda'] / 10
el_list_Eu['loggf'] = el_Eu['loggf']
el_list_Eu['lower_state_eV'] = el_Eu['Elow']
el_list_Eu['lower_state_cm1'] = el_Eu['Elow'] * 8065.54429
el_list_Eu['lower_j'] = el_Eu['Jlow']
el_list_Eu['upper_state_eV'] = el_Eu['Eup']
el_list_Eu['upper_state_cm1'] = el_Eu['Eup'] * 8065.54429
el_list_Eu['upper_j'] = el_Eu['Jup']
el_list_Eu['ion'] = el_Eu['Ion'].astype(int)
el_list_Eu['moog_support'] = 'T'
el_list_Eu['waals'] = el_Eu['Vdw-damp']
el_list_Eu['stark'] = -7.0
el_list_Eu['rad'] = 8
el_list_Eu['spectrum_moog_species'] = round(63.1 + (el_Eu['Isotope']  * 0.0001), 4).astype(float)
el_list_Eu['theoretical_depth'] = 0.01

# then when we are done we want to set any other value to be zero because they arent needed

#do the same formatting for Barium
'''Pull out the Ba lines I want to reformat them to match iSpec'''
el_Ba = lines[(lines["Element"] == "Ba")].reset_index(drop=True)

# Add the correct values to the new DataFrame for Barium
el_list_Ba = pd.DataFrame()

el_list_Ba['element'] = el_Ba['Element']+ " " + el_Ba['Ion'].astype(str)
el_list_Ba['wave_A'] = el_Ba['lambda']
el_list_Ba['wave_nm'] = el_Ba['lambda'] / 10
el_list_Ba['loggf'] = el_Ba['loggf']
el_list_Ba['lower_state_eV'] = el_Ba['Elow']
el_list_Ba['lower_state_cm1'] = el_Ba['Elow'] * 8065.54429
el_list_Ba['lower_j'] = el_Ba['Jlow']
el_list_Ba['upper_state_eV'] = el_Ba['Eup']
el_list_Ba['upper_state_cm1'] = el_Ba['Eup'] * 8065.54429
el_list_Ba['upper_j'] = el_Ba['Jup']
el_list_Ba['ion'] = el_Ba['Ion'].astype(int)
el_list_Ba['moog_support'] = 'T'
el_list_Ba['waals'] = el_Ba['Vdw-damp']
el_list_Ba['stark'] = -7.0
el_list_Ba['rad'] = 8
el_list_Ba['spectrum_moog_species'] = round(56.1 + (el_Ba['Isotope'] * 0.0001), 4).astype(float)
el_list_Ba['theoretical_depth'] = 0.01


#%%
# Concatenate the two DataFrames, keeping the original data intact and adding zeros for missing values
merged_df_new = pd.concat([merged_df, el_list_Eu, el_list_Ba], axis=0, ignore_index=True)

merged_df_new = merged_df_new.fillna(0)

#Sort the list by wavelength
merged_df_new = merged_df_new.sort_values(by='wave_A')

#Create folder for the new linelist in iSpec directories
if not os.path.exists(ispec_dir + '/input/linelists/transitions/Quin_GES_LIST.420_920nm/'):
    os.makedirs(ispec_dir +'/input/linelists/transitions/Quin_GES_LIST.420_920nm/')
    
#Save the new linelist as a tsv file called atomic_lines.tsv
merged_df_new.to_csv(ispec_dir + '/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv', sep='\t', index=False)

# %%
