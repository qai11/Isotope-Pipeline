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
molecules = read_molecular_symbols('/home/users/qai11/iSpec_v20201001/input/abundances/molecular_symbols.dat')
elements = read_chemical_elements('/home/users/qai11/iSpec_v20201001/input/abundances/chemical_elements_symbols.dat')
elements = elements.filled()
elements = np.array(elements, dtype=[('id', 'i4'), ('symbol', 'U2'), ('name', 'U20'), ('group', 'i4'), ('period', 'i4'), ('block', 'i4'), ('atomic_num', 'i4'), ('ion', 'i4')])

# Load both linelists
linelist_iSpec = pd.read_csv(ispec_dir + '/input/linelists/transitions/GESv6_atom_hfs_iso.420_920nm/atomic_lines.tsv', delimiter='\t', engine='python')
linelist_quin = pd.read_csv('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/quinlist.MgH', delimiter=r'\s+', engine='python')

# Convert the molecule array to a dictionary for easy mapping
molecule_dict = {item['atomic_num']: item['symbol'] for item in molecules}

linelist_quin.columns = ['wave_A', 'element', 'lower_state_eV', 'loggf', 'waals', 'd0', 'equivalent_width', 'comment']

# Concatenate the two DataFrames, keeping the original data intact
merged_df = pd.concat([linelist_iSpec, linelist_quin], axis=0, ignore_index=True)

# Fill in missing values with 0 for the non-overlapping columns
merged_df = merged_df.fillna(0)

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

#Sort the list by wavelength
merged_df = merged_df.sort_values(by='wave_A')

#Create folder for the new linelist in iSpec directories
if not os.path.exists(ispec_dir + '/input/linelists/transitions/Quin_GES_LIST.420_920nm/'):
    os.makedirs(ispec_dir +'/input/linelists/transitions/Quin_GES_LIST.420_920nm/')
    
#Save the new linelist as a tsv file called atomic_lines.tsv
merged_df.to_csv(ispec_dir + '/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv', sep='\t', index=False)

# %%
 