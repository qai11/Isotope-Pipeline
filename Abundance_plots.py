""" Abundance_plots.py"""

# %%  for X/Fe

import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XFe(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/Fe]', 'e[Eu/Fe]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/Fe]', 'e[Ba/Fe]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/Fe]', 'e[Mg/Fe]'])
    
    #Open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    
    for star_name in star_list:
        #Extract the FEH value from the masters stars csv
        feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        #Take away feh from each X/H value
        eu_row['[X/H]'] = eu_row['[X/H]'] - feh
        ba_row['[X/H]'] = ba_row['[X/H]'] - feh
        mg_row['[X/H]'] = mg_row['[X/H]'] - feh

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/Fe]', 'e[Eu/Fe]']
        ba_row.columns = ['star_name', '[Ba/Fe]', 'e[Ba/Fe]']
        mg_row.columns = ['star_name', '[Mg/Fe]', 'e[Mg/Fe]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        
    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the 'mg_fe', 'd_mg_fe', 'mg24_fe', 'd_mg24_fe', 'mg25_fe', 'd_mg25_fe', 'mg26_fe', 'd_mg26_fe columns
    iso_mg = isotope[['mg_fe', 'd_mg_fe', 'mg_fe24', 'd_mg_fe24', 'mg_fe25', 'd_mg_fe25', 'mg_fe26', 'd_mg_fe26']]
    # print(iso_mg['mg_fe']-isotope['mg'])
    # plot with X/Fe and IS Mg
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,15))
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'])}")
    ax[0,0].errorbar(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[0,0].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[0,0].set_ylim(-0.2,1.1)
    #Plot Ba vs Mg
    print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'])}")
    ax[0,1].errorbar(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[0,1].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Mg vs iso_mg
    print(f"Mg vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'])}")
    ax[1,2].errorbar(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o')
    ax[1,2].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[1,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Eu vs iso_mg
    print(f"Eu vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Eu_values['[Eu/Fe]'])}")
    ax[1,0].errorbar(iso_mg['mg_fe'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o')
    ax[1,0].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[1,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[1,0].set_ylim(-0.2,1.1)
    #Plot Ba vs iso_mg
    print(f"Ba vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Ba_values['[Ba/Fe]'])}")
    ax[1,1].errorbar(iso_mg['mg_fe'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o')
    ax[1,1].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[1,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Eu vs mg24
    print(f"Eu vs mg24: {scipy.stats.pearsonr(iso_mg['mg_fe24'], Eu_values['[Eu/Fe]'])}")
    ax[2,0].errorbar(iso_mg['mg_fe24'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_mg_fe24'], fmt='o')
    ax[2,0].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[2,0].set_ylim(-0.2,1.1)
    #Plot Ba vs mg24
    print(f"Ba vs mg24: {scipy.stats.pearsonr(iso_mg['mg_fe24'], Ba_values['[Ba/Fe]'])}")
    ax[2,1].errorbar(iso_mg['mg_fe24'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_mg_fe24'], fmt='o')
    ax[2,1].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Eu vs mg25
    print(f"Eu vs mg25: {scipy.stats.pearsonr(iso_mg['mg_fe25'], Eu_values['[Eu/Fe]'])}")
    ax[3,0].errorbar(iso_mg['mg_fe25'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_mg_fe25'], fmt='o')
    ax[3,0].set_xlabel('IS[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[3,0].set_ylim(-0.2,1.1)
    #Plot Ba vs mg25
    print(f"Ba vs mg25 {scipy.stats.pearsonr(iso_mg['mg_fe25'], Ba_values['[Ba/Fe]'])}")
    ax[3,1].errorbar(iso_mg['mg_fe25'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_mg_fe25'], fmt='o')
    ax[3,1].set_xlabel('IS[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Eu vs mg26
    print(f"Eu vs mg26: {scipy.stats.pearsonr(iso_mg['mg_fe26'], Eu_values['[Eu/Fe]'])}")
    ax[4,0].errorbar(iso_mg['mg_fe26'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_mg_fe26'], fmt='o')
    ax[4,0].set_xlabel('IS[$^{26}$Mg/Fe]',fontsize=12)
    ax[4,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[4,0].set_ylim(-0.2,1.1)
    #Plot Ba vs mg26
    print(f"Ba vs mg26 {scipy.stats.pearsonr(iso_mg['mg_fe26'], Ba_values['[Ba/Fe]'])}")
    ax[4,1].errorbar(iso_mg['mg_fe26'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_mg_fe26'], fmt='o')
    ax[4,1].set_xlabel('IS[$^{26}$Mg/Fe]',fontsize=12)
    ax[4,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg_fe24'], Mg_values['[Mg/Fe]'])}")
    ax[2,2].errorbar(iso_mg['mg_fe24'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe24'], fmt='o')
    ax[2,2].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Mg vs mg25
    print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['mg_fe25'], Mg_values['[Mg/Fe]'])}")
    ax[3,2].errorbar(iso_mg['mg_fe25'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe25'], fmt='o')
    ax[3,2].set_xlabel('IS[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Mg vs mg26
    print(f"Mg vs mg26: {scipy.stats.pearsonr(iso_mg['mg_fe26'], Mg_values['[Mg/Fe]'])}")
    ax[4,2].errorbar(iso_mg['mg_fe26'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe26'], fmt='o')
    ax[4,2].set_xlabel('IS[$^{26}$Mg/Fe]',fontsize=12)
    ax[4,2].set_ylabel('[Mg/Fe]',fontsize=12)

    #plot difference plot for Mg vs iso mg
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'] - iso_mg['mg_fe'])}")
    ax[0,2].errorbar(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'] - iso_mg['mg_fe'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o')
    ax[0,2].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/Fe] - IS[Mg/Fe]',fontsize=12)
    
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_Fe.png', dpi=300, bbox_inches='tight')

element_plots_XFe(star_list)

#%%  plots for X/H
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XH(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,15))
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'])}")
    ax[0,0].errorbar(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o')
    ax[0,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,0].set_ylim(-0.7,0.7)
    #Plot Ba vs Mg
    print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'])}")
    ax[0,1].errorbar(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o')
    ax[0,1].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Mg vs iso_mg
    print(f"Mg vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg'], Mg_values['[Mg/H]'])}")
    ax[1,2].errorbar(iso_mg['mg'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg'], fmt='o')
    ax[1,2].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[1,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Eu vs iso_mg
    print(f"Eu vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg'], Eu_values['[Eu/H]'])}")
    ax[1,0].errorbar(iso_mg['mg'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_mg'], fmt='o')
    ax[1,0].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[1,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[1,0].set_ylim(-0.7,0.7)
    #Plot Ba vs iso_mg
    print(f"Ba vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg'], Ba_values['[Ba/H]'])}")
    ax[1,1].errorbar(iso_mg['mg'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_mg'], fmt='o')
    ax[1,1].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[1,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg24
    print(f"Eu vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Eu_values['[Eu/H]'])}")
    ax[2,0].errorbar(iso_mg['mg24'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_mg24'], fmt='o')
    ax[2,0].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[2,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[2,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg24
    print(f"Ba vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Ba_values['[Ba/H]'])}")
    ax[2,1].errorbar(iso_mg['mg24'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_mg24'], fmt='o')
    ax[2,1].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[2,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg25
    print(f"Eu vs mg25: {scipy.stats.pearsonr(iso_mg['mg25'], Eu_values['[Eu/H]'])}")
    ax[3,0].errorbar(iso_mg['mg25'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_mg25'], fmt='o')
    ax[3,0].set_xlabel('IS[$^{25}$Mg/H]',fontsize=12)
    ax[3,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[3,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg25
    print(f"Ba vs mg25: {scipy.stats.pearsonr(iso_mg['mg25'], Ba_values['[Ba/H]'])}")
    ax[3,1].errorbar(iso_mg['mg25'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_mg25'], fmt='o')
    ax[3,1].set_xlabel('IS[$^{25}$Mg/H]',fontsize=12)
    ax[3,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg26
    print(f"Eu vs mg26: {scipy.stats.pearsonr(iso_mg['mg26'], Eu_values['[Eu/H]'])}")
    ax[4,0].errorbar(iso_mg['mg26'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_mg26'], fmt='o')
    ax[4,0].set_xlabel('IS[$^{26}$Mg/H]',fontsize=12)
    ax[4,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[4,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg26
    print(f"Ba vs mg26: {scipy.stats.pearsonr(iso_mg['mg26'], Ba_values['[Ba/H]'])}")
    ax[4,1].errorbar(iso_mg['mg26'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_mg26'], fmt='o')
    ax[4,1].set_xlabel('IS[$^{26}$Mg/H]',fontsize=12)
    ax[4,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Mg_values['[Mg/H]'])}")
    ax[2,2].errorbar(iso_mg['mg24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg24'], fmt='o')
    ax[2,2].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg vs mg25
    print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['mg25'], Mg_values['[Mg/H]'])}")
    ax[3,2].errorbar(iso_mg['mg25'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg25'], fmt='o')
    ax[3,2].set_xlabel('IS[$^{25}$Mg/H]',fontsize=12)
    ax[3,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg vs mg26
    print(f"Mg vs mg26: {scipy.stats.pearsonr(iso_mg['mg26'], Mg_values['[Mg/H]'])}")
    ax[4,2].errorbar(iso_mg['mg26'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg26'], fmt='o')
    ax[4,2].set_xlabel('IS[$^{26}$Mg/H]',fontsize=12)
    ax[4,2].set_ylabel('[Mg/H]',fontsize=12)

    #plot difference plot for Mg vs iso mg
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(iso_mg['mg'], Mg_values['[Mg/H]'] - iso_mg['mg'])}")
    ax[0,2].errorbar(iso_mg['mg'], Mg_values['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]']-iso_mg['d_mg'], xerr=iso_mg['d_mg'], fmt='o')
    ax[0,2].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/H] - IS[Mg/H]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H.png', dpi=300, bbox_inches='tight')



element_plots_XH(star_list)

#%%  X/Fe with Teff logg and Bavs Eu plots
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def single_el_vs_params(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    Eu_values_H = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values_H = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values_H = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    
    #Open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #remove the 10th row
    star_info = star_info.drop(10)
    #reset index
    star_info = star_info.reset_index(drop=True)
    
    for star_name in star_list:
        #Extract the FEH value from the masters stars csv
        feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        #EUH
        eu_row_H = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row_H = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row_H = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        
        #Take away feh from each X/H value
        eu_row['[X/H]'] = eu_row['[X/H]'] - feh
        ba_row['[X/H]'] = ba_row['[X/H]'] - feh
        mg_row['[X/H]'] = mg_row['[X/H]'] - feh
        
        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)
        eu_row_H.insert(0, 'star_name', star_name)
        ba_row_H.insert(0, 'star_name', star_name)
        mg_row_H.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/Fe]', 'e[Eu/Fe]']
        ba_row.columns = ['star_name', '[Ba/Fe]', 'e[Ba/Fe]']
        mg_row.columns = ['star_name', '[Mg/Fe]', 'e[Mg/Fe]']
        eu_row_H.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row_H.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row_H.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        Eu_values_H = pd.concat([Eu_values_H, eu_row_H], ignore_index=True)
        Ba_values_H = pd.concat([Ba_values_H, ba_row_H], ignore_index=True)
        Mg_values_H = pd.concat([Mg_values_H, mg_row_H], ignore_index=True)
    
    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg_fe', 'd_mg_fe','mg','d_mg']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[iso_mg['mg'] < 0.25]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values_H[Mg_values_H['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values_H[Eu_values_H['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values_H[Ba_values_H['star_name'].isin(small_mg['Unnamed: 0'])]
    small_star_teff = star_info[star_info['ID2'].isin(small_mg['Unnamed: 0'])]
    small_star_logg = star_info[star_info['ID2'].isin(small_mg['Unnamed: 0'])]


    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(3, 4, figsize=(16,12),constrained_layout=True)

    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    # print(scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]']))
    ax[0,0].errorbar(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o', elinewidth=0.5)
    ax[0,0].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[0,0].set_ylim(-0.4,0.9)
    #Plot Eu vs Ba
    print(f"Eu vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], Eu_values['[Eu/Fe]'])}")
    ax[0,1].errorbar(Ba_values['[Ba/Fe]'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[0,1].set_xlabel('[Ba/Fe]',fontsize=12)
    ax[0,1].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[0,1].set_ylim(-0.4,1)    
    #Plot Ba vs Mg
    ax[0,2].errorbar(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o', elinewidth=0.5)
    ax[0,2].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,2].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Eu vs Teff
    print(f"Eu vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Eu_values['[Eu/Fe]'])}")
    ax[1,0].errorbar(star_info['TEFF'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    ax[1,0].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[1,0].set_ylim(-0.4,0.9)
    #Plot Ba vs Teff
    print(f"Ba vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Ba_values['[Ba/Fe]'])}")
    ax[1,1].errorbar(star_info['TEFF'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[1,1].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Mg difference vs Teff
    print(f"Mg diff vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Mg_values['[Mg/Fe]'])}")
    ax[1,3].errorbar(star_info['TEFF'], Mg_values_H['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values_H['e[Mg/H]']-iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[1,3].errorbar(small_star_teff['TEFF'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,3].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,3].set_ylabel('[Mg/H] - IS[Mg/H]',fontsize=12)
    
    #Plot Eu vs logg
    print(f"Eu vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Eu_values['[Eu/Fe]'])}")
    ax[2,0].errorbar(star_info['LOGG'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    ax[2,0].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[2,0].set_ylim(-0.4,0.9)
    #Plot Ba vs logg
    print(f"Ba vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Ba_values['[Ba/Fe]'])}")
    ax[2,1].errorbar(star_info['LOGG'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[2,1].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,1].set_ylabel('[Ba/Fe]',fontsize=12)
    
    #Plot Mg differecne vs logg
    print(f"Mg diff vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Mg_values['[Mg/Fe]'] - iso_mg['mg_fe'])}")
    ax[2,3].errorbar(star_info['LOGG'], Mg_values_H['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values_H['e[Mg/H]']-iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[2,3].errorbar(small_star_logg['LOGG'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,3].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,3].set_ylabel('[Mg/H] - IS[Mg/H]',fontsize=12)
    
    # Plot Mg diffvs isoMg
    ax[0,3].errorbar(iso_mg['mg'],  Mg_values_H['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]']-iso_mg['d_mg'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[0,3].errorbar(small_mg['mg'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    ax[0,3].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[0,3].set_ylabel('[Mg/H]-IS[Mg/H]',fontsize=12)
    #Plot Mg vs Teff
    print(f"Mg vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Mg_values['[Mg/Fe]'])}")
    ax[1,2].errorbar(star_info['TEFF'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], fmt='o', elinewidth=0.5)
    ax[1,2].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Mg vs logg
    print(f"Mg vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Mg_values['[Mg/Fe]'])}")
    ax[2,2].errorbar(star_info['LOGG'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], fmt='o', elinewidth=0.5)
    ax[2,2].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_vFe_params.png', dpi=300, bbox_inches='tight')



single_el_vs_params(star_list) 
#%% X/H with ratios made from lbl Mg
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XFe_LMG(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26','MgH','MgH24','MgH25','MgH26',
                      'd_MgH','d_MgH24','d_MgH25','d_MgH26']]

    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,15))
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'])}")
    ax[0,0].errorbar(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o')
    ax[0,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,0].set_ylim(-0.7,0.7)
    #Plot Ba vs Mg
    print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'])}")
    ax[0,1].errorbar(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o')
    ax[0,1].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Mg vs iso_mg
    print(f"Mg vs iso_mg: {scipy.stats.pearsonr(iso_mg['MgH'], Mg_values['[Mg/H]'])}")
    ax[1,2].errorbar(iso_mg['MgH'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH'], fmt='o')
    ax[1,2].set_xlabel('[Mg/H]',fontsize=12)
    ax[1,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Eu vs iso_mg
    print(f"Eu vs iso_mg: {scipy.stats.pearsonr(iso_mg['MgH'], Eu_values['[Eu/H]'])}")
    ax[1,0].errorbar(iso_mg['MgH'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH'], fmt='o')
    ax[1,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[1,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[1,0].set_ylim(-0.7,0.7)
    #Plot Ba vs iso_mg
    print(f"Ba vs iso_mg: {scipy.stats.pearsonr(iso_mg['MgH'], Ba_values['[Ba/H]'])}")
    ax[1,1].errorbar(iso_mg['MgH'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_MgH'], fmt='o')
    ax[1,1].set_xlabel('[Mg/H]',fontsize=12)
    ax[1,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg24
    print(f"Eu vs mg24: {scipy.stats.pearsonr(iso_mg['MgH24'], Eu_values['[Eu/H]'])}")
    ax[2,0].errorbar(iso_mg['MgH24'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH24'], fmt='o')
    ax[2,0].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[2,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[2,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg24
    print(f"Ba vs mg24: {scipy.stats.pearsonr(iso_mg['MgH24'], Ba_values['[Ba/H]'])}")
    ax[2,1].errorbar(iso_mg['MgH24'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_MgH24'], fmt='o')
    ax[2,1].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[2,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg25
    print(f"Eu vs mg25: {scipy.stats.pearsonr(iso_mg['MgH25'], Eu_values['[Eu/H]'])}")
    ax[3,0].errorbar(iso_mg['MgH25'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH25'], fmt='o')
    ax[3,0].set_xlabel('[$^{25}$Mg/H]',fontsize=12)
    ax[3,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[3,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg25
    print(f"Ba vs mg25: {scipy.stats.pearsonr(iso_mg['MgH25'], Ba_values['[Ba/H]'])}")
    ax[3,1].errorbar(iso_mg['MgH25'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_MgH25'], fmt='o')
    ax[3,1].set_xlabel('[$^{25}$Mg/H]',fontsize=12)
    ax[3,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg26
    print(f"Eu vs mg26: {scipy.stats.pearsonr(iso_mg['MgH26'], Eu_values['[Eu/H]'])}")
    ax[4,0].errorbar(iso_mg['MgH26'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH26'], fmt='o')
    ax[4,0].set_xlabel('[$^{26}$Mg/H]',fontsize=12)
    ax[4,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[4,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg26
    print(f"Ba vs mg26: {scipy.stats.pearsonr(iso_mg['MgH26'], Ba_values['[Ba/H]'])}")
    ax[4,1].errorbar(iso_mg['MgH26'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_MgH26'], fmt='o')
    ax[4,1].set_xlabel('[$^{26}$Mg/H]',fontsize=12)
    ax[4,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['MgH24'], Mg_values['[Mg/H]'])}")
    ax[2,2].errorbar(iso_mg['MgH24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH24'], fmt='o')
    ax[2,2].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg vs mg25
    print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['MgH25'], Mg_values['[Mg/H]'])}")
    ax[3,2].errorbar(iso_mg['MgH25'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH25'], fmt='o')
    ax[3,2].set_xlabel('[$^{25}$Mg/H]',fontsize=12)
    ax[3,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg vs mg26
    print(f"Mg vs mg26: {scipy.stats.pearsonr(iso_mg['MgH26'], Mg_values['[Mg/H]'])}")
    ax[4,2].errorbar(iso_mg['MgH26'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH26'], fmt='o')
    ax[4,2].set_xlabel('[$^{26}$Mg/H]',fontsize=12)
    ax[4,2].set_ylabel('[Mg/H]',fontsize=12)

    #plot difference plot for Mg vs iso mg
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(iso_mg['MgH'], Mg_values['[Mg/H]'] - iso_mg['MgH'])}")
    ax[0,2].errorbar(iso_mg['MgH'], Mg_values['[Mg/H]'] - iso_mg['MgH'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH'], fmt='o')
    ax[0,2].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/H] - [Mg/H]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_LblMg.png', dpi=300, bbox_inches='tight')

element_plots_XFe_LMG(star_list)

#%% X/Fe with ratios made from lbl Mg
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XFe_LMg(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/Fe]', 'e[Eu/Fe]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/Fe]', 'e[Ba/Fe]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/Fe]', 'e[Mg/Fe]'])

    #Open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    
    for star_name in star_list:
        #Extract the FEH value from the masters stars csv
        feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        #Take away feh from each X/H value
        eu_row['[X/H]'] = eu_row['[X/H]'] - feh
        ba_row['[X/H]'] = ba_row['[X/H]'] - feh
        mg_row['[X/H]'] = mg_row['[X/H]'] - feh

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        
        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/Fe]', 'e[Eu/Fe]']
        ba_row.columns = ['star_name', '[Ba/Fe]', 'e[Ba/Fe]']
        mg_row.columns = ['star_name', '[Mg/Fe]', 'e[Mg/Fe]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
    # print(Mg_values)
    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26','MgH','MgH24','MgH25','MgH26',
                      'd_MgH','d_MgH24','d_MgH25','d_MgH26','MgFe','MgFe24','MgFe25','MgFe26',
                      'd_MgFe','d_MgFe24','d_MgFe25','d_MgFe26']]
    
    # print(Mg_values['[Mg/Fe]']*(isotope['R_26']/100))
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,15))
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'])}")
    ax[0,0].errorbar(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[0,0].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[0,0].set_ylim(-0.2,1.1)
    #Plot Ba vs Mg
    print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'])}")
    ax[0,1].errorbar(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[0,1].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,1].set_ylabel('[Ba/Fe]',fontsize=12)
    # ax[0,1].set_ylim(-0.2,1.1)
    # ax[0,1].set_xlim(-0.2,0.75)
    #Plot Mg vs iso_mg
    print(f"Mg vs iso_mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Mg_values['[Mg/Fe]'])}")
    ax[1,2].errorbar(Mg_values['[Mg/Fe]'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[1,2].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[1,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Eu vs iso_mg
    print(f"Eu vs iso_mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'])}")
    ax[1,0].errorbar(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[1,0].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[1,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[1,0].set_ylim(-0.2,1.1)
    #Plot Ba vs iso_mg
    print(f"Ba vs iso_mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'])}")
    ax[1,1].errorbar(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[1,1].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[1,1].set_ylabel('[Ba/Fe]',fontsize=12)
    # ax[1,1].set_ylim(-0.2,1.1)
    # ax[1,1].set_xlim(-0.2,0.75)
    #Plot Eu vs mg24
    print(f"Eu vs mg24: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_24']/100), Eu_values['[Eu/Fe]'])}")
    ax[2,0].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_24']/100), Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_MgFe24'], fmt='o')
    ax[2,0].set_xlabel('[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[2,0].set_ylim(-0.2,1.1)
    #Plot Ba vs mg24
    print(f"Ba vs mg24: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_24']/100), Ba_values['[Ba/Fe]'])}")
    ax[2,1].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_24']/100), Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_MgFe24'], fmt='o')
    ax[2,1].set_xlabel('[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,1].set_ylabel('[Ba/Fe]',fontsize=12)
    # ax[2,1].set_ylim(-0.2,1.1)
    # ax[2,1].set_xlim(-0.2,0.4)
    #Plot Eu vs mg25
    print(f"Eu vs mg25: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_25']/100), Eu_values['[Eu/Fe]'])}")
    ax[3,0].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_25']/100), Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_MgFe25'], fmt='o')
    ax[3,0].set_xlabel('[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[3,0].set_ylim(-0.2,1.1)
    #Plot Ba vs mg25
    print(f"Ba vs mg25: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_25']/100), Ba_values['[Ba/Fe]'])}")
    ax[3,1].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_MgFe25'], fmt='o')
    ax[3,1].set_xlabel('[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,1].set_ylabel('[Ba/Fe]',fontsize=12)
    # ax[3,1].set_ylim(-0.2,1.1)
    # ax[3,1].set_xlim(-0.05,0.10)
    #Plot Eu vs mg26
    print(f"Eu vs mg26: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Eu_values['[Eu/Fe]'])}")
    ax[4,0].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_MgFe26'], fmt='o')
    ax[4,0].set_xlabel('[$^{26}$Mg/Fe]',fontsize=12)
    ax[4,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[4,0].set_ylim(-0.2,1.1)
    #Plot Ba vs mg26
    print(f"Ba vs mg26: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Ba_values['[Ba/Fe]'])}")
    ax[4,1].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_MgFe26'], fmt='o')
    ax[4,1].set_xlabel('[$^{26}$Mg/Fe]',fontsize=12)
    ax[4,1].set_ylabel('[Ba/Fe]',fontsize=12)
    # ax[4,1].set_ylim(-0.2,1.1)
    # ax[4,1].set_xlim(-0.05,0.10)
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_24']/100), Mg_values['[Mg/Fe]'])}")
    ax[2,2].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_24']/100), Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_MgFe24'], fmt='o')
    ax[2,2].set_xlabel('[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Mg vs mg25
    print(f"Mg vs mg25: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_25']/100), Mg_values['[Mg/Fe]'])}")
    ax[3,2].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_25']/100), Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_MgFe25'], fmt='o')
    ax[3,2].set_xlabel('[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Mg vs mg26
    print(f"Mg vs mg26: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Mg_values['[Mg/Fe]'])}")
    ax[4,2].errorbar(Mg_values['[Mg/Fe]']*(isotope['R_26']/100), Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_MgFe26'], fmt='o')
    ax[4,2].set_xlabel('[$^{26}$Mg/Fe]',fontsize=12)
    ax[4,2].set_ylabel('[Mg/Fe]',fontsize=12)

    #plot difference plot for Mg vs iso mg
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Mg_values['[Mg/Fe]'] - Mg_values['[Mg/Fe]'])}")
    ax[0,2].errorbar(Mg_values['[Mg/Fe]'], Mg_values['[Mg/Fe]'] - Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o')
    ax[0,2].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/Fe] - [Mg/Fe]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_Fe_LblMg.png', dpi=300, bbox_inches='tight')


element_plots_XFe_LMg(star_list)

# %% isotope abundance vs wavelength

import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import ast

#Define the regions with thier wavelengths first
def region_wavelengths(region):
    if region == 1:    
        region = 5134.570
    '''Region 2'''
    if region == 2:
        region = 5138.710
    '''Region 3'''
    if region == 3:
        region = 5140.229
    '''Region 4'''
    if region == 4:
        region = 5134.208
    '''Region 5'''
    if region == 5:
        region = 5135.111
    '''Region 6'''
    if region == 6:
        region = 5135.999
    '''Region 7'''
    if region == 7:
        region = 5136.439
    '''Region 8'''
    if region == 8:
        region = 5138.486
    '''Region 9'''
    if region == 9:
        region = 5141.234
    '''Region 10'''
    if region == 10:
        region = 5133.174
    return region

def get_region(r):
    if r == 0:
        lw = 5134.42
        uw = 5140.46
    elif r == 1:
        lw = 5134.42
        uw = 5134.85
    elif r == 2:
        lw = 5138.55
        uw = 5138.95
    elif r == 3:
        lw = 5140.00
        uw = 5140.46
    elif r == 4:
        lw = 5134.0
        uw = 5134.4
    elif r == 5:
        lw = 5134.9
        uw = 5135.3
    elif r == 6:
        lw = 5135.9
        uw = 5136.3
    elif r == 7:
        lw = 5136.2
        uw = 5136.6
    elif r == 8:
        lw = 5138.2
        uw = 5138.6
    elif r == 9:
        lw = 5141.0
        uw = 5141.45
    elif r == 10:
        lw = 5133.0
        uw = 5133.4
    else:
        print('wavelength region error')
        lw = 0
        uw = 1
    return lw, uw


star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']

# star_list = ['hd_11695']
#Make an empty df to hold the isotope information
isotope_df = pd.DataFrame(columns=['star_name','s','mg','d_mg', 'i_24', 'i_25', 'i_26','R_24','R_25',
                                   'R_26','d_i_24', 'd_i_25', 'd_i_26','d_R_24','d_R_25','d_R_26',
                                   'pass','region','ratio'])

#Open the Masters_stars csv
star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')

Mg24_df = pd.DataFrame(columns=['star_name','R_24', 'd_R_24','region','wavelength'])
Mg25_df = pd.DataFrame(columns=['star_name','R_25', 'd_R_25','region','wavelength'])
Mg26_df = pd.DataFrame(columns=['star_name','R_26', 'd_R_26','region','wavelength'])

for star_name in star_list:
    # Read all the isotope abundance files
    # isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/par_unc_{star_name}.csv', delimiter=',', index_col=0)
    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/All_isotope_ratios_pre_avg.csv', delimiter=',')
    # /home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/par_unc_{star_name}.csv'
    #Get the star regions
    regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0] 
    # wavelength = np.array([5134.570, 5138.710, 5140.229, 5134.208, 5135.111, 5135.999, 5136.439, 5138.486, 5141.234, 5133.174]) 
    
    for region in regions:
        #Create a new df to add in the isotope and wavelength information for mg 24
        Mg24 = isotope[(isotope['star_name'] == star_name) & (isotope['region'] == region)][['star_name','R_24', 'd_R_24','region']]
        Mg25 = isotope[(isotope['star_name'] == star_name) & (isotope['region'] == region)][['star_name','R_25', 'd_R_25','region']]
        Mg26 = isotope[(isotope['star_name'] == star_name) & (isotope['region'] == region)][['star_name','R_26', 'd_R_26','region']]
        # Match the wavelength to the region
        wavelength = region_wavelengths(region)
        # print(wavelength)
        #Add the wavelength to the Mg24, Mg25, Mg26 df
        Mg24['wavelength'] = wavelength
        Mg25['wavelength'] = wavelength
        Mg26['wavelength'] = wavelength
        #Append the Mg24, Mg25, Mg26 df to the main df
        Mg24_df = pd.concat([Mg24_df, Mg24], ignore_index=True)
        Mg25_df = pd.concat([Mg25_df, Mg25], ignore_index=True)
        Mg26_df = pd.concat([Mg26_df, Mg26], ignore_index=True)
    #Take the average of the Mg24, Mg25, Mg26 values for each region
    Mg24_avg = Mg24_df.groupby('wavelength').mean().reset_index()
    Mg25_avg = Mg25_df.groupby('wavelength').mean().reset_index()
    Mg26_avg = Mg26_df.groupby('wavelength').mean().reset_index()
    # find the standard deviation of the Mg24, Mg25, Mg26 values for each region
    Mg24_std = Mg24_df.groupby('wavelength').std().reset_index()
    Mg25_std = Mg25_df.groupby('wavelength').std().reset_index()
    Mg26_std = Mg26_df.groupby('wavelength').std().reset_index()


        
wavelength = np.array([5134.570, 5138.710, 5140.229, 5134.208, 5135.111, 5135.999, 5136.439, 5138.486, 5141.234, 5133.174]) 
#Plot wavelength from 5134.570 to 5133.174 and each isotope datapoint with error bars in the region it was found
#Plot the isotope mg abundance with error bars
plt.figure(figsize=(8,4))
#plot Mg24
plt.errorbar(Mg24_avg['wavelength'], Mg24_avg['R_24'], yerr=Mg24_std['R_24'], fmt='o', color='#e41a1c', label='Mg24', capsize=3)
plt.plot(Mg24_avg['wavelength'], Mg24_avg['R_24'], color='#e41a1c')
#plot Mg25
plt.errorbar(Mg25_avg['wavelength'], Mg25_avg['R_25'], yerr=Mg25_std['R_25'], fmt='o', color='#377eb8', label='Mg25', capsize=3)
plt.plot(Mg25_avg['wavelength'], Mg25_avg['R_25'], color='#377eb8')
#Plot Mg26
plt.errorbar(Mg26_avg['wavelength'], Mg26_avg['R_26'], yerr=Mg26_std['R_26'], fmt='o', color='#4daf4a', label='Mg26', capsize=3) 
plt.plot(Mg26_avg['wavelength'], Mg26_avg['R_26'], color='#4daf4a')
plt.xlim(wavelength.min()-0.5, wavelength.max()+0.5)
plt.ylim(-10,100)
plt.xlabel('Wavelength ($\AA$)',fontsize=12)
plt.ylabel('Isotope Percentage',fontsize=12)
plt.legend(loc='upper right')
#Draw a box around each region using get_region
for region in range(1,11):
    lw, uw = get_region(region)
    plt.fill_between([lw, uw], 0.3, 120, facecolor = '#CCDBFD', alpha = 0.5)
    # put a label at the top of the plot with the region number
    plt.text((lw+uw)/2, -8, f'R{region}', fontsize=12, ha='center')

#Save the plot
plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Isotope_Percentage_vs_wavelength.png', dpi=300, bbox_inches='tight')


# %% X/H with shuffling of figures
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XH_new(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        
        # # Find the stars with the Mg/H values less than 0.25 and add then to a seperate variable
        # if Mg_values['[Mg/H]'] < 0.25:
        #     Mg_small = Mg_values          
            

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,20), constrained_layout=True)
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    
    # #Plot Eu vs Mg
    # print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'])}")
    # ax[0,0].errorbar(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5)
    # ax[0,0].errorbar(small_mg['mg'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    # ax[0,0].set_xlabel('[Mg/H]',fontsize=12)
    # ax[0,0].set_ylabel('[Eu/H]',fontsize=12)
    # ax[0,0].set_ylim(-0.7,0.7)
    # ax[0,0].set_xlim(-0.3,1)
    # #Plot Ba vs Mg
    # print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'])}")
    # ax[0,1].errorbar(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5)
    # ax[0,1].errorbar(small_mg['mg'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    # ax[0,1].set_xlabel('[Mg/H]',fontsize=12)
    # ax[0,1].set_ylabel('[Ba/H]',fontsize=12)
    # #Plot Mg vs iso_mg
    # print(f"Mg vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg'], Mg_values['[Mg/H]'])}")
    # ax[1,2].errorbar(iso_mg['mg'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    # ax[1,2].errorbar(small_mg['mg'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5,color='orange')
    # ax[1,2].set_xlabel('IS[Mg/H]',fontsize=12)
    # ax[1,2].set_ylabel('[Mg/H]',fontsize=12)
    # #Plot Eu vs iso_mg
    # print(f"Eu vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg'], Eu_values['[Eu/H]'])}")
    # ax[1,0].errorbar(iso_mg['mg'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    # ax[1,0].errorbar(small_mg['mg'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    # ax[1,0].set_xlabel('IS[Mg/H]',fontsize=12)
    # ax[1,0].set_ylabel('[Eu/H]',fontsize=12)
    # ax[1,0].set_ylim(-0.7,0.7)
    # ax[1,0].set_xlim(-0.2,1)
    # #Plot Ba vs iso_mg
    # print(f"Ba vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg'], Ba_values['[Ba/H]'])}")
    # ax[1,1].errorbar(iso_mg['mg'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    # ax[1,1].errorbar(small_mg['mg'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    # ax[1,1].set_xlabel('IS[Mg/H]',fontsize=12)
    # ax[1,1].set_ylabel('[Ba/H]',fontsize=12)
    
    
    # #Plot mg24 vs Eu
    # print(f"mg24 vs Eu: {scipy.stats.pearsonr(iso_mg['mg24'] ,Eu_values['[Eu/H]'])}")
    # ax[2,0].errorbar( iso_mg['mg24'], Eu_values['[Eu/H]'],yerr=iso_mg['d_mg24'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[2,0].errorbar(small_mg['mg24'], small_eu_values['[Eu/H]'],  yerr=small_mg['d_mg24'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[2,0].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    # ax[2,0].set_ylabel('[Eu/H]',fontsize=12)
    # ax[2,0].set_ylim(-0.5,0.6)
    # ax[2,0].set_xlim(-1.4,1.4)
    # #Plot mg25 vs Eu
    # print(f"mg25 vs Eu: {scipy.stats.pearsonr( Eu_values['[Eu/H]'],iso_mg['mg25'])}")
    # ax[3,0].errorbar(Eu_values['[Eu/H]'],iso_mg['mg25'],  yerr=iso_mg['d_mg25'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[3,0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[3,0].set_xlabel('[Eu/H]',fontsize=12)
    # ax[3,0].set_ylabel('IS[$^{25}$Mg/H]',fontsize=12)
    # ax[3,0].set_ylim(-0.75,0.9)
    # ax[3,0].set_xlim(-0.5,0.75)
    # #Plot mg26 vs Eu
    # print(f"mg26 vs Eu: {scipy.stats.pearsonr( Eu_values['[Eu/H]'],iso_mg['mg26'])}")
    # ax[4,0].errorbar( Eu_values['[Eu/H]'],iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[4,0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[4,0].set_xlabel('[Eu/H]',fontsize=12)
    # ax[4,0].set_ylabel('IS[$^{26}$Mg/H]',fontsize=12)
    # # ax[4,0].set_ylim(-0.1,1.6)
    # ax[4,0].set_xlim(-0.5,0.75)
    
    # #Plot mg24 vs Ba
    # print(f"mg24 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg24'])}")
    # ax[2,1].errorbar(iso_mg['mg24'], Ba_values['[Ba/H]'], yerr=iso_mg['d_mg24'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[2,1].errorbar(small_mg['mg24'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    # ax[2,1].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    # ax[2,1].set_ylabel('[Ba/H]',fontsize=12)
    # #Plot mg25 vs Ba
    # print(f"mg25 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg25'])}")
    # ax[3,1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[3,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[3,1].set_xlabel('[Ba/H]',fontsize=12)
    # ax[3,1].set_ylabel('IS[$^{25}$Mg/H]',fontsize=12)
    # #Plot mg26 vs Ba
    # print(f"mg26 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'],iso_mg['mg26'])}")
    # ax[4,1].errorbar(Ba_values['[Ba/H]'],iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[4,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[4,1].set_xlabel('[Ba/H]',fontsize=12)
    # ax[4,1].set_ylabel('IS[$^{26}$Mg/H]',fontsize=12)
    
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Mg_values['[Mg/H]'])}")
    ax[2,2].errorbar(iso_mg['mg24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5)
    ax[2,2].errorbar(small_mg['mg24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,2].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg24 vs mg25
    print(f"Mg24 vs mg25: {scipy.stats.pearsonr(iso_mg['mg24'], iso_mg['mg25'])}")
    ax[3,2].errorbar(iso_mg['mg24'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5)
    ax[3,2].errorbar(small_mg['mg24'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    ax[3,2].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[3,2].set_ylabel('IS[$^{25}$Mg/H]',fontsize=12)
    #Plot Mg24 vs mg26
    print(f"Mg24 vs mg26: {scipy.stats.pearsonr(iso_mg['mg24'], iso_mg['mg26'])}")
    ax[4,2].errorbar(iso_mg['mg24'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5)
    ax[4,2].errorbar(small_mg['mg24'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    ax[4,2].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[4,2].set_ylabel('IS[$^{26}$Mg/H]',fontsize=12)
    #plot difference plot for Mg vs iso mg
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(iso_mg['mg'], Mg_values['[Mg/H]'] - iso_mg['mg'])}")
    ax[0,2].errorbar(iso_mg['mg'], Mg_values['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]']-iso_mg['d_mg'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[0,2].errorbar(small_mg['mg'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    ax[0,2].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/H] - IS[Mg/H]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_new.png', dpi=300, bbox_inches='tight')



element_plots_XH_new(star_list)

# %% shuffling parameters teff, logg, feh
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def single_el_vs_params(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    Eu_values_H = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values_H = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values_H = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    
    #Open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #remove the 10th row
    star_info = star_info.drop(10)
    #reset index
    star_info = star_info.reset_index(drop=True)
    
    for star_name in star_list:
        #Extract the FEH value from the masters stars csv
        feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        #EUH
        eu_row_H = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row_H = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row_H = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        
        #Take away feh from each X/H value
        eu_row['[X/H]'] = eu_row['[X/H]'] - feh
        ba_row['[X/H]'] = ba_row['[X/H]'] - feh
        mg_row['[X/H]'] = mg_row['[X/H]'] - feh
        
        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)
        eu_row_H.insert(0, 'star_name', star_name)
        ba_row_H.insert(0, 'star_name', star_name)
        mg_row_H.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/Fe]', 'e[Eu/Fe]']
        ba_row.columns = ['star_name', '[Ba/Fe]', 'e[Ba/Fe]']
        mg_row.columns = ['star_name', '[Mg/Fe]', 'e[Mg/Fe]']
        eu_row_H.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row_H.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row_H.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        Eu_values_H = pd.concat([Eu_values_H, eu_row_H], ignore_index=True)
        Ba_values_H = pd.concat([Ba_values_H, ba_row_H], ignore_index=True)
        Mg_values_H = pd.concat([Mg_values_H, mg_row_H], ignore_index=True)
    
    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg_fe', 'd_mg_fe','mg','d_mg']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values_H[Mg_values_H['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values_H[Eu_values_H['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values_H[Ba_values_H['star_name'].isin(small_mg['Unnamed: 0'])]
    small_star_teff = star_info[star_info['ID2'].isin(small_mg['Unnamed: 0'])]
    small_star_logg = star_info[star_info['ID2'].isin(small_mg['Unnamed: 0'])]


    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(3, 4, figsize=(16,12),constrained_layout=True)

    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    # print(scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]']))
    ax[0,0].errorbar(Mg_values_H['[Mg/H]'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], xerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    ax[0,0].errorbar(small_mg_values['[Mg/H]'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[0,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,0].set_ylim(-0.4,0.9)
    #Plot Eu vs Ba
    print(f"Eu vs Ba: {scipy.stats.pearsonr(Ba_values_H['[Ba/H]'], Eu_values_H['[Eu/H]'])}")
    ax[0,1].errorbar(Ba_values_H['[Ba/H]'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], xerr=Ba_values_H['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[0,1].errorbar(small_ba_values['[Ba/H]'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[0,1].set_xlabel('[Ba/H]',fontsize=12)
    ax[0,1].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,1].set_ylim(-0.4,1)    
    #Plot Ba vs Mg
    ax[0,2].errorbar(Mg_values_H['[Mg/H]'], Ba_values_H['[Ba/H]'], yerr=Ba_values_H['e[Ba/H]'], xerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    ax[0,2].errorbar(small_mg_values['[Mg/H]'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[0,2].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,2].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs Teff
    print(f"Eu vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Eu_values_H['[Eu/H]'])}")
    ax[1,0].errorbar(star_info['TEFF'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], fmt='o', elinewidth=0.5)
    ax[1,0].errorbar(small_star_teff['TEFF'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,0].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[1,0].set_ylim(-0.4,0.9)
    #Plot Ba vs Teff
    print(f"Ba vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Ba_values_H['[Ba/H]'])}")
    ax[1,1].errorbar(star_info['TEFF'], Ba_values_H['[Ba/H]'], yerr=Ba_values_H['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[1,1].errorbar(small_star_teff['TEFF'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,1].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,1].set_ylabel('[Ba/H]',fontsize=12)
    
    
    #Plot Eu vs logg
    print(f"Eu vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Eu_values_H['[Eu/H]'])}")
    ax[2,0].errorbar(star_info['LOGG'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], fmt='o', elinewidth=0.5)
    ax[2,0].errorbar(small_star_logg['LOGG'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[2,0].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[2,0].set_ylim(-0.4,0.9)
    #Plot Ba vs logg
    print(f"Ba vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Ba_values_H['[Ba/H]'])}")
    ax[2,1].errorbar(star_info['LOGG'], Ba_values_H['[Ba/H]'], yerr=Ba_values_H['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[2,1].errorbar(small_star_logg['LOGG'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[2,1].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,1].set_ylabel('[Ba/H]',fontsize=12)
    
    #Plot Mg difference vs Teff
    print(f"Mg diff vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Mg_values_H['[Mg/H]'])}")
    ax[1,3].errorbar(star_info['TEFF'], Mg_values_H['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values_H['e[Mg/H]']-iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[1,3].errorbar(small_star_teff['TEFF'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,3].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,3].set_ylabel('[Mg/H] - IS[Mg/H]',fontsize=12)
    
    #Plot Mg differecne vs logg
    print(f"Mg diff vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Mg_values_H['[Mg/H]'] - iso_mg['mg_fe'])}")
    ax[2,3].errorbar(star_info['LOGG'], Mg_values_H['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values_H['e[Mg/H]']-iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[2,3].errorbar(small_star_logg['LOGG'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,3].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,3].set_ylabel('[Mg/H] - IS[Mg/H]',fontsize=12)
    
    # Plot Mg diffvs isoMg
    ax[0,3].errorbar(iso_mg['mg'],  Mg_values_H['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values_H['e[Mg/H]']-iso_mg['d_mg'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5)
    ax[0,3].errorbar(small_mg['mg'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]']-small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange')
    ax[0,3].set_xlabel('IS[Mg/H]',fontsize=12)
    ax[0,3].set_ylabel('[Mg/H]-IS[Mg/H]',fontsize=12)
    
    #Plot Mg vs Teff
    print(f"Mg vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Mg_values_H['[Mg/H]'])}")
    ax[1,2].errorbar(star_info['TEFF'], Mg_values_H['[Mg/H]'], yerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    ax[1,2].errorbar(small_star_teff['TEFF'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,2].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg vs logg
    print(f"Mg vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Mg_values_H['[Mg/H]'])}")
    ax[2,2].errorbar(star_info['LOGG'], Mg_values_H['[Mg/H]'], yerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    ax[2,2].errorbar(small_star_logg['LOGG'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[2,2].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,2].set_ylabel('[Mg/H]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_vH_params.png', dpi=300, bbox_inches='tight')



single_el_vs_params(star_list) 

# %% shuffling x/h lbl
#%% X/H with ratios made from lbl Mg
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XFe_LMG(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26','MgH','MgH24','MgH25','MgH26',
                      'd_MgH','d_MgH24','d_MgH25','d_MgH26']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,20), constrained_layout=True)
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    
    #Plot Eu vs Mg
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'])}")
    ax[0,0].errorbar(Mg_values['[Mg/H]'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5)
    ax[0,0].errorbar(small_mg_values['[Mg/H]'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg_values['e[Mg/H]'], fmt='o', color='orange', elinewidth=0.5)
    ax[0,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,0].set_ylim(-0.65,0.65)
    ax[0,0].set_xlim(-0.3,0.5)
    #Plot Ba vs Mg
    print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'])}")
    ax[0,1].errorbar(Mg_values['[Mg/H]'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=Mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5)
    ax[0,1].errorbar(small_mg_values['[Mg/H]'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg_values['e[Mg/H]'], fmt='o', color='orange', elinewidth=0.5)
    ax[0,1].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Mg vs iso_mg
    print(f"Mg vs iso_mg: {scipy.stats.pearsonr(iso_mg['MgH'], Mg_values['[Mg/H]'])}")
    ax[1,2].errorbar(iso_mg['MgH'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH'], fmt='o', elinewidth=0.5)
    ax[1,2].errorbar(small_mg['MgH'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_MgH'], fmt='o', color='orange', elinewidth=0.5)
    ax[1,2].set_xlabel('[Mg/H]',fontsize=12)
    ax[1,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Eu vs iso_mg
    print(f"Eu vs iso_mg: {scipy.stats.pearsonr(iso_mg['MgH'], Eu_values['[Eu/H]'])}")
    ax[1,0].errorbar(iso_mg['MgH'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH'], fmt='o', elinewidth=0.5)
    ax[1,0].errorbar(small_mg['MgH'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg['d_MgH'], fmt='o', color='orange', elinewidth=0.5)
    ax[1,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[1,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[1,0].set_ylim(-0.65,0.65)
    ax[1,0].set_xlim(-0.3,0.5)
    #Plot Ba vs iso_mg
    print(f"Ba vs iso_mg: {scipy.stats.pearsonr(iso_mg['MgH'], Ba_values['[Ba/H]'])}")
    ax[1,1].errorbar(iso_mg['MgH'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_MgH'], fmt='o', elinewidth=0.5)
    ax[1,1].errorbar(small_mg['MgH'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg['d_MgH'], fmt='o', color='orange', elinewidth=0.5)
    ax[1,1].set_xlabel('[Mg/H]',fontsize=12)
    ax[1,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg24
    print(f"Eu vs mg24: {scipy.stats.pearsonr(iso_mg['MgH24'], Eu_values['[Eu/H]'])}")
    ax[2,0].errorbar(iso_mg['MgH24'], Eu_values['[Eu/H]'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH24'], fmt='o', elinewidth=0.5)
    ax[2,0].errorbar(small_mg['MgH24'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg['d_MgH24'], fmt='o', color='orange', elinewidth=0.5)
    ax[2,0].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[2,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[2,0].set_ylim(-0.7,0.7)
    #Plot Ba vs mg24
    print(f"Ba vs mg24: {scipy.stats.pearsonr(iso_mg['MgH24'], Ba_values['[Ba/H]'])}")
    ax[2,1].errorbar(iso_mg['MgH24'], Ba_values['[Ba/H]'], yerr=Ba_values['e[Ba/H]'], xerr=iso_mg['d_MgH24'], fmt='o', elinewidth=0.5)
    ax[2,1].errorbar(small_mg['MgH24'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg['d_MgH24'], fmt='o', color='orange', elinewidth=0.5)
    ax[2,1].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[2,1].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs mg25
    print(f"Eu vs mg25: {scipy.stats.pearsonr( Eu_values['[Eu/H]'],iso_mg['MgH25'],)}")
    ax[3,0].errorbar(Eu_values['[Eu/H]'], iso_mg['MgH25'], yerr=Eu_values['e[Eu/H]'], xerr=iso_mg['d_MgH25'], fmt='o', elinewidth=0.5)
    ax[3,0].errorbar(small_eu_values['[Eu/H]'],small_mg['MgH25'], yerr=small_mg['d_MgH25'], xerr=small_eu_values['e[Eu/H]'], fmt='o', color='orange', elinewidth=0.5)
    ax[3,0].set_ylabel('[$^{25}$Mg/H]',fontsize=12)
    ax[3,0].set_xlabel('[Eu/H]',fontsize=12)
    ax[3,0].set_ylim(-0.05,0.11)
    ax[3,0].set_xlim(-0.5,0.6)
    
    #Plot Ba vs mg25
    print(f"Ba vs mg25: {scipy.stats.pearsonr(Ba_values['[Ba/H]'],iso_mg['MgH25'])}")
    ax[3,1].errorbar(Ba_values['[Ba/H]'],iso_mg['MgH25'],  yerr=iso_mg['d_MgH25'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[3,1].errorbar(small_ba_values['[Ba/H]'],small_mg['MgH25'], yerr=small_mg['d_MgH25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', color='orange', elinewidth=0.5)
    ax[3,1].set_ylabel('[$^{25}$Mg/H]',fontsize=12)
    ax[3,1].set_xlabel('[Ba/H]',fontsize=12)
    ax[3,1].set_ylim(-0.05,0.11)
    #Plot Eu vs mg26
    print(f"Eu vs mg26: {scipy.stats.pearsonr(Eu_values['[Eu/H]'],iso_mg['MgH26'])}")
    ax[4,0].errorbar( Eu_values['[Eu/H]'],iso_mg['MgH26'], yerr=iso_mg['d_MgH26'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    ax[4,0].errorbar(small_eu_values['[Eu/H]'],small_mg['MgH26'], yerr=small_mg['d_MgH26'], xerr=small_eu_values['e[Eu/H]'], fmt='o', color='orange', elinewidth=0.5)
    ax[4,0].set_ylabel('[$^{26}$Mg/H]',fontsize=12)
    ax[4,0].set_xlabel('[Eu/H]',fontsize=12)
    ax[4,0].set_ylim(-0.05,0.11)
    ax[4,0].set_xlim(-0.5,0.6)
    #Plot Ba vs mg26
    print(f"Ba vs mg26: {scipy.stats.pearsonr(Ba_values['[Ba/H]'],iso_mg['MgH26'])}")
    ax[4,1].errorbar( Ba_values['[Ba/H]'],iso_mg['MgH26'], yerr=iso_mg['d_MgH26'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[4,1].errorbar(small_ba_values['[Ba/H]'],small_mg['MgH26'], yerr=small_mg['d_MgH26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', color='orange', elinewidth=0.5)
    ax[4,1].set_ylabel('[$^{26}$Mg/H]',fontsize=12)
    ax[4,1].set_xlabel('[Ba/H]',fontsize=12)
    ax[4,1].set_ylim(-0.05,0.11)
    
    
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['MgH24'], Mg_values['[Mg/H]'])}")
    ax[2,2].errorbar(iso_mg['MgH24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH24'], fmt='o', elinewidth=0.5)
    ax[2,2].errorbar(small_mg['MgH24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_MgH24'], fmt='o', color='orange', elinewidth=0.5)
    ax[2,2].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/H]',fontsize=12)
    #Plot Mg vs mg25
    print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['MgH24'],iso_mg['MgH25'])}")
    ax[3,2].errorbar(iso_mg['MgH24'],iso_mg['MgH25'],  yerr=iso_mg['d_MgH25'], xerr=iso_mg['d_MgH24'], fmt='o', elinewidth=0.5)
    ax[3,2].errorbar(small_mg['MgH24'], small_mg['MgH25'], yerr=small_mg['d_MgH25'], xerr=small_mg['d_MgH24'], fmt='o', color='orange', elinewidth=0.5)
    ax[3,2].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[3,2].set_ylabel('[$^{25}$Mg/H]',fontsize=12)
    ax[3,2].set_ylim(-0.05,0.11)
    #Plot Mg vs mg26
    print(f"Mg vs mg26: {scipy.stats.pearsonr(iso_mg['MgH24'],iso_mg['MgH26'],)}")
    ax[4,2].errorbar(iso_mg['MgH24'],iso_mg['MgH26'], yerr=iso_mg['d_MgH26'], xerr=iso_mg['d_MgH24'], fmt='o', elinewidth=0.5)
    ax[4,2].errorbar(small_mg['MgH24'], small_mg['MgH26'], yerr=small_mg['d_MgH26'], xerr=small_mg['d_MgH24'], fmt='o', color='orange', elinewidth=0.5)
    ax[4,2].set_xlabel('[$^{24}$Mg/H]',fontsize=12)
    ax[4,2].set_ylabel('[$^{26}$Mg/H]',fontsize=12)
    ax[4,2].set_ylim(-0.05,0.11)
    #plot difference plot for Mg vs iso mg
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(iso_mg['MgH'], Mg_values['[Mg/H]'] - iso_mg['MgH'])}")
    ax[0,2].errorbar(iso_mg['MgH'], Mg_values['[Mg/H]'] - iso_mg['MgH'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_MgH'], fmt='o', elinewidth=0.5)
    ax[0,2].errorbar(small_mg['MgH'], small_mg_values['[Mg/H]'] - small_mg['MgH'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_MgH'], fmt='o', color='orange', elinewidth=0.5)
    ax[0,2].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/H] - [Mg/H]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_LblMg_new.png', dpi=300, bbox_inches='tight')

element_plots_XFe_LMG(star_list)

# %% X/Fe with shuffling of figures
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XH_new(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    #Open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    
    for star_name in star_list:
        #Extract the FEH value from the masters stars csv
        feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        #Take away feh from each X/H value
        eu_row['[X/H]'] = eu_row['[X/H]'] - feh
        ba_row['[X/H]'] = ba_row['[X/H]'] - feh
        mg_row['[X/H]'] = mg_row['[X/H]'] - feh

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/Fe]', 'e[Eu/Fe]']
        ba_row.columns = ['star_name', '[Ba/Fe]', 'e[Ba/Fe]']
        mg_row.columns = ['star_name', '[Mg/Fe]', 'e[Mg/Fe]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        
        # # Find the stars with the Mg/H values less than 0.25 and add then to a seperate variable
        # if Mg_values['[Mg/H]'] < 0.25:
        #     Mg_small = Mg_values          
            

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    # iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    iso_mg = isotope[['Unnamed: 0','mg_fe', 'd_mg_fe', 'mg_fe24', 'd_mg_fe24', 'mg_fe25', 'd_mg_fe25', 'mg_fe26', 'd_mg_fe26','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(5, 3, figsize=(12,20), constrained_layout=True)
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    
    #Plot Eu vs Mg
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'])}")
    ax[0,0].errorbar(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o', elinewidth=0.5)
    ax[0,0].errorbar(small_mg['mg_fe'], small_eu_values['[Eu/Fe]'], yerr=small_eu_values['e[Eu/Fe]'], xerr=small_mg['d_mg_fe'], fmt='o', elinewidth=0.5, color='orange')
    ax[0,0].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[0,0].set_ylim(-0.7,0.7)
    ax[0,0].set_xlim(-0.3,1)
    #Plot Ba vs Mg
    print(f"Ba vs Mg: {scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'])}")
    ax[0,1].errorbar(Mg_values['[Mg/Fe]'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=Mg_values['e[Mg/Fe]'], fmt='o', elinewidth=0.5)
    ax[0,1].errorbar(small_mg['mg_fe'], small_ba_values['[Ba/Fe]'], yerr=small_ba_values['e[Ba/Fe]'], xerr=small_mg['d_mg_fe'], fmt='o', elinewidth=0.5, color='orange')
    ax[0,1].set_xlabel('[Mg/Fe]',fontsize=12)
    ax[0,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot Mg vs iso_mg
    print(f"Mg vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'])}")
    ax[1,2].errorbar(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o', elinewidth=0.5)
    ax[1,2].errorbar(small_mg['mg_fe'], small_mg_values['[Mg/Fe]'], yerr=small_mg_values['e[Mg/Fe]'], xerr=small_mg['d_mg_fe'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,2].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[1,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Eu vs iso_mg
    print(f"Eu vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Eu_values['[Eu/Fe]'])}")
    ax[1,0].errorbar(iso_mg['mg_fe'], Eu_values['[Eu/Fe]'], yerr=Eu_values['e[Eu/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o', elinewidth=0.5)
    # ax[1,0].errorbar(small_mg['mg_fe'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg['d_mg_fe'], fmt='o', elinewidth=0.5, color='orange')
    ax[1,0].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[1,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[1,0].set_ylim(-0.7,0.7)
    ax[1,0].set_xlim(-0.2,1)
    #Plot Ba vs iso_mg
    print(f"Ba vs iso_mg: {scipy.stats.pearsonr(iso_mg['mg_fe'], Ba_values['[Ba/Fe]'])}")
    ax[1,1].errorbar(iso_mg['mg_fe'], Ba_values['[Ba/Fe]'], yerr=Ba_values['e[Ba/Fe]'], xerr=iso_mg['d_mg_fe'], fmt='o', elinewidth=0.5)
    ax[1,1].errorbar(small_mg['mg_fe'], small_ba_values['[Ba/Fe]'], yerr=small_ba_values['e[Ba/Fe]'], xerr=small_mg['d_mg_fe'], fmt='o', elinewidth=0.5, color='orange')
    ax[1,1].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[1,1].set_ylabel('[Ba/Fe]',fontsize=12)
    
    
    #Plot mg24 vs Eu
    print(f"mg24 vs Eu: {scipy.stats.pearsonr(iso_mg['mg_fe'] ,Eu_values['[Eu/Fe]'])}")
    ax[2,0].errorbar( iso_mg['mg_fe24'], Eu_values['[Eu/Fe]'],yerr=iso_mg['d_mg_fe24'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    ax[2,0].errorbar(small_mg['mg_fe24'], small_eu_values['[Eu/Fe]'],  yerr=small_mg['d_mg_fe24'], xerr=small_eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,0].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,0].set_ylabel('[Eu/Fe]',fontsize=12)
    ax[2,0].set_ylim(-0.5,0.6)
    ax[2,0].set_xlim(-1.4,1.4)
    #Plot mg25 vs Eu
    print(f"mg25 vs Eu: {scipy.stats.pearsonr( Eu_values['[Eu/Fe]'],iso_mg['mg_fe25'])}")
    ax[3,0].errorbar(Eu_values['[Eu/Fe]'],iso_mg['mg_fe25'],  yerr=iso_mg['d_mg_fe25'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    ax[3,0].errorbar(small_eu_values['[Eu/Fe]'], small_mg['mg_fe25'], yerr=small_mg['d_mg_fe25'], xerr=small_eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3,0].set_xlabel('[Eu/Fe]',fontsize=12)
    ax[3,0].set_ylabel('IS[$^{25}$Mg/Fe]',fontsize=12)
    ax[3,0].set_ylim(-0.75,0.9)
    ax[3,0].set_xlim(-0.5,0.75)
    #Plot mg26 vs Eu
    print(f"mg26 vs Eu: {scipy.stats.pearsonr( Eu_values['[Eu/Fe]'],iso_mg['mg_fe26'])}")
    ax[4,0].errorbar( Eu_values['[Eu/Fe]'],iso_mg['mg_fe26'], yerr=iso_mg['d_mg_fe26'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    ax[4,0].errorbar(small_eu_values['[Eu/Fe]'], small_mg['mg_fe26'], yerr=small_mg['d_mg_fe26'], xerr=small_eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5, color='orange')
    ax[4,0].set_xlabel('[Eu/Fe]',fontsize=12)
    ax[4,0].set_ylabel('IS[$^{26}$Mg/Fe]',fontsize=12)
    # ax[4,0].set_ylim(-0.1,1.6)
    ax[4,0].set_xlim(-0.5,0.75)
    
    #Plot mg24 vs Ba
    print(f"mg24 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], iso_mg['mg_fe24'])}")
    ax[2,1].errorbar(iso_mg['mg_fe24'], Ba_values['[Ba/Fe]'], yerr=iso_mg['d_mg_fe24'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[2,1].errorbar(small_mg['mg_fe24'], small_ba_values['[Ba/Fe]'], yerr=small_ba_values['e[Ba/Fe]'], xerr=small_mg['d_mg_fe24'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,1].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,1].set_ylabel('[Ba/Fe]',fontsize=12)
    #Plot mg25 vs Ba
    print(f"mg25 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], iso_mg['mg_fe25'])}")
    ax[3,1].errorbar(Ba_values['[Ba/Fe]'], iso_mg['mg_fe25'], yerr=iso_mg['d_mg_fe25'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[3,1].errorbar(small_ba_values['[Ba/Fe]'], small_mg['mg_fe25'], yerr=small_mg['d_mg_fe25'], xerr=small_ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3,1].set_xlabel('[Ba/Fe]',fontsize=12)
    ax[3,1].set_ylabel('IS[$^{25}$Mg/Fe]',fontsize=12)
    #Plot mg26 vs Ba
    print(f"mg26 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'],iso_mg['mg_fe26'])}")
    ax[4,1].errorbar(Ba_values['[Ba/Fe]'],iso_mg['mg_fe26'], yerr=iso_mg['d_mg_fe26'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[4,1].errorbar(small_ba_values['[Ba/Fe]'], small_mg['mg_fe26'], yerr=small_mg['d_mg_fe26'], xerr=small_ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5, color='orange')
    ax[4,1].set_xlabel('[Ba/Fe]',fontsize=12)
    ax[4,1].set_ylabel('IS[$^{26}$Mg/Fe]',fontsize=12)
    
    #Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg_fe24'], Mg_values['[Mg/Fe]'])}")
    ax[2,2].errorbar(iso_mg['mg_fe24'], Mg_values['[Mg/Fe]'], yerr=Mg_values['e[Mg/Fe]'], xerr=iso_mg['d_mg_fe24'], fmt='o', elinewidth=0.5)
    ax[2,2].errorbar(small_mg['mg_fe24'], small_mg_values['[Mg/Fe]'], yerr=small_mg_values['e[Mg/Fe]'], xerr=small_mg['d_mg_fe24'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,2].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[2,2].set_ylabel('[Mg/Fe]',fontsize=12)
    #Plot Mg24 vs mg25
    print(f"Mg24 vs mg25: {scipy.stats.pearsonr(iso_mg['mg_fe24'], iso_mg['mg_fe25'])}")
    ax[3,2].errorbar(iso_mg['mg_fe24'], iso_mg['mg_fe25'], yerr=iso_mg['d_mg_fe25'], xerr=iso_mg['d_mg_fe24'], fmt='o', elinewidth=0.5)
    ax[3,2].errorbar(small_mg['mg_fe24'], small_mg['mg_fe25'], yerr=small_mg['d_mg_fe25'], xerr=small_mg['d_mg_fe24'], fmt='o', elinewidth=0.5, color='orange')
    ax[3,2].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[3,2].set_ylabel('IS[$^{25}$Mg/Fe]',fontsize=12)
    #Plot Mg24 vs mg26
    print(f"Mg24 vs mg26: {scipy.stats.pearsonr(iso_mg['mg_fe24'], iso_mg['mg_fe26'])}")
    ax[4,2].errorbar(iso_mg['mg_fe24'], iso_mg['mg_fe26'], yerr=iso_mg['d_mg_fe26'], xerr=iso_mg['d_mg_fe24'], fmt='o', elinewidth=0.5)
    ax[4,2].errorbar(small_mg['mg_fe24'], small_mg['mg_fe26'], yerr=small_mg['d_mg_fe26'], xerr=small_mg['d_mg_fe24'], fmt='o', elinewidth=0.5, color='orange')
    ax[4,2].set_xlabel('IS[$^{24}$Mg/Fe]',fontsize=12)
    ax[4,2].set_ylabel('IS[$^{26}$Mg/Fe]',fontsize=12)
    #plot difference plot for Mg vs iso mg_fe
    print(f"Mg vs diff iso mg_fe: {scipy.stats.pearsonr(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'] - iso_mg['mg_fe'])}")
    ax[0,2].errorbar(iso_mg['mg_fe'], Mg_values['[Mg/Fe]'] - iso_mg['mg_fe'], yerr=Mg_values['e[Mg/Fe]']-iso_mg['d_mg_fe'], xerr=iso_mg['d_mg_fe'], fmt='o', elinewidth=0.5)
    ax[0,2].errorbar(small_mg['mg_fe'], small_mg_values['[Mg/Fe]'] - small_mg['mg_fe'], yerr=small_mg_values['e[Mg/Fe]']-small_mg['d_mg_fe'], xerr=small_mg['d_mg_fe'], fmt='o', elinewidth=0.5, color='orange')
    ax[0,2].set_xlabel('IS[Mg/Fe]',fontsize=12)
    ax[0,2].set_ylabel('[Mg/Fe] - IS[Mg/Fe]',fontsize=12)
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_Fe_new.png', dpi=300, bbox_inches='tight')



element_plots_XH_new(star_list)

# %% X/H with mg only
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XH_new_small(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    #Open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #remove the 10th row
    star_info = star_info.drop(10)
    #reset index
    star_info = star_info.reset_index(drop=True)
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        
        # # Find the stars with the Mg/H values less than 0.25 and add then to a seperate variable
        # if Mg_values['[Mg/H]'] < 0.25:
        #     Mg_small = Mg_values          
            

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    # add hd_100407 to the small_mg variable
    small_mg = small_mg.append(iso_mg[iso_mg['Unnamed: 0'] == 'hd_100407'])
    # print(small_mg)
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_star_teff = star_info[star_info['ID2'].isin(small_mg['Unnamed: 0'])]
    small_star_logg = star_info[star_info['ID2'].isin(small_mg['Unnamed: 0'])]

    #Calculate and print the median spread of small_mg_values
    print(f"Median spread of Orange values: {np.median(small_mg_values['[Mg/H]'])}")
    #print the median spread of the Mg_values
    print(f"Median spread of: {np.median(Mg_values['[Mg/H]'])}")
    # print the MAd of the small_mg_values
    print(f"MAD of Orange values: {np.median(np.abs(small_mg_values['[Mg/H]'] - np.median(small_mg_values['[Mg/H]'])))}")
    #print the median diffeence between the small_mg_values and the Mg_values
    # print(f"Median difference between Orange and Blue values: {np.median(small_mg_values['e[Mg/H]']) - Mg_values['e[Mg/H]']}")
    print(f"MAD of values: {np.median(np.abs(Mg_values['[Mg/H]'] - np.median(Mg_values['[Mg/H]'])))}")
    
    #Make the plot for X/H and mg from IS.
    fig, ax = plt.subplots(2, 4, figsize=(12,6), constrained_layout=True)
    
    # Plot Mg vs iso mg difference
    print(f"Mg vs diff iso mg: {scipy.stats.pearsonr(iso_mg['mg'], Mg_values['[Mg/H]'] - iso_mg['mg'])}")
    ax[0, 0].errorbar(iso_mg['mg'], Mg_values['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]'] - iso_mg['d_mg'], xerr=iso_mg['d_mg'], fmt='o', elinewidth=0.5,markersize=5)
    ax[0, 0].errorbar(small_mg['mg'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]'] - small_mg['d_mg'], xerr=small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    #draw a dotted line along at 0 on x axis
    ax[0, 0].axhline(0, color='black', linestyle='--', linewidth=0.7)
    ax[0, 0].set_xlabel('IS[Mg/H]', fontsize=12)
    ax[0, 0].set_ylabel('[Mg/H] - IS[Mg/H]', fontsize=12)

    # Plot Mg difference vs Teff
    print(f"Mg diff vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Mg_values['[Mg/H]'])}")
    ax[0, 1].errorbar(star_info['TEFF'], Mg_values['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]'] - iso_mg['d_mg'], fmt='o', elinewidth=0.5,markersize=5)
    ax[0, 1].errorbar(small_star_teff['TEFF'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]'] - small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    ax[0,1].axhline(0, color='black', linestyle='--', linewidth=0.7)
    ax[0, 1].set_xlabel('$T_{eff}$', fontsize=12)
    ax[0, 1].set_ylabel('[Mg/H] - IS[Mg/H]', fontsize=12)

    # Plot Mg difference vs logg
    print(f"Mg diff vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Mg_values['[Mg/H]'] - iso_mg['mg'])}")
    ax[0, 2].errorbar(star_info['LOGG'], Mg_values['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]'] - iso_mg['d_mg'], fmt='o', elinewidth=0.5,markersize=5)
    ax[0, 2].errorbar(small_star_logg['LOGG'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]'] - small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    ax[0,2].axhline(0, color='black', linestyle='--', linewidth=0.7)
    ax[0, 2].set_xlabel('$\log g (cm/s^{2}$)', fontsize=12)
    ax[0, 2].set_ylabel('[Mg/H] - IS[Mg/H]', fontsize=12)
    
    #plot Mg difference vs [Fe/H]
    print(f"Mg diff vs [Fe/H]: {scipy.stats.pearsonr(star_info['FEH'], Mg_values['[Mg/H]'] - iso_mg['mg'])}")
    ax[0, 3].errorbar(star_info['FEH'], Mg_values['[Mg/H]'] - iso_mg['mg'], yerr=Mg_values['e[Mg/H]'] - iso_mg['d_mg'], fmt='o', elinewidth=0.5,markersize=5)
    ax[0, 3].errorbar(small_star_logg['FEH'], small_mg_values['[Mg/H]'] - small_mg['mg'], yerr=small_mg_values['e[Mg/H]'] - small_mg['d_mg'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    ax[0,3].axhline(0, color='black', linestyle='--', linewidth=0.7)
    ax[0, 3].set_xlabel('[Fe/H]', fontsize=12)
    ax[0, 3].set_ylabel('[Mg/H] - IS[Mg/H]', fontsize=12)
    

    # Plot Mg vs mg24
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Mg_values['[Mg/H]'])}")
    ax[1, 0].errorbar(iso_mg['mg24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5,markersize=5)
    ax[1, 0].errorbar(small_mg['mg24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    ax[1, 0].set_xlabel('IS[$^{24}$Mg/H]', fontsize=12)
    ax[1, 0].set_ylabel('[Mg/H]', fontsize=12)
    ax[1, 0].set_xlim(-1.3, 1.25)

    # Plot Mg24 vs mg25
    print(f"Mg24 vs mg25: {scipy.stats.pearsonr(iso_mg['mg24'], iso_mg['mg25'])}")
    ax[1,1].errorbar(iso_mg['mg24'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5,markersize=5)
    ax[1,1].errorbar(small_mg['mg24'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    ax[1,1].set_xlabel('IS[$^{24}$Mg/H]', fontsize=12)
    ax[1,1].set_ylabel('IS[$^{25}$Mg/H]', fontsize=12)
    ax[1,1].set_xlim(-1.3, 1.25)

    # Plot Mg24 vs mg26
    print(f"Mg24 vs mg26: {scipy.stats.pearsonr(iso_mg['mg24'], iso_mg['mg26'])}")
    ax[1,2].errorbar(iso_mg['mg24'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5,markersize=5)
    ax[1,2].errorbar(small_mg['mg24'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange',markersize=5)
    ax[1,2].set_xlabel('IS[$^{24}$Mg/H]', fontsize=12)
    ax[1,2].set_ylabel('IS[$^{26}$Mg/H]', fontsize=12)
    ax[1,2].set_xlim(-1.3, 1.25)
    
    #plot an Hr diagram
    ax[1,3].errorbar(star_info['TEFF'], star_info['LOGG'], label='All Stars', fmt='o',markersize=5)
    ax[1,3].errorbar(small_star_teff['TEFF'], small_star_logg['LOGG'], c='orange', label='Mg < 0.25', fmt='o',markersize=5)
    ax[1,3].set_ylabel('$\log g (cm/s^{2}$)', fontsize=12)
    ax[1,3].set_xlabel('$T_{eff}$', fontsize=12)
    ax[1,3].invert_yaxis()
    ax[1,3].invert_xaxis()
    # Remove the last subplot
    # fig.delaxes(ax[1, 3])
    
    

    # Save the plot
    # plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_small.png', dpi=300, bbox_inches='tight')


#Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_small_rotated.png', dpi=300, bbox_inches='tight')

element_plots_XH_new_small(star_list)

#%%

import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np

star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

def element_plots_XH_new(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        
        # # Find the stars with the Mg/H values less than 0.25 and add then to a seperate variable
        # if Mg_values['[Mg/H]'] < 0.25:
        #     Mg_small = Mg_values          
            

    #Open the isotope mg abundance file
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(3, 2, figsize=(9, 12), constrained_layout=True)
    # increase the space between the plots
    # fig.subplots_adjust(hspace=0.35, wspace=0.35)

    # Plot mg24 vs Eu
    print(f"mg24 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg24'])}")
    ax[0, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    ax[0, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[0, 0].set_ylabel('IS[$^{24}$Mg/H]', fontsize=12)
    ax[0, 0].set_xlabel('[Eu/H]', fontsize=12)
    # ax[0, 0].set_ylim(-0.75, 0.9)
    ax[0, 0].set_xlim(-0.5, 0.75)

    # Plot mg25 vs Eu
    print(f"mg25 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg25'])}")
    ax[1, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    ax[1, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[1, 0].set_ylabel('IS[$^{25}$Mg/H]', fontsize=12)
    ax[1, 0].set_ylim(-0.75, 0.9)
    ax[1, 0].set_xlim(-0.5, 0.75)

    # Plot mg26 vs Eu
    print(f"mg26 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg26'])}")
    ax[2, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    ax[2, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[2, 0].set_ylabel('IS[$^{26}$Mg/H]', fontsize=12)
    ax[2, 0].set_xlim(-0.5, 0.75)

    # Plot mg24 vs Ba
    print(f"mg24 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg24'])}")
    ax[0, 1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[0, 1].set_ylabel('IS[$^{24}$Mg/H]', fontsize=12)
    ax[0, 1].set_xlabel('[Ba/H]', fontsize=12)

    # Plot mg25 vs Ba
    print(f"mg25 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg25'])}")
    ax[1,1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[1,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1,1].set_xlabel('[Ba/H]', fontsize=12)
    ax[1,1].set_ylabel('IS[$^{25}$Mg/H]', fontsize=12)

    # Plot mg26 vs Ba
    print(f"mg26 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg26'])}")
    ax[2, 1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    ax[2, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2, 1].set_xlabel('[Ba/H]', fontsize=12)
    ax[2, 1].set_ylabel('IS[$^{26}$Mg/H]', fontsize=12)

    # Save the plot
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_new.png', dpi=300, bbox_inches='tight')
    
    #Save the plot
    # plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Element_fits_X_H_new.png', dpi=300, bbox_inches='tight')



element_plots_XH_new(star_list)

# %%
