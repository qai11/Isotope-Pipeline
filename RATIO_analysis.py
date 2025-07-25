"""
Title: RATIO_analysis.py
Author: Quin Aicken Davies
Date: 12/06/2025

Description: Run adter Ratio_uncertainties to calculate the final abundances into csv forms.
"""
# %% """Make a table for the isotopic ratios"""
# """Makes the All_isotope_ratios_pre_avg file"""

vpass = '6_2'

import pandas as pd
import numpy as np
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# star_list = ['moon','hd_18907']
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407']
star_list = ['hd_10700']

#Make an empty df to hold the isotope information
isotope_df = pd.DataFrame(columns=['star_name','s','mg','d_mg', 'i_24', 'i_25', 'i_26','R_24','R_25',
                                   'R_26','d_i_24', 'd_i_25', 'd_i_26','d_R_24','d_R_25','d_R_26',
                                   'pass','region','ratio'])
for star_name in star_list:
    # Read all the abundance data for the star (across all regions)
    iso_abund_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/par_unc_{star_name}_paper_vpass_{vpass}.csv', delimiter=',', index_col=0)
   
    # Add the information to the isotope_df
    for i in range(len(iso_abund_params)):
        new_row = {
            'star_name': star_name,
            's': iso_abund_params['s'][i],
            'd_s': iso_abund_params['d_s'][i],
            'mg': iso_abund_params['mg'][i],
            # 'd_mg': iso_abund_params['d_mg'][i],
            'd_mg': 1,
            'i_24': iso_abund_params['i_24'][i],
            'd_i_24': iso_abund_params['d_i_24'][i],
            'i_25': iso_abund_params['i_25'][i],
            'd_i_25': iso_abund_params['d_i_25'][i],
            'i_26': iso_abund_params['i_26'][i],
            'd_i_26': iso_abund_params['d_i_26'][i],
            'R_24': iso_abund_params['R_24'][i],
            'd_R_24': iso_abund_params['d_R_24'][i],
            'R_25': iso_abund_params['R_25'][i],
            'd_R_25': iso_abund_params['d_R_25'][i],
            'R_26': iso_abund_params['R_26'][i],
            'd_R_26': iso_abund_params['d_R_26'][i],
            'pass': iso_abund_params['pass'][i],
            'region': iso_abund_params['region'][i],
            'ratio': f'{iso_abund_params["R_24"][i]}$\pm${round(iso_abund_params["d_R_24"][i],2)}:{iso_abund_params["R_25"][i]}$\pm${round(iso_abund_params["d_R_25"][i],2)}:{iso_abund_params["R_26"][i]}$\pm${round(iso_abund_params["d_R_26"][i],2)}'
        }
        isotope_df = isotope_df.append(pd.Series(new_row), ignore_index=True)
        
#save the isotope_df to a csv file
isotope_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/All_isotope_ratios_pre_avg_paper_vpass_{vpass}.csv')




# %%"""Caclulate the weighted averages of isotopic abundance ratios for each star"""
#Makes the weighted_avg_iso file

import pandas as pd
import numpy as np
import ast
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# star_list = ['moon','hd_18907']
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
# star_list = ['hd_11695']
# star_list = ['hd_10700']
# star_list = ['hd_23249','hd_128621','hd_10700']
# star_list = ['hd_18884']
# Initialize an empty dictionary to hold abundance data for each star
abundance_dict = {}

vpass = '6'

def calc_ratio(i_24, i_25, i_26):
    i24_percentage=1/(0.01*i_24)
    i25_percentage=1/(0.01*i_25)
    i26_percentage=1/(0.01*i_26)

    isotope_sum = i24_percentage + i25_percentage + i26_percentage
    print(f"sum {isotope_sum}")

    i24_ratio = (i24_percentage/isotope_sum) * 100
    i25_ratio = (i25_percentage/isotope_sum) * 100
    i26_ratio = (i26_percentage/isotope_sum) * 100

    return round(i24_ratio,2), round(i25_ratio,2), round(i26_ratio,2)
    
for star_name in star_list:
    # Read all the abundance data for the star (across all regions)
    # iso_abund_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/par_unc_{star_name}_paper_vpass_{vpass}.csv', delimiter=',', index_col=0)
    # print(iso_abund_params)
    iso_abund_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/par_unc_{star_name}_paper_vpass_{vpass}.csv', delimiter=',', index_col=0)
    
    #open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #get the star regions
    regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0]
    
    #only take the average of the given regions
    iso_abund_params = iso_abund_params[iso_abund_params['region'].isin(regions)]
    # Separate the error and parameter columns
    errors = iso_abund_params.filter(like='d_')
    params = iso_abund_params.drop(columns=errors.columns).drop(columns=['region', 'pass'])
    print(params)
    print(errors)
    # Calculate the weighted average of isotopic abundance ratios for the entire star
    weighted_avg = np.average(params, axis=0, weights=1/(errors**2))
    
    
    # # Re Normalize so ratios sum to 100%
    # ratio_columns = [-1, -2, -3]  # Columns that should sum to 100%
    # ratio_sum = weighted_avg[ratio_columns].sum()
    # if ratio_sum != 0:
    #     normalization_factor = 100 / ratio_sum
    #     weighted_avg[ratio_columns] *= normalization_factor  # Normalize abundances

    # Calculate propagated uncertainties
    summed_inv_error_sq = np.nansum(1 / (errors**2), axis=0)  
    weighted_avg_error = np.sqrt(1 / summed_inv_error_sq)
    
    # # **Normalize errors using the same factor applied to abundances**
    # if ratio_sum != 0:
    #     weighted_avg_error[ratio_columns] *= normalization_factor
    
    # Add the weighted average and error to the dictionary
    abundance_dict[star_name] = {
        'abundance': weighted_avg,
        'error': weighted_avg_error
    }
    
    

# Convert the abundance dictionary into a DataFrame
abundance_df = pd.DataFrame(abundance_dict).T


# Create a structured DataFrame to hold the final results
final_df = pd.DataFrame(columns=['s','i_24', 'i_25','i_26','R_24','R_25','R_26','d_s', 
                                 'd_i_24', 'd_i_25','d_i_26','d_R_24','d_R_25','d_R_26','mg24','mg25','mg26'
                                 ,'d_mg24','d_mg25','d_mg26','MgH','d_MgH','MgH24','MgH25','MgH26','d_MgH24','d_MgH25','d_MgH26',
                                 'MgFe','d_MgFe','MgFe24','MgFe25','MgFe26','d_MgFe24','d_MgFe25','d_MgFe26'])
#Open masters stars csv
star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
#remove the 10th row
star_info = star_info.drop(10)
#reset index
star_info = star_info.reset_index(drop=True)

# Iterate over each star and add the results to the structured DataFrame
for star_name in star_list:
    feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
    # Extract the abundance and error data for the current star
    abundances = abundance_df.loc[star_name, 'abundance']
    
    errors = abundance_df.loc[star_name, 'error']
    #Open summary abundances file
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    summary_abundances_Fe = summary_abundances.sort_values(by=['[X/Fe]', 'e[X/Fe]'], ascending=False)
    summary_abundances_Fe = summary_abundances.drop_duplicates(subset=['element'], keep='first')
    summary_abundances_H = summary_abundances.sort_values(by=['[X/H]', 'e[X/H]'], ascending=False)
    summary_abundances_H = summary_abundances.drop_duplicates(subset=['element'], keep='first')
    #Extract the Mg [X/H] and error
    MgFe = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/Fe]','e[X/Fe]']]
    MgH = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    print(MgH)
    # MgH['[X/H]'] = MgH['[X/H]'] - feh
    # print(MgH)
    # print(MgH)
    mg_24,mg_25,mg_26 = calc_ratio(abundances[1], abundances[2], abundances[3])
    # print(mg_24)
    # print(abundances[4])
    # print(errors[4], errors[5], errors[6])
    
    #Open summary abundances file for Mg abundance
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    #Extract the Mg [X/H] and error
    Mg_lbl_abundance = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    Mg_lbl_abundance = Mg_lbl_abundance['[X/H]'].values[0]
    
    #calculate the log abundances
    log_solar_mg = Mg_lbl_abundance+7.53
    log_solar_i_24 = log_solar_mg * (mg_24/100)
    log_solar_i_25 = log_solar_mg * (mg_25/100)
    log_solar_i_26 = log_solar_mg * (mg_26/100)
    #convert to Mg/H
    mg_24_H_is = log_solar_i_24 - (7.53*0.7899)
    print(mg_24_H_is)
    mg_25_H_is = log_solar_i_25 - (7.53*0.1000)
    mg_26_H_is = log_solar_i_26 - (7.53*0.1101) 
    
    #Do the same for the errors
    #calculate the log abundances
    log_solar_mg = errors[1]+7.53
    log_solar_i_24 = log_solar_mg * (mg_24/100)
    log_solar_i_25 = log_solar_mg * (mg_25/100)
    log_solar_i_26 = log_solar_mg * (mg_26/100)
    #convert to Mg/H
    mg_24_H_is_err = log_solar_i_24 - (7.53*0.7899)
    mg_25_H_is_err = log_solar_i_25 - (7.53*0.1000)
    mg_26_H_is_err = log_solar_i_26 - (7.53*0.1101)
    
    #line by line
    #Convert the line by line anundances
    log_solar_mg_lbl = MgH['[X/H]'].values[0] + 7.53
    log_solar_i_24_lbl = log_solar_mg_lbl * (mg_24/100)
    log_solar_i_25_lbl = log_solar_mg_lbl * (mg_25/100)
    log_solar_i_26_lbl = log_solar_mg_lbl * (mg_26/100)
    #convert to Mg/H
    mg_24_H_is_lbl = log_solar_i_24_lbl - (7.53*0.7899)
    mg_25_H_is_lbl = log_solar_i_25_lbl - (7.53*0.1000)
    mg_26_H_is_lbl = log_solar_i_26_lbl - (7.53*0.1101)
    #covert to feh
    log_solar_mg_lbl = (MgH['[X/H]'].values[0]-feh) + 7.53
    log_solar_i_24_lbl = log_solar_mg_lbl * (mg_24/100)
    log_solar_i_25_lbl = log_solar_mg_lbl * (mg_25/100)
    log_solar_i_26_lbl = log_solar_mg_lbl * (mg_26/100)
    #convert to Mg/H
    mg_24_fe_is_lbl = log_solar_i_24_lbl - (7.53*0.7899)
    mg_25_fe_is_lbl = log_solar_i_25_lbl - (7.53*0.1000)
    mg_26_fe_is_lbl = log_solar_i_26_lbl - (7.53*0.1101)
    #do the same for the errors
    #calculate the log abundances
    log_solar_mg = MgH['e[X/H]'].values[0]+7.53
    log_solar_i_24 = log_solar_mg * (mg_24/100)
    log_solar_i_25 = log_solar_mg * (mg_25/100)
    log_solar_i_26 = log_solar_mg * (mg_26/100)
    #convert to Mg/H
    mg_24_H_is_lbl_err = log_solar_i_24 - (7.53*0.7899)
    mg_25_H_is_lbl_err = log_solar_i_25 - (7.53*0.1000)
    mg_26_H_is_lbl_err = log_solar_i_26 - (7.53*0.1101)

    
    # Create a new row for the structured DataFrame
    new_row = {
        's': round(abundances[0], 4), 'd_s': round(errors[0], 4),
        # 'mg': round(abundances[1], 4), 'd_mg': round(errors[1], 4),
        # 'mg_fe': round(abundances[1]-feh, 4), 'd_mg_fe': round(errors[1]-feh, 4),
        'mg_fe24': round(mg_24_H_is-feh, 4), 'd_mg_fe24': round( mg_24_H_is_err, 4),
        'mg_fe25': round(mg_25_H_is-feh, 4), 'd_mg_fe25': round( mg_25_H_is_err, 4),
        'mg_fe26': round(mg_26_H_is-feh, 4), 'd_mg_fe26': round( mg_26_H_is_err, 4),
        # 'i_24': round(abundances[2], 4), 'd_i_24': round(errors[2], 4),
        # 'i_25': round(abundances[3], 4), 'd_i_25': round(errors[3], 4),
        # 'i_26': round(abundances[4], 4), 'd_i_26': round(errors[4], 4),
        # 'R_24': mg_24, 'd_R_24': round(errors[5], 4),
        # 'R_25': mg_25, 'd_R_25': round(errors[6], 4),
        # 'R_26': mg_26, 'd_R_26': round(errors[7], 4),
        'i_24': round(abundances[1], 4), 'd_i_24': round(errors[1], 4),
        'i_25': round(abundances[2], 4), 'd_i_25': round(errors[2], 4),
        'i_26': round(abundances[3], 4), 'd_i_26': round(errors[3], 4),
        'R_24': mg_24, 'd_R_24': round(errors[4], 4),
        'R_25': mg_25, 'd_R_25': round(errors[5], 4),
        'R_26': mg_26, 'd_R_26': round(errors[6], 4),
        'mg24': round(mg_24_H_is,4), 'd_mg24': round(mg_24_H_is_err,4),
        'mg25': round(mg_25_H_is,4), 'd_mg25': round(mg_25_H_is_err,4),
        'mg26': round(mg_26_H_is,4), 'd_mg26': round(mg_26_H_is_err,4),
        'MgH': MgH['[X/H]'].values[0], 'd_MgH': MgH['e[X/H]'].values[0], #The lbl tests with H
        'MgH24': round(mg_24_H_is_lbl,4), 'd_MgH24': round(mg_24_H_is_lbl_err,4),
        'MgH25': round(mg_25_H_is_lbl,4), 'd_MgH25': round(mg_25_H_is_lbl_err,4),
        'MgH26': round(mg_26_H_is_lbl,4), 'd_MgH26': round(mg_26_H_is_lbl_err,4),
        'MgFe': MgH['[X/H]'].values[0]-feh, 'd_MgFe': MgH['e[X/H]'].values[0], #The lbl tests with Fe
        'MgFe24': round(mg_24_fe_is_lbl,4), 'd_MgFe24': round(mg_24_H_is_lbl_err,4),
        'MgFe25': round(mg_25_fe_is_lbl,4), 'd_MgFe25': round(mg_25_H_is_lbl_err,4),
        'MgFe26': round(mg_26_fe_is_lbl,4), 'd_MgFe26': round(mg_26_H_is_lbl_err,4),
    }
    
    # Append the new row to the structured DataFrame with star_name as the index
    final_df = final_df.append(pd.Series(new_row, name=star_name))

# print(final_df)
#save the final_df to a csv file
final_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund_paper_vpass_{vpass}_test.csv')

#%% Plot the things here to not haveto look later
# ----------------------------------------------------------------------------------------------

#%% """Define the region plot lines"""

def region_plots(region, raw,ax):
    if region == 1:    
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
            
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.35, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 1', fontsize=12)
        

        #mg24
        ax.axvline(x=5134.570, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5134.656, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5134.734, ymin=0, color='black',lw=1,alpha=0.5)
        
    '''Region 2'''
    if region == 2:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.26, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 2', fontsize=12)
        
        #mg24
        ax.axvline(x=5138.710, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5138.768, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        ax.axvline(x=5138.785, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5138.826, ymin=0, color='black',lw=1,alpha=0.5)
        ax.axvline(x=5138.862, ymin=0, color='black',lw=1,alpha=0.5)

    '''Region 3'''
    if region == 3:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.45, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 3', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5140.229, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5140.286, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5140.359, ymin=0, color='black',lw=1,alpha=0.5)

    '''Region 4'''
    if region == 4:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.45, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 4', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5134.208, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5134.295, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5134.376, ymin=0, color='black',lw=1,alpha=0.5)
        
    '''Region 5'''
    if region == 5:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.37, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 5', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5135.111, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5135.160, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5135.240, ymin=0, color='black',lw=1,alpha=0.5)

    '''Region 6'''
    if region == 6:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.41, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 6', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5136.123, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5136.087, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5136.144 , ymin=0, color='black',lw=1,alpha=0.5)

    '''Region 7'''
    if region == 7:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.5, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 7', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5136.439, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5136.502, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5136.560, ymin=0, color='black',lw=1,alpha=0.5)

    '''Region 8'''
    if region == 8:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
         #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.33, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 8', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5138.486, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5138.427, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5138.501, ymin=0, color='black',lw=1,alpha=0.5)


    '''Region 9'''
    if region == 9:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
        #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.5, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 9', fontsize=12)
        
        
        #mg24
        ax.axvline(x=5141.234, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5141.288, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5141.338, ymin=0, color='black',lw=1,alpha=0.5)

    '''Region 10'''
    if region == 10:
        try:
            # Force plain numbers on the x-axis
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
            ax.ticklabel_format(style='plain', axis='x')
        except:
            None
        #Find min wavelength
        lw, uw = get_region(region)
        ax.set_xlim(lw - 0.4, uw + 0.5)
        # ax.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        ax.set_ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        ax.fill_between([lw, uw], 0.5, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        ax.set_title('Region 10', fontsize=12)
        
        #mg24
        ax.axvline(x=5133.174, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        ax.axvline(x=5133.231, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        ax.axvline(x=5133.292, ymin=0, color='black',lw=1,alpha=0.5)

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
#%%   """Plot all of the regions best fits in thier own plot"""

import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import glob
import os
from astropy.io import fits
import pandas as pd
import time
import scipy as sp
import sys
import matplotlib
import matplotlib.ticker as mticker
#--- iSpec directory -------------------------------------------------------------
if os.path.exists('/home/users/qai11/iSpec_v20201001'):
    "Location of the files on Uni computer"
    ispec_dir = '/home/users/qai11/iSpec_v20201001'
else:
    "location of data on Mac"
    ispec_dir = '/Users/quin/Desktop/2024_Data/iSpec_v20230804'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec



# Create empty DataFrames to store parameters and errors separately
def isotope_regions(star_name,regions):
    all_params = pd.DataFrame()
    #Read in star_colour information
    star_colour = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    star_colour = star_colour[star_colour['ID2'] == star_name]
    star_colour = star_colour['colour'].values[0]
    #set up plots
    y_sub = 2
    x_sub = int(np.ceil(len(regions) / y_sub))  # Ensure we have enough rows
    #Read in raw spectrum
    # Spectrum from uni computer
    raw = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_5100-5200.txt')
    if len(regions) == 4:
        fig, ax = plt.subplots(x_sub, 2, figsize=(12, 10))
        ax=ax.flatten()  
    elif len(regions) == 3:
        fig, ax = plt.subplots(x_sub, 2, figsize=(12, 10))
        ax=ax.flatten()
    elif len(regions) == 2:
        fig, ax = plt.subplots(x_sub, 2, figsize=(12, 6))
        ax=ax.flatten()
    else:
        fig, ax = plt.subplots(x_sub, 2, figsize=(12, 15))
        ax=ax.flatten()
    #Set the iteration to 0

    iteration = 0
    for region in regions:
        i = iteration
        #Load the best fit values for the region
        try:
            fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_pass_{vpass}_fine.csv', sep=',')
        except:
            fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}.csv', sep=',')
        
        #Create a dataframe with the name of the best fit file
        best_fit = fit_pass.loc[fit_pass['chi_squared'].idxmin()]['filename']
        print(f'Best fit for region {region} is {best_fit}')
        #Open the best fit file
        model_spectra = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/{best_fit}', sep="     ", header=None, skiprows = [0,1])
        # Plot each region in subsequent subplots
        #Call region_plots
        region_plots(region, raw, ax[i])
        
        # plot the synthetic spectrum
        ax[i].plot(model_spectra[0], model_spectra[1], label='Synthetic Spectrum')
        # plot the observed spectrum
        ax[i].plot(raw['waveobs'], raw['flux'] , label='Observed Spectrum', c=star_colour)
        ax[i].set_xlabel('Wavelength ($\AA$)',fontsize=12)
        ax[i].set_ylabel('Flux',fontsize=12)
        ax[i].legend(loc='upper right')
        iteration +=1
        
    #still plot if less than an evan numer of regions
    fig.tight_layout()
    diff_axes = len(ax) - len(regions)
    for i in range(diff_axes):
        ax[-int(i+1)].set_axis_off()

        
    #Save the plot
    # plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    # plt.close()
    #Save for the papers
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Paper_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

vpass = 12
    
# regions = [3,4,5,6,7,8]
# All stars
# import ast
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','moon','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588'] #removed the ones with only 1
# star_list = ['hd_157244']
# star_list = ['hd_18884'] #there is a problem here
# star_list = ['hd_157244'] #Same here
# star_list = ['moon']
star_list = ['hd_10700']
for star in star_list:
    #open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #get the star regions
    # regions = star_info[star_info['ID2'] == star]['regions'].apply(ast.literal_eval).values[0]
    regions = [1]
    # print(regions)
    isotope_regions(star,regions)
    print(f'{star} Done')
    
# %% """replot bc regions with only one dont work"""


def region_1_5(star_name,regions):
    if region == 1:    
        #Find min wavelength
        lw, uw = get_region(region)
        plt.xlim(lw - 0.4, uw + 0.5)
        # plt.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        plt.ylim(min_flux-0.05,1.01)
            
        #Plot the box where the fitting region is
        plt.fill_between([lw, uw], 0.35, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        plt.title('Region 1', fontsize=12)
        

        #mg24
        plt.axvline(x=5134.570, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        plt.axvline(x=5134.656, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        plt.axvline(x=5134.734, ymin=0, color='black',lw=1,alpha=0.5)
        
        '''Region 5'''
    if region == 5:
        #Find min wavelength
        lw, uw = get_region(region)
        plt.xlim(lw - 0.4, uw + 0.5)
        # plt.set_ylim(0.3,1.01)
        cropped_flux = raw[(raw['waveobs'] > lw) & (raw['waveobs'] < uw)]['flux']
        min_flux = cropped_flux.min()
        plt.ylim(min_flux-0.05,1.01)
        
        #Plot the box where the fitting region is
        plt.fill_between([lw, uw], 0.37, 1, facecolor = '#CCDBFD', alpha = 0.3)
        #set a plot title
        plt.title('Region 5', fontsize=12)
        
        
        #mg24
        plt.axvline(x=5135.111, ymin=0, linestyle='-.', color='black',lw=1,alpha=0.5)
        #mg25
        plt.axvline(x=5135.160, ymin=0, linestyle='--', color='black',lw=1,alpha=0.5)
        #mg26
        plt.axvline(x=5135.240, ymin=0, color='black',lw=1,alpha=0.5)

# star_list = ['hd_128620','hd_156098','hd_160691']
# star_list = ['hd_156098']
# star_list = ['hd_160691']

for star_name in star_list:
    #open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #get colours
    star_colour = star_info[star_info['ID2'] == star_name]
    star_colour = star_colour['colour'].values[0]
    #get the star regions
    regions = star_info[star_info['ID2'] == star]['regions'].apply(ast.literal_eval).values[0]
    region = regions[0]
    #Load the best fit values for the region
    try:
        # fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/all_fits_region_{region}_pass_{vpass}.csv', sep=',')
        fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_pass_{vpass}_fine.csv', sep=',')
    except:
        # fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/all_fits_region_{region}.csv', sep=',')
        # fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}.csv', sep=',')
        None
 
    #Create a dataframe with the name of the best fit file
    best_fit = fit_pass.loc[fit_pass['chi_squared'].idxmin()]['filename']
    # print(best_fit)
    
    
    #Open the best fit file
    # model_spectra = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/{best_fit}', sep="     ", header=None, skiprows = [0,1])
    model_spectra = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/{best_fit}', sep="     ", header=None, skiprows = [0,1])
    
    # Spectrum from uni computer
    raw = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_5100-5200.txt')
        
    plt.figure(figsize=(8,4))
    # plot the synthetic spectrum
    plt.plot(model_spectra[0], model_spectra[1], label='Synthetic Spectrum')
    # plot the observed spectrum
    plt.plot(raw['waveobs'], raw['flux'] , label='Observed Spectrum', c=star_colour)
    plt.xlabel('Wavelength ($\AA$)',fontsize=12)
    plt.ylabel('Flux',fontsize=12)
    region_1_5(star,region)
    plt.legend(loc='upper right')
    
    # plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Paper_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    plt.close()
    
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
# plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Results/Isotope_Percentage_vs_wavelength.png', dpi=300, bbox_inches='tight')
plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Paper_Figures/Results/Isotope_Percentage_vs_wavelength_{vpass}.png', dpi=300, bbox_inches='tight')
