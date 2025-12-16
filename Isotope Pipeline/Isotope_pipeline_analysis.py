"""
Title: Isotope_pipeline_analysis.py
Author: Quin Aicken Davies
Date: 12/06/2025

Description: Run after Isotope_pipeline_uncertainties to calculate the final abundances into csv forms.
"""
# %% """Make a table for the isotopic ratios"""
# """Makes the All_isotope_ratios_pre_avg file"""

vpass = 24

import pandas as pd
import numpy as np
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# star_list = ['moon','hd_18907']
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407']
# star_list = ['hd_10700']

# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407'] 

# star_list = ['hd_18884','hd_157244']
# star_list = ['hd_157244']
# star_list = ['hd_18884']
#all paper stars
star_list = ['hd_11695','hd_18884','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407'] 
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
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407']
# star_list = ['hd_11695']
# star_list = ['hd_10700']
# star_list = ['hd_23249','hd_128621','hd_10700']
# star_list = ['hd_18884']
# Initialize an empty dictionary to hold abundance data for each star
abundance_dict = {}
# star_list = ['hd_10700']
star_list = ['hd_18884'] #use vpass 37
# star_list = ['hd_11695','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407'] #use vpass 24
# vpass = '6'
vpass = 24
# star_list = ['hd_18884','hd_157244']
# star_list = ['hd_157244']
# star_list = ['hd_18884']
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
    star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
    #get the star regions
    regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0]
    
    #only take the average of the given regions
    iso_abund_params = iso_abund_params[iso_abund_params['region'].isin(regions)]
    # Separate the error and parameter columns
    errors = iso_abund_params.filter(like='d_')
    params = iso_abund_params.drop(columns=errors.columns).drop(columns=['region', 'pass'])
    # print(params)
    # print(errors)
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
star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
#remove the 10th row
star_info = star_info.drop(10)
#reset index
star_info = star_info.reset_index(drop=True)

# Iterate over each star and add the results to the structured DataFrame
for star_name in star_list:
    feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
    # Extract the abundance and error data for the current star
    abundances = abundance_df.loc[star_name, 'abundance']
    # print(abundances)
    errors = abundance_df.loc[star_name, 'error']
    #Open summary abundances file
    # summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    # summary_abundances_Fe = summary_abundances.sort_values(by=['[X/Fe]', 'e[X/Fe]'], ascending=False)
    # summary_abundances_Fe = summary_abundances.drop_duplicates(subset=['element'], keep='first')
    # summary_abundances_H = summary_abundances.sort_values(by=['[X/H]', 'e[X/H]'], ascending=False)
    # summary_abundances_H = summary_abundances.drop_duplicates(subset=['element'], keep='first')
    # #Extract the Mg [X/H] and error
    # MgFe = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/Fe]','e[X/Fe]']]
    # MgH = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    # print(MgH)
    # MgH['[X/H]'] = MgH['[X/H]'] - feh
    # print(MgH)
    # print(MgH)
    mg_24,mg_25,mg_26 = calc_ratio(abundances[2], abundances[3], abundances[4])
    # print(mg_24)
    # print(abundances[4])
    # print(errors[4], errors[5], errors[6])
    
    #Open summary abundances file for Mg abundance
    summary_abundances = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep='\s+', engine='python')
    summary_abundances = summary_abundances.sort_values(by=['[X/H]', 'e[X/H]'], ascending=False)
    #Extract the Mg [X/H] and error
    Mg_lbl_abundance = summary_abundances.loc[summary_abundances['element']=='Mg',['[X/H]','e[X/H]']]
    Mg_lbl_abundance = Mg_lbl_abundance['[X/H]'].values[0]
    print(f"Mg_lbl_abundance={Mg_lbl_abundance}")
    
    #New calculation for paper including check
    import numpy as np

    def calculate_mg_isotope_abundances(mg_h, solar_log_eps_mg, isotope_fractions_star, isotope_fractions_solar,isotope_fraction_errors_star):
        """
        Calculate [X/H], relative abundance to solar, and ε for Mg isotopes.
        Also performs a consistency check that log ε(Mg) equals log10(sum of isotopic ε values).

        Parameters:
        - mg_h: [Mg/H] value for the star
        - solar_log_eps_mg: solar log abundance of Mg (e.g., 7.53)
        - isotope_fractions_star: dict with isotope fractions in the star (keys: '24Mg', '25Mg', '26Mg')
        - isotope_fractions_solar: dict with isotope fractions in the Sun (keys: '24Mg', '25Mg', '26Mg')
        - isotope_fraction_errors_star: dict with fractional uncertainties in the star's isotope fraction
        Returns:
        - results: dict with isotope data and consistency check
        """
        results = {}
        log_eps_mg_star = mg_h + solar_log_eps_mg
        total_epsilon_star = 0
        total_epsilon_error_squared = 0
        
        for isotope in ['24Mg', '25Mg', '26Mg']:
            f_star = isotope_fractions_star[isotope]
            f_solar = isotope_fractions_solar[isotope]
            f_star_err_frac = isotope_fraction_errors_star[isotope]

            log_eps_isotope_star = log_eps_mg_star + np.log10(f_star)
            log_eps_isotope_solar = solar_log_eps_mg + np.log10(f_solar)

            isotope_xh = log_eps_isotope_star - log_eps_isotope_solar
            relative_abundance = 10 ** isotope_xh
            epsilon_isotope = 10 ** log_eps_isotope_star
            
            # Error propagation
            d_log_eps_isotope_star = f_star_err_frac / (f_star * np.log(10))
            d_epsilon_isotope = epsilon_isotope * f_star_err_frac
            total_epsilon_star += epsilon_isotope
            total_epsilon_error_squared += d_epsilon_isotope ** 2

            results[isotope] = {
                '[X/H]': isotope_xh,
                '[X/H] error': d_log_eps_isotope_star,
                'Relative to Solar': relative_abundance,
                'Relative error': relative_abundance * d_log_eps_isotope_star,
                'Epsilon': epsilon_isotope,
                'Epsilon error': d_epsilon_isotope
            }
        
        # Consistency check
        log_eps_mg_from_isotopes = np.log10(total_epsilon_star)
        d_log_eps_mg_from_isotopes = (1 / (total_epsilon_star * np.log(10))) * np.sqrt(total_epsilon_error_squared)
        log_eps_mg_expected = log_eps_mg_star
        consistent = np.isclose(log_eps_mg_from_isotopes, log_eps_mg_expected, atol=1e-4)

        results['Consistency Check'] = {
            'log ε(Mg) from isotopes': log_eps_mg_from_isotopes,
            'log ε(Mg) error': d_log_eps_mg_from_isotopes,
            'Expected log ε(Mg)': log_eps_mg_expected,
            'Consistent': consistent
        }

        return results

    
    mg_h = Mg_lbl_abundance
    solar_log_eps_mg = 7.53
    isotope_fractions_star = {'24Mg': mg_24/100, '25Mg': mg_25/100, '26Mg': mg_26/100}
    isotope_fractions_solar = {'24Mg': 0.7899, '25Mg': 0.1000, '26Mg': 0.1101}
    isotope_fractions_errors_star = {'24Mg': errors[2]/100, '25Mg': errors[3]/100, '26Mg': errors[4]/100}
    
    results = calculate_mg_isotope_abundances(mg_h, solar_log_eps_mg, isotope_fractions_star,
                                              isotope_fractions_solar,isotope_fractions_errors_star)

    for iso in ['24Mg', '25Mg', '26Mg']:
        r = results[iso]
        print(f"{iso}: [X/H] = {r['[X/H]']:.4f}, Relative = {r['Relative to Solar']:.4f}, ε = {r['Epsilon']:.4e}")

    check = results['Consistency Check']
    print(f"\nConsistency Check:")
    print(f"  log ε(Mg) from isotopes = {check['log ε(Mg) from isotopes']:.4f}")
    print(f"  Expected log ε(Mg) = {check['Expected log ε(Mg)']:.4f}")
    print(f"  Consistent: {check['Consistent']}")
    
    mg_24_H_is=results['24Mg']['[X/H]']
    mg_24_H_is_err=results['24Mg']['[X/H] error']
    mg_25_H_is=results['25Mg']['[X/H]']
    mg_25_H_is_err=results['25Mg']['[X/H] error']
    mg_26_H_is=results['26Mg']['[X/H]']
    mg_26_H_is_err=results['26Mg']['[X/H] error']

    """Old calculation method"""
    # #calculate the log abundances
    # log_solar_mg = Mg_lbl_abundance+7.53
    # print(f"log_solar_mg={log_solar_mg}")
    # log_solar_i_24 = log_solar_mg * np.log10(mg_24/100)
    # print(f"mg24={mg_24/100}")
    # print(f"log_solar_i24={log_solar_i_24}")
    # print(f"{log_solar_i_24/log_solar_mg}") #correct to here
    # log_solar_i_25 = log_solar_mg * np.log10(mg_25/100)
    # log_solar_i_26 = log_solar_mg * np.log10(mg_26/100)
    # #convert to Mg/H
    # # mg_24_H_is = log_solar_i_24 - (7.53*0.7899)
    # mg_24_H_is = log_solar_i_24 - (7.53 + np.log10(0.7899))
    # print(f"mg24H={mg_24_H_is}")
    # mg_25_H_is = log_solar_i_25 - (7.53*0.1000)
    # mg_26_H_is = log_solar_i_26 - (7.53*0.1101) 
    
    # print(f"the ratio of the global and isotope={mg_24_H_is/Mg_lbl_abundance}")
    # #Do the same for the errors
    # #calculate the log abundances
    # log_solar_mg = errors[1]+7.53
    # log_solar_i_24 = log_solar_mg * (mg_24/100)
    # log_solar_i_25 = log_solar_mg * (mg_25/100)
    # log_solar_i_26 = log_solar_mg * (mg_26/100)
    # #convert to Mg/H
    # mg_24_H_is_err = log_solar_i_24 - (7.53*0.7899)
    # mg_25_H_is_err = log_solar_i_25 - (7.53*0.1000)
    # mg_26_H_is_err = log_solar_i_26 - (7.53*0.1101)
    
    #line by line
    #Convert the line by line anundances
    # log_solar_mg_lbl = MgH['[X/H]'].values[0] + 7.53
    # log_solar_i_24_lbl = log_solar_mg_lbl * (mg_24/100)
    # log_solar_i_25_lbl = log_solar_mg_lbl * (mg_25/100)
    # log_solar_i_26_lbl = log_solar_mg_lbl * (mg_26/100)
    # #convert to Mg/H
    # mg_24_H_is_lbl = log_solar_i_24_lbl - (7.53*0.7899)
    # mg_25_H_is_lbl = log_solar_i_25_lbl - (7.53*0.1000)
    # mg_26_H_is_lbl = log_solar_i_26_lbl - (7.53*0.1101)
    # #covert to feh
    # log_solar_mg_lbl = (MgH['[X/H]'].values[0]-feh) + 7.53
    # log_solar_i_24_lbl = log_solar_mg_lbl * (mg_24/100)
    # log_solar_i_25_lbl = log_solar_mg_lbl * (mg_25/100)
    # log_solar_i_26_lbl = log_solar_mg_lbl * (mg_26/100)
    # #convert to Mg/H
    # mg_24_fe_is_lbl = log_solar_i_24_lbl - (7.53*0.7899)
    # mg_25_fe_is_lbl = log_solar_i_25_lbl - (7.53*0.1000)
    # mg_26_fe_is_lbl = log_solar_i_26_lbl - (7.53*0.1101)
    # #do the same for the errors
    # #calculate the log abundances
    # log_solar_mg = MgH['e[X/H]'].values[0]+7.53
    # log_solar_i_24 = log_solar_mg * (mg_24/100)
    # log_solar_i_25 = log_solar_mg * (mg_25/100)
    # log_solar_i_26 = log_solar_mg * (mg_26/100)
    # #convert to Mg/H
    # mg_24_H_is_lbl_err = log_solar_i_24 - (7.53*0.7899)
    # mg_25_H_is_lbl_err = log_solar_i_25 - (7.53*0.1000)
    # mg_26_H_is_lbl_err = log_solar_i_26 - (7.53*0.1101)
    
    """Calculate the Mg/Fe abundances"""
    mg_24_fe_is = mg_24_H_is - feh
    mg_25_fe_is = mg_25_H_is - feh
    mg_26_fe_is = mg_26_H_is - feh
    mg_24_fe_is_err = mg_24_H_is_err
    mg_25_fe_is_err = mg_25_H_is_err
    mg_26_fe_is_err = mg_26_H_is_err

    
    # Create a new row for the structured DataFrame
    new_row = {
        's': round(abundances[0], 4), 'd_s': round(errors[0], 4),
        'mg': round(abundances[1], 4), 'd_mg': round(errors[1], 4),
        # 'mg_fe': round(abundances[1]-feh, 4), 'd_mg_fe': round(errors[1]-feh, 4),
        # 'mg_fe24': round(mg_24_H_is-feh, 4), 'd_mg_fe24': round( mg_24_H_is_err, 4),
        # 'mg_fe25': round(mg_25_H_is-feh, 4), 'd_mg_fe25': round( mg_25_H_is_err, 4),
        # 'mg_fe26': round(mg_26_H_is-feh, 4), 'd_mg_fe26': round( mg_26_H_is_err, 4),
        'i_24': round(abundances[2], 4), 'd_i_24': round(errors[2], 4),
        'i_25': round(abundances[3], 4), 'd_i_25': round(errors[3], 4),
        'i_26': round(abundances[4], 4), 'd_i_26': round(errors[4], 4),
        'R_24': mg_24, 'd_R_24': round(errors[5], 4),
        'R_25': mg_25, 'd_R_25': round(errors[6], 4),
        'R_26': mg_26, 'd_R_26': round(errors[7], 4),
        # 'i_24': round(abundances[1], 4), 'd_i_24': round(errors[1], 4),
        # 'i_25': round(abundances[2], 4), 'd_i_25': round(errors[2], 4),
        # 'i_26': round(abundances[3], 4), 'd_i_26': round(errors[3], 4),
        # 'R_24': mg_24, 'd_R_24': round(errors[4], 4),
        # 'R_25': mg_25, 'd_R_25': round(errors[5], 4),
        # 'R_26': mg_26, 'd_R_26': round(errors[6], 4),
        'mg24': round(mg_24_H_is,4), 'd_mg24': round(mg_24_H_is_err,4),
        'mg25': round(mg_25_H_is,4), 'd_mg25': round(mg_25_H_is_err,4),
        'mg26': round(mg_26_H_is,4), 'd_mg26': round(mg_26_H_is_err,4),
        # 'MgH': MgH['[X/H]'].values[0], 'd_MgH': MgH['e[X/H]'].values[0], #The lbl tests with H
        # 'MgH24': round(mg_24_H_is_lbl,4), 'd_MgH24': round(mg_24_H_is_lbl_err,4),
        # 'MgH25': round(mg_25_H_is_lbl,4), 'd_MgH25': round(mg_25_H_is_lbl_err,4),
        # 'MgH26': round(mg_26_H_is_lbl,4), 'd_MgH26': round(mg_26_H_is_lbl_err,4),
        # 'MgFe': MgH['[X/H]'].values[0]-feh, 'd_MgFe': MgH['e[X/H]'].values[0], #The lbl tests with Fe
        # 'MgFe24': round(mg_24_fe_is_lbl,4), 'd_MgFe24': round(mg_24_H_is_lbl_err,4),
        # 'MgFe25': round(mg_25_fe_is_lbl,4), 'd_MgFe25': round(mg_25_H_is_lbl_err,4),
        # 'MgFe26': round(mg_26_fe_is_lbl,4), 'd_MgFe26': round(mg_26_H_is_lbl_err,4),
    }
    
    # Append the new row to the structured DataFrame with star_name as the index
    final_df = final_df.append(pd.Series(new_row, name=star_name))

# print(final_df)
#save the final_df to a csv file
final_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund_paper_New_hd_18884.csv')

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
import ast
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
    star_colour = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
    star_colour = star_colour[star_colour['ID2'] == star_name]
    star_colour = star_colour['colour'].values[0]
    #set up plots
    y_sub = 2
    x_sub = int(np.ceil(len(regions) / y_sub))  # Ensure we have enough rows
    #Read in raw spectrum
    # Spectrum from uni computer
    """Spectrum for non giants"""
    # raw = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/{star_name}_5100-5200.txt')
    """Spectrum for giants hd18884 and hd157244"""
    raw = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/{star_name}_5100-5200.txt')
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
    #Open the best fit files for each region
    # w_avg_file = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/w_avg_models_vpass_{vpass}.csv', sep=',')
    # #open the file with the model in it
    # w_avg_star_name = w_avg_file[w_avg_file['star_name'] == star_name]['filename'].values[0]
    # w_avg_model = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/{w_avg_star_name}', sep="     ", header=None, skiprows = [0,1])
    # print(f'Weighted average model for {star_name} is {w_avg_star_name}')
    # print('w_avg_model shape:', w_avg_model.shape)
    iteration = 0
    #open the file with the weighted average information
    # weighted_avg_file = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', sep=',')
    #scrape the ratios out of it
    for region in regions:
        i = iteration
        #Load the best fit values for the region
        try:
            fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_pass_{vpass}_fine.csv', sep=',')
            # fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_pass_{vpass}.csv', sep=',')
        except:
            fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}_fine.csv', sep=',')
            # fit_pass = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/all_fits_region_{region}.csv', sep=',')
        #Create a dataframe with the name of the best fit file
        best_fit = fit_pass.loc[fit_pass['chi_squared'].idxmin()]['filename']
        # print(f'Best fit for region {region} is {best_fit}')
        """Open the best fit file for all but giants"""
        model_spectra = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/{best_fit}', sep="     ", header=None, skiprows = [0,1])
        """Open for giants"""
        # model_spectra = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/{best_fit}', sep="     ", header=None, skiprows = [0,1])
        # Plot each region in subsequent subplots
        #Call region_plots
        region_plots(region, raw, ax[i])
        #paste the best fit values for the isotopes onto the plot in the bottom left corner
        ratio = fit_pass.loc[fit_pass['chi_squared'].idxmin()]['ratio']
        ax[i].text(0.05, 0.90, f'Best Fit {ratio}', transform=ax[i].transAxes, fontsize=12, verticalalignment='bottom')
        #paste the weighted average values for the isotopes onto the plot in the bottom left corner
        # ax[i].text(0.05, 0.80, f"Weighted Avg {weighted_avg_file['R_24'].values[0]:.2f}_{weighted_avg_file['R_25'].values[0]:.2f}_{weighted_avg_file['R_26'].values[0]:.2f}", transform=ax[i].transAxes, fontsize=12, verticalalignment='bottom')
        # plot the synthetic spectrum
        ax[i].plot(model_spectra[0], model_spectra[1], label='Synth Spectrum')
        # plot the observed spectrum
        ax[i].plot(raw['waveobs'], raw['flux'] , label='Obs Spectrum', c=star_colour)
        #plot the weighted average specturm overtop
        # ax[i].plot(w_avg_model[0], w_avg_model[1], label='W_Avg Spectrum', c='black', linestyle='--')
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
    # plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    # plt.close()
    #Save for the papers
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    plt.close()

# vpass = 12
    
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
# star_list = ['hd_10700']

#stars less than 5300K
# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    # 'hd_10700','hd_100407']
vpass = 24
# star_list = ['hd_18884','hd_157244']
# star_list = ['hd_157244']
# star_list = ['hd_18884']
star_list = ['hd_18907']

# star_list = ['hd_18884'] #there is a problem here
for star in star_list:
    #open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
    #get the star regions
    regions = star_info[star_info['ID2'] == star]['regions'].apply(ast.literal_eval).values[0]
    # regions = [1,4,5]
    #test regions aug 14th 2025
    # regions = star_info[star_info['ID2'] == star]['old regions'].apply(ast.literal_eval).values[0]
    # regions = [1]
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
    star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
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
    
    # plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/all_fits/Isotope_fits_{star_name}.png', dpi=300, bbox_inches='tight')
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
star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')

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
# plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_Figures/Results/Isotope_Percentage_vs_wavelength.png', dpi=300, bbox_inches='tight')
plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/Isotope_Percentage_vs_wavelength_{vpass}.png', dpi=300, bbox_inches='tight')

#%% Plot element vs isotopes abundances
vpass = 24
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
star_list = ['hd_11695','hd_18884','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
# # star_list = ['hd_11695']
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
    # isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', delimiter=',')
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_New.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    # print(f"star {iso_mg['Unnamed: 0']} with Mg/H {iso_mg['mg']}")
    # print(f"star {Mg_values['star_name']} with Mg/H {Mg_values['[Mg/H]']}")
    # #Make a variable when the mg values are less than 0.25
    # small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    # #for the related stars in the small_mg variable, find the Mg values
    # small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    # small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    # small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(4, 2, figsize=(9, 12), constrained_layout=True)
    # increase the space between the plots
    # fig.subplots_adjust(hspace=0.35, wspace=0.35)

    #plot Mg vs Eu
    print(f"mg vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], Mg_values['[Mg/H]'])}")
    ax[0, 0].errorbar(Eu_values['[Eu/H]'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg'], yerr=small_mg['d_mg'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[0, 0].set_ylabel('[Mg/H]', fontsize=12)
    ax[0, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[0, 0].set_ylim(-0.3, 0.6)
    ax[0, 0].set_xlim(-0.5, 0.5)
    #plot Mg vs Ba
    print(f"mg vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], Mg_values['[Mg/H]'])}")
    ax[0, 1].errorbar(Ba_values['[Ba/H]'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # print(Ba_values['[Ba/H]'])
    # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg'], yerr=small_mg['d_mg'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[0, 1].set_ylabel('[Mg/H]', fontsize=12)
    ax[0, 1].set_xlabel('[Ba/H]', fontsize=12)
    ax[0, 1].set_ylim(-0.3, 0.6)
    
    
    # Plot mg24 vs Eu
    print(f"mg24 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg24'])}")
    ax[1, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 0].set_ylabel('[$^{24}$Mg/H]', fontsize=12)
    ax[1, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[1, 0].set_ylim(-0.3, 0.6)
    ax[1, 0].set_xlim(-0.5, 0.5)

    # Plot mg25 vs Eu
    print(f"mg25 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg25'])}")
    ax[2, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[1, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[2, 0].set_ylabel('[$^{25}$Mg/H]', fontsize=12)
    # ax[2, 0].set_ylim(-0.5, 0.5)
    ax[2, 0].set_xlim(-0.5, 0.5)

    # Plot mg26 vs Eu
    print(f"mg26 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg26'])}")
    ax[3, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[2, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[3, 0].set_ylabel('[$^{26}$Mg/H]', fontsize=12)
    ax[3, 0].set_ylim(-0.3, 0.7)
    ax[3, 0].set_xlim(-0.5, 0.5)

    # Plot mg24 vs Ba
    print(f"mg24 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg24'])}")
    ax[1, 1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 1].set_ylabel('[$^{24}$Mg/H]', fontsize=12)
    ax[1, 1].set_xlabel('[Ba/H]', fontsize=12)
    ax[1, 1].set_ylim(-0.3, 0.6)

    # Plot mg25 vs Ba
    print(f"mg25 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg25'])}")
    ax[2,1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[1,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,1].set_xlabel('[Ba/H]', fontsize=12)
    ax[2,1].set_ylabel('[$^{25}$Mg/H]', fontsize=12)
    # ax[2, 1].set_ylim(-0.5, 0.5)

    # Plot mg26 vs Ba
    print(f"mg26 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg26'])}")
    ax[3, 1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[2, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3, 1].set_xlabel('[Ba/H]', fontsize=12)
    ax[3, 1].set_ylabel('[$^{26}$Mg/H]', fontsize=12)
    ax[3, 1].set_ylim(-0.3, 0.7)

    # Save the plot
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/Element_fits_X_H_new_calc.png', dpi=300, bbox_inches='tight')

element_plots_XH_new(star_list)
# %% shuffling parameters teff, logg, feh
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
# # star_list = ['hd_11695']
# element = ["Eu", "Ba", "Mg"]
vpass = 24
star_list = ['hd_11695','hd_18884','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
# # star_list = ['hd_11695']
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
    star_info = pd.read_csv(f'/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv', sep=',')
    #remove the 10th row
    star_info = star_info.drop(10)
    #reset index
    star_info = star_info.reset_index(drop=True)
    
    for star_name in star_list:
        #Extract the FEH value from the masters stars csv
        feh = star_info[star_info['ID2'] == star_name]['FEH'].values[0]
        #make a filter to only take the stars in the star_list and override the star_info variable
        star_info = star_info[star_info['ID2'].isin(star_list)]
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
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg','d_mg']]
    
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
    fig, ax = plt.subplots(3, 3, figsize=(16,12),constrained_layout=True)

    #increase the space between the plots
    fig.subplots_adjust(hspace=0.3, wspace=0.3)
    
    #Plot Eu vs Mg
    # print(scipy.stats.pearsonr(Mg_values['[Mg/Fe]'], Eu_values['[Eu/Fe]']))
    print(f"Eu vs Mg: {scipy.stats.pearsonr(Mg_values_H['[Mg/H]'], Eu_values_H['[Eu/H]'])}")
    ax[0,0].errorbar(Mg_values_H['[Mg/H]'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], xerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    # ax[0,0].errorbar(small_mg_values['[Mg/H]'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[0,0].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,0].set_ylim(-0.4,0.9)
    
    #Plot Eu vs Ba
    print(f"Eu vs Ba: {scipy.stats.pearsonr(Ba_values_H['[Ba/H]'], Eu_values_H['[Eu/H]'])}")
    ax[0,1].errorbar(Ba_values_H['[Ba/H]'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], xerr=Ba_values_H['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[0,1].errorbar(small_ba_values['[Ba/H]'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[0,1].set_xlabel('[Ba/H]',fontsize=12)
    ax[0,1].set_ylabel('[Eu/H]',fontsize=12)
    ax[0,1].set_ylim(-0.4,1)  
      
    #Plot Ba vs Mg
    print(f"Mg vs Ba: {scipy.stats.pearsonr(Mg_values_H['[Mg/H]'], Ba_values_H['[Ba/H]'])}")
    ax[0,2].errorbar(Mg_values_H['[Mg/H]'], Ba_values_H['[Ba/H]'], yerr=Ba_values_H['e[Ba/H]'], xerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    # ax[0,2].errorbar(small_mg_values['[Mg/H]'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], xerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[0,2].set_xlabel('[Mg/H]',fontsize=12)
    ax[0,2].set_ylabel('[Ba/H]',fontsize=12)
    #Plot Eu vs Teff

    print(f"Eu vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Eu_values_H['[Eu/H]'])}")
    ax[1,0].errorbar(star_info['TEFF'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[1,0].errorbar(small_star_teff['TEFF'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,0].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[1,0].set_ylim(-0.4,0.9)
    
    #Plot Ba vs Teff
    print(f"Ba vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Ba_values_H['[Ba/H]'])}")
    ax[1,1].errorbar(star_info['TEFF'], Ba_values_H['[Ba/H]'], yerr=Ba_values_H['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[1,1].errorbar(small_star_teff['TEFF'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,1].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,1].set_ylabel('[Ba/H]',fontsize=12)
    
    
    #Plot Eu vs logg
    print(f"Eu vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Eu_values_H['[Eu/H]'])}")
    ax[2,0].errorbar(star_info['LOGG'], Eu_values_H['[Eu/H]'], yerr=Eu_values_H['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[2,0].errorbar(small_star_logg['LOGG'], small_eu_values['[Eu/H]'], yerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[2,0].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,0].set_ylabel('[Eu/H]',fontsize=12)
    ax[2,0].set_ylim(-0.4,0.9)
    
    #Plot Ba vs logg
    print(f"Ba vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Ba_values_H['[Ba/H]'])}")
    ax[2,1].errorbar(star_info['LOGG'], Ba_values_H['[Ba/H]'], yerr=Ba_values_H['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # ax[2,1].errorbar(small_star_logg['LOGG'], small_ba_values['[Ba/H]'], yerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[2,1].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,1].set_ylabel('[Ba/H]',fontsize=12)
    
    #Plot Mg vs Teff
    print(f"Mg vs Teff: {scipy.stats.pearsonr(star_info['TEFF'], Mg_values_H['[Mg/H]'])}")
    ax[1,2].errorbar(star_info['TEFF'], Mg_values_H['[Mg/H]'], yerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    # ax[1,2].errorbar(small_star_teff['TEFF'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[1,2].set_xlabel('$T_{eff}$',fontsize=12)
    ax[1,2].set_ylabel('[Mg/H]',fontsize=12)
    
    #Plot Mg vs logg
    print(f"Mg vs logg: {scipy.stats.pearsonr(star_info['LOGG'], Mg_values_H['[Mg/H]'])}")
    ax[2,2].errorbar(star_info['LOGG'], Mg_values_H['[Mg/H]'], yerr=Mg_values_H['e[Mg/H]'], fmt='o', elinewidth=0.5)
    # ax[2,2].errorbar(small_star_logg['LOGG'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], fmt='o', elinewidth=0.5,color='orange')
    ax[2,2].set_xlabel('$\log g (cm/s^{2}$)',fontsize=12)
    ax[2,2].set_ylabel('[Mg/H]',fontsize=12)
    
    #Save the plot
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/Element_fits_vH_params.png', dpi=300, bbox_inches='tight')



single_el_vs_params(star_list) 

# %% Three synthetic models
"""Plot three synthetic spectra ontop of the observed spectrum One with 
Mg24, one with Mg25 and one with Mg26 and one with all three of a best fit result"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
star_name = 'hd_18907'
#open Masters_stars.csv for finding the colour of the star
stars = pd.read_csv('/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv')
star_colour = stars[stars['ID2'] == star_name]['colour'].values[0]
#Open the files for plotting
smoothed_24 = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/out_s57_mg041_i23_100_100_rv0', sep="     ", header=None, skiprows = [0,1])
smoothed_25 = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/out_s57_mg041_i100_105_100_rv0', sep="     ", header=None, skiprows = [0,1])
smoothed_26 = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/out_s57_mg041_i100_100_13_rv0', sep="     ", header=None, skiprows = [0,1])
smoothed_all = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/out_s57_mg041_i23_105_130_rv0', sep="     ", header=None, skiprows = [0,1])
raw = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests/{star_name}_5100-5200.txt', sep="	", header=None)
#make a subfigure 
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
#Plot the observed spectrum
lw, uw = get_region(region)
cropped_flux = raw[(raw[0] > lw) & (raw[0] < uw)][1]
min_flux = cropped_flux.min()

# # Force plain numbers on the x-axis
# axs.xaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=False))
# axs.ticklabel_format(style='plain', axis='x')
# Plot the 24 isotope
axs[0, 0].plot(smoothed_24[0], smoothed_24[1])
axs[0, 0].plot(raw[0], raw[1], c = star_colour)
axs[0, 0].text(5134.208, 0.72, 'Mg24 \n5134.208', fontsize=12, color='black',horizontalalignment='center')
axs[0, 0].text(5134.570, 0.72, 'Mg24 \n5134.570', fontsize=12, color='black',horizontalalignment='center')
axs[0, 0].text(5135.111, 0.72, 'Mg24 \n5135.111', fontsize=12, color='black',horizontalalignment='center')
axs[0, 0].annotate('', xy=(5134.208, 0.83), xytext=(5134.208, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[0, 0].annotate('', xy=(5134.570, 0.775), xytext=(5134.570, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[0, 0].annotate('', xy=(5135.111, 0.915), xytext=(5135.111, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[0, 0].set_xlim(lw - 0.4, uw + 0.5)
axs[0, 0].set_ylim(min_flux-0.05,1.01)
axs[0, 0].set_ylabel('Flux', fontsize=12)
axs[0, 0].set_xlabel('Wavelength ($\AA$)', fontsize=12)
axs[0, 0].set_title('Mg24')
# axs[0, 0].savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Analysis/{star_name}_Model_Mg26.png', dpi=300, bbox_inches='tight')
# Plot the 25 isotope
axs[0, 1].plot(smoothed_25[0], smoothed_25[1])
axs[0, 1].plot(raw[0], raw[1], c = star_colour)
axs[0, 1].text(5134.295, 0.72, 'Mg25 \n5134.295', fontsize=12, color='black',horizontalalignment='center')
axs[0, 1].text(5134.656, 0.72, 'Mg25 \n5134.656', fontsize=12, color='black',horizontalalignment='center')
axs[0, 1].text(5135.160, 0.72, 'Mg25 \n5135.160', fontsize=12, color='black',horizontalalignment='center')
axs[0, 1].annotate('', xy=(5134.295, 0.95), xytext=(5134.295, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[0, 1].annotate('', xy=(5134.656, 0.935), xytext=(5134.656, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[0, 1].annotate('', xy=(5135.160, 0.955), xytext=(5135.160, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[0, 1].set_xlim(lw - 0.4, uw + 0.5)
axs[0, 1].set_ylim(min_flux-0.05,1.01)
axs[0, 1].set_xlabel('Wavelength ($\AA$)', fontsize=12)
axs[0, 1].set_ylabel('Flux', fontsize=12)
axs[0, 1].set_title('Mg25')
# plt.savefig('/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Analysis/{star_name}_Model_Mg25.png', dpi=300, bbox_inches='tight')
# Plot the 26 isotope
axs[1, 0].plot(smoothed_26[0], smoothed_26[1])
axs[1,0].plot(raw[0], raw[1], c = star_colour)
axs[1, 0].text(5134.376, 0.72, 'Mg26 \n5134.376', fontsize=12, color='black',horizontalalignment='center')
axs[1, 0].text(5134.734, 0.72, 'Mg26 \n5134.734', fontsize=12, color='black',horizontalalignment='center')
axs[1, 0].text(5135.24, 0.72, 'Mg26 \n5135.24', fontsize=12, color='black',horizontalalignment='center')
axs[1, 0].annotate('', xy=(5134.376, 0.96), xytext=(5134.376, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[1, 0].annotate('', xy=(5134.734, 0.955), xytext=(5134.734, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[1,0].annotate('', xy=(5135.24, 0.97), xytext=(5135.24, 0.75), arrowprops=dict(arrowstyle='->', color='black'))
axs[1, 0].set_xlim(lw - 0.4, uw + 0.5)
axs[1, 0].set_ylim(min_flux-0.05,1.01)
axs[1, 0].set_ylabel('Flux', fontsize=12)
axs[1, 0].set_xlabel('Wavelength ($\AA$)', fontsize=12)
axs[1, 0].set_title('Mg26')
# axs[1, 0].savefig('/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Analysis/{star_name}_Model_Mg24.png', dpi=300, bbox_inches='tight')
axs[1, 1].plot(smoothed_all[0], smoothed_all[1])
axs[1, 1].plot(raw[0], raw[1], c = star_colour)
axs[1, 1].set_xlim(lw - 0.4, uw + 0.5)
axs[1, 1].set_ylim(min_flux-0.05,1.01)
axs[1, 1].set_ylabel('Flux', fontsize=12)
axs[1, 1].set_xlabel('Wavelength ($\AA$)', fontsize=12)
axs[1,1].set_title('Mg24, Mg25 and Mg26')
#Draw a square around the fitting region
axs[1, 1].fill_between([lw, uw], min_flux - 0.01, 1, facecolor = '#CCDBFD', alpha = 0.5)
fig.subplots_adjust(hspace=0.3)  # Increase spacing
#Save figure
plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Analysis/{star_name}_Model_Mg24_Mg25_Mg26.png', dpi=300, bbox_inches='tight')


#%% Print the plot with Mg against Mg24
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

vpass = 'New'
star_list = ['hd_11695','hd_18884','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
# # star_list = ['hd_11695']
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
    # isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', delimiter=',')
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_New.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    
    #Make a variable when the mg values are less than 0.25
    small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    #for the related stars in the small_mg variable, find the Mg values
    small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    # #Plot Mg vs mg24
    # print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Mg_values['[Mg/H]'])}")
    # plt.errorbar(iso_mg['mg24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5)
    # # plt.errorbar(small_mg['mg24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    # plt.xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    # plt.ylabel('[Mg/H]',fontsize=12)
    
    # #Plot Mg vs mg25
    # print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['mg25'], Mg_values['[Mg/H]'])}")
    # plt.errorbar(iso_mg['mg25'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg25'], fmt='o', elinewidth=0.5)
    # # plt.errorbar(small_mg['mg24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    # plt.xlabel('IS[$^{25}$Mg/H]',fontsize=12)
    # plt.ylabel('[Mg/H]',fontsize=12)

    # #Plot Mg vs mg26
    # print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['mg26'], Mg_values['[Mg/H]'])}")
    # plt.errorbar(iso_mg['mg26'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg25'], fmt='o', elinewidth=0.5)
    # # plt.errorbar(small_mg['mg24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    # plt.xlabel('IS[$^{26}$Mg/H]',fontsize=12)
    # plt.ylabel('[Mg/H]',fontsize=12)
    
    # Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(3, 1, figsize=(12,20), constrained_layout=True)
    #increase the space between the plots
    fig.subplots_adjust(hspace=0.35, wspace=0.35)
    
    print(f"Mg vs mg24: {scipy.stats.pearsonr(iso_mg['mg24'], Mg_values['[Mg/H]'])}")
    ax[0].errorbar(iso_mg['mg24'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg24'], fmt='o', elinewidth=0.5)
    # ax[1].errorbar(small_mg['mg24'], small_mg_values['[Mg/H]'], yerr=small_mg_values['e[Mg/H]'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    ax[0].set_xlabel('IS[$^{24}$Mg/H]',fontsize=12)
    ax[0].set_ylabel('[Mg/H]',fontsize=12)
    
    #Plot Mg vs mg25
    print(f"Mg vs mg25: {scipy.stats.pearsonr(iso_mg['mg25'], Mg_values['[Mg/H]'])}")
    ax[1].errorbar(iso_mg['mg25'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg25'], fmt='o', elinewidth=0.5)
    # ax[2].errorbar(small_mg['mg24'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    ax[1].set_xlabel('IS[$^{25}$Mg/H]',fontsize=12)
    ax[1].set_ylabel('[Mg/H]',fontsize=12)
    
    #Plot Mg vs mg26
    print(f"Mg vs mg26: {scipy.stats.pearsonr(iso_mg['mg24'], Mg_values['[Mg/H]'])}")
    ax[2].errorbar(iso_mg['mg26'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=iso_mg['d_mg26'], fmt='o', elinewidth=0.5)
    # ax[3].errorbar(small_mg['mg24'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_mg['d_mg24'], fmt='o', elinewidth=0.5, color='orange')
    ax[2].set_xlabel('IS[$^{26}$Mg/H]',fontsize=12)
    ax[2].set_ylabel('[Mg/H]',fontsize=12)
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/Mg_vs_mg_iso_new_calc.png', dpi=300, bbox_inches='tight')

element_plots_XH_new(star_list)
#%% Plot element vs isotopes abundances with Y
vpass = 24
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
star_list = ['hd_11695','hd_18884','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
# # star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg", "Y"]

def element_plots_XH_new(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/H]', 'e[Eu/H]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/H]', 'e[Ba/H]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/H]', 'e[Mg/H]'])
    Y_values = pd.DataFrame(columns=['star_name', '[Y/H]', 'e[Y/H]'])
    
    for star_name in star_list:
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}_with_Y.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/H]', 'e[X/H]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/H]', 'e[X/H]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/H]', 'e[X/H]']].copy()
        y_row = elements[elements['element'] == 'Y'][['[X/H]', 'e[X/H]']].copy()

        # Add the star_name column
        eu_row.insert(0, 'star_name', star_name)
        ba_row.insert(0, 'star_name', star_name)
        mg_row.insert(0, 'star_name', star_name)
        y_row.insert(0, 'star_name', star_name)

        # Rename columns to reflect the element name
        eu_row.columns = ['star_name', '[Eu/H]', 'e[Eu/H]']
        ba_row.columns = ['star_name', '[Ba/H]', 'e[Ba/H]']
        mg_row.columns = ['star_name', '[Mg/H]', 'e[Mg/H]']
        y_row.columns = ['star_name', '[Y/H]', 'e[Y/H]']

        # Append to the main DataFrames
        Eu_values = pd.concat([Eu_values, eu_row], ignore_index=True)
        Ba_values = pd.concat([Ba_values, ba_row], ignore_index=True)
        Mg_values = pd.concat([Mg_values, mg_row], ignore_index=True)
        Y_values = pd.concat([Y_values, y_row], ignore_index=True)
    print(len(Y_values))
    print(len(Mg_values))
        # # Find the stars with the Mg/H values less than 0.25 and add then to a seperate variable
        # if Mg_values['[Mg/H]'] < 0.25:
        #     Mg_small = Mg_values          
            

    #Open the isotope mg abundance file
    # isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', delimiter=',')
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_New.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26']]
    # print(f"star {iso_mg['Unnamed: 0']} with Mg/H {iso_mg['mg']}")
    # print(f"star {Mg_values['star_name']} with Mg/H {Mg_values['[Mg/H]']}")
    # #Make a variable when the mg values are less than 0.25
    # small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    # #for the related stars in the small_mg variable, find the Mg values
    # small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    # small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    # small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(4, 2, figsize=(9, 12), constrained_layout=True)
    # increase the space between the plots
    # fig.subplots_adjust(hspace=0.35, wspace=0.35)

    #plot Mg vs Eu
    print(f"mg vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], Mg_values['[Mg/H]'])}")
    ax[0, 0].errorbar(Eu_values['[Eu/H]'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg'], yerr=small_mg['d_mg'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[0, 0].set_ylabel('[Mg/H]', fontsize=12)
    ax[0, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[0, 0].set_ylim(-0.3, 0.6)
    ax[0, 0].set_xlim(-0.5, 0.5)
    #plot Mg vs Ba
    # print(f"mg vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], Mg_values['[Mg/H]'])}")
    # ax[0, 1].errorbar(Ba_values['[Ba/H]'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # # print(Ba_values['[Ba/H]'])
    # # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg'], yerr=small_mg['d_mg'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[0, 1].set_ylabel('[Mg/H]', fontsize=12)
    # ax[0, 1].set_xlabel('[Ba/H]', fontsize=12)
    # ax[0, 1].set_ylim(-0.3, 0.6)
    
    #plot Mg vs Y
    print(f"mg vs Y: {scipy.stats.pearsonr(Y_values['[Y/H]'], Mg_values['[Mg/H]'])}")
    ax[0, 1].errorbar(Y_values['[Y/H]'], Mg_values['[Mg/H]'], yerr=Mg_values['e[Mg/H]'], xerr=Y_values['e[Y/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg'], yerr=small_mg['d_mg'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[0, 1].set_ylabel('[Mg/H]', fontsize=12)
    ax[0, 1].set_xlabel('[Y/H]', fontsize=12)
    ax[0, 1].set_ylim(-0.3, 0.6)
    ax[0, 1].set_xlim(-1, 0.5)
    
    # Plot mg24 vs Eu
    print(f"mg24 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg24'])}")
    ax[1, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 0].set_ylabel('[$^{24}$Mg/H]', fontsize=12)
    ax[1, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[1, 0].set_ylim(-0.3, 0.6)
    ax[1, 0].set_xlim(-0.5, 0.5)

    # Plot mg25 vs Eu
    print(f"mg25 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg25'])}")
    ax[2, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[1, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[2, 0].set_ylabel('[$^{25}$Mg/H]', fontsize=12)
    # ax[2, 0].set_ylim(-0.5, 0.5)
    ax[2, 0].set_xlim(-0.5, 0.5)

    # Plot mg26 vs Eu
    print(f"mg26 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/H]'], iso_mg['mg26'])}")
    ax[3, 0].errorbar(Eu_values['[Eu/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5)
    # ax[2, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3, 0].set_xlabel('[Eu/H]', fontsize=12)
    ax[3, 0].set_ylabel('[$^{26}$Mg/H]', fontsize=12)
    ax[3, 0].set_ylim(-0.3, 0.7)
    ax[3, 0].set_xlim(-0.5, 0.5)

    # # Plot mg24 vs Ba
    # print(f"mg24 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg24'])}")
    # ax[1, 1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[1, 1].set_ylabel('[$^{24}$Mg/H]', fontsize=12)
    # ax[1, 1].set_xlabel('[Ba/H]', fontsize=12)
    # ax[1, 1].set_ylim(-0.3, 0.6)
    
    # # Plot mg25 vs Ba
    # print(f"mg25 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg25'])}")
    # ax[2,1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # # ax[1,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[2,1].set_xlabel('[Ba/H]', fontsize=12)
    # ax[2,1].set_ylabel('[$^{25}$Mg/H]', fontsize=12)
    # # ax[2, 1].set_ylim(-0.5, 0.5)
    
    # Plot mg26 vs Ba
    # print(f"mg26 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/H]'], iso_mg['mg26'])}")
    # ax[3, 1].errorbar(Ba_values['[Ba/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5)
    # # ax[2, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    # ax[3, 1].set_xlabel('[Ba/H]', fontsize=12)
    # ax[3, 1].set_ylabel('[$^{26}$Mg/H]', fontsize=12)
    # ax[3, 1].set_ylim(-0.3, 0.7)
    
    #Plot mg24 vs Y
    print(f"mg24 vs Y: {scipy.stats.pearsonr(Y_values['[Y/H]'], iso_mg['mg24'])}")
    ax[1, 1].errorbar(Y_values['[Y/H]'], iso_mg['mg24'], yerr=iso_mg['d_mg24'], xerr=Y_values['e[Y/H]'], fmt='o', elinewidth=0.5)
    # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 1].set_ylabel('[$^{24}$Mg/H]', fontsize=12)
    ax[1, 1].set_xlabel('[Y/H]', fontsize=12)
    ax[1, 1].set_ylim(-0.3, 0.6)
    ax[1, 1].set_xlim(-1, 0.5)
    
    #Plot mg25 vs Y
    print(f"mg25 vs Y: {scipy.stats.pearsonr(Y_values['[Y/H]'], iso_mg['mg25'])}")
    ax[2,1].errorbar(Y_values['[Y/H]'], iso_mg['mg25'], yerr=iso_mg['d_mg25'], xerr=Y_values['e[Y/H]'], fmt='o', elinewidth=0.5)
    # ax[1,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,1].set_xlabel('[Y/H]', fontsize=12)
    ax[2,1].set_ylabel('[$^{25}$Mg/H]', fontsize=12)
    # ax[2, 1].set_ylim(-0.5, 0.5)
    ax[2, 1].set_xlim(-0.5, 0.5)
    ax[2, 1].set_xlim(-1, 0.5)
    
    #Plot mg26 vs Y
    print(f"mg26 vs Y: {scipy.stats.pearsonr(Y_values['[Y/H]'], iso_mg['mg26'])}")
    ax[3, 1].errorbar(Y_values['[Y/H]'], iso_mg['mg26'], yerr=iso_mg['d_mg26'], xerr=Y_values['e[Y/H]'], fmt='o', elinewidth=0.5)
    # ax[2, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3, 1].set_xlabel('[Y/H]', fontsize=12)
    ax[3, 1].set_ylabel('[$^{26}$Mg/H]', fontsize=12)
    ax[3, 1].set_ylim(-0.3, 0.7)
    ax[3, 1].set_xlim(-1, 0.5)

    # Save the plot
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/Element_fits_X_H_new_calc_with_Y.png', dpi=300, bbox_inches='tight')

element_plots_XH_new(star_list)

# %% Do the same plot but for Mg/Fe vs isotopes
vpass = 24
import pandas as pd
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import numpy as np
import scipy

# star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
#     'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
#     'hd_102870','hd_45588','hd_156098']
star_list = ['hd_11695','hd_18884','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407']
# # star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]


def element_plots_XFe_new(star_name):
    """Create a plot for Eu, Ba, Mg vs Mg"""
    # Initialize empty DataFrames with star_name as the first column
    Eu_values = pd.DataFrame(columns=['star_name', '[Eu/Fe]', 'e[Eu/Fe]'])
    Ba_values = pd.DataFrame(columns=['star_name', '[Ba/Fe]', 'e[Ba/Fe]'])
    Mg_values = pd.DataFrame(columns=['star_name', '[Mg/Fe]', 'e[Mg/Fe]'])
    star_teff = pd.DataFrame(columns=['star_name', 'TEFF'])
    
    for star_name in star_list:
        teff = pd.read_csv('/home/users/qai11/Documents/Isotope-Pipeline/Masters_stars.csv')
        
        star_teff_value = teff.loc[teff['ID2'] == star_name, 'TEFF'].values[0]
        # Append to the star_teff DataFrame
        star_teff_row = pd.DataFrame({'star_name': [star_name], 'TEFF': [star_teff_value]})
        
        # print(star_teff)
        # Open the lbl abundances
        file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
        elements = pd.read_csv(file_path, delimiter=' ')

        # Extract [X/H] and e[X/H] for each element and store them in a DataFrame
        eu_row = elements[elements['element'] == 'Eu_2'][['[X/Fe]', 'e[X/Fe]']].copy()
        ba_row = elements[elements['element'] == 'Ba'][['[X/Fe]', 'e[X/Fe]']].copy()
        mg_row = elements[elements['element'] == 'Mg'][['[X/Fe]', 'e[X/Fe]']].copy()

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
    # isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_vpass_{vpass}.csv', delimiter=',')
    isotope = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/Isotope_abund_files/weighted_avg_iso_abund_paper_New.csv', delimiter=',')
    #Extract the mg, d_mg, mg24, mg_24_err, mg25, d_mg25, mg26, d_mg26 columns
    iso_mg = isotope[['Unnamed: 0','mg', 'd_mg', 'mg24', 'd_mg24', 'mg25', 'd_mg25', 'mg26', 'd_mg26','MgFe'
                      ,'d_MgFe','MgFe24', 'd_MgFe24','MgFe25', 'd_MgFe25','MgFe26', 'd_MgFe26']]
    # print(f"star {iso_mg['Unnamed: 0']} with Mg/H {iso_mg['mg']}")
    # print(f"star {Mg_values['star_name']} with Mg/H {Mg_values['[Mg/H]']}")
    # #Make a variable when the mg values are less than 0.25
    # small_mg = iso_mg[np.logical_and(-0.25 < iso_mg['mg'], iso_mg['mg'] < 0.25)]
    # #for the related stars in the small_mg variable, find the Mg values
    # small_mg_values = Mg_values[Mg_values['star_name'].isin(small_mg['Unnamed: 0'])]
    # small_eu_values = Eu_values[Eu_values['star_name'].isin(small_mg['Unnamed: 0'])]
    # small_ba_values = Ba_values[Ba_values['star_name'].isin(small_mg['Unnamed: 0'])]
    
    #Make the plot for X/H and mg from IS.
    #Plot the elements vs each Mg and the isotope mg, mg24, mg25, mg26
    fig, ax = plt.subplots(4, 2, figsize=(9, 12), constrained_layout=True)
    # increase the space between the plots
    # fig.subplots_adjust(hspace=0.35, wspace=0.35)

    #plot Mg vs Eu
    print(f"mg vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/Fe]'], iso_mg['MgFe'])}")
    ax[0, 0].errorbar(Eu_values['[Eu/Fe]'], iso_mg['MgFe'], yerr=iso_mg['d_MgFe'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    ax[0, 0].set_ylabel('[Mg/Fe]', fontsize=12)
    ax[0, 0].set_xlabel('[Eu/Fe]', fontsize=12)
    ax[0, 0].set_ylim(-0.3, 0.6)
    ax[0, 0].set_xlim(-0.2, 0.5)
    
    
    #plot Mg vs Ba
    print(f"mg vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], iso_mg['MgFe'])}")
    ax[0, 1].errorbar(Ba_values['[Ba/Fe]'], iso_mg['MgFe'], yerr=iso_mg['d_MgFe'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    ax[0, 1].set_ylabel('[Mg/Fe]', fontsize=12)
    ax[0, 1].set_xlabel('[Ba/Fe]', fontsize=12)
    ax[0, 1].set_ylim(-0.3, 0.6)

    
    
    # Plot mg24 vs Eu
    print(f"mg24 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/Fe]'], iso_mg['MgFe24'])}")
    ax[1, 0].errorbar(Eu_values['[Eu/Fe]'], iso_mg['MgFe24'], yerr=iso_mg['d_MgFe24'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    # ax[0, 0].errorbar(small_eu_values['[Eu/Fe]'], small_mg['MgFe24'], yerr=small_mg['d_mg24'], xerr=small_eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 0].set_ylabel('[$^{24}$Mg/Fe]', fontsize=12)
    ax[1, 0].set_xlabel('[Eu/Fe]', fontsize=12)
    ax[1, 0].set_ylim(-0.2, 0.05)
    ax[1, 0].set_xlim(-0.2, 0.5)


    # Plot mg25 vs Eu
    print(f"mg25 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/Fe]'], iso_mg['MgFe25'])}")
    ax[2, 0].errorbar(Eu_values['[Eu/Fe]'], iso_mg['MgFe25'], yerr=iso_mg['d_MgFe25'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    # ax[1, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2, 0].set_xlabel('[Eu/Fe]', fontsize=12)
    ax[2, 0].set_ylabel('[$^{25}$Mg/Fe]', fontsize=12)
    # ax[2, 0].set_ylim(-0.5, 0.5)
    ax[2, 0].set_xlim(-0.2, 0.5)


    # Plot mg26 vs Eu
    print(f"mg26 vs Eu: {scipy.stats.pearsonr(Eu_values['[Eu/Fe]'], iso_mg['MgFe26'])}")
    ax[3, 0].errorbar(Eu_values['[Eu/Fe]'], iso_mg['MgFe26'], yerr=iso_mg['d_MgFe26'], xerr=Eu_values['e[Eu/Fe]'], fmt='o', elinewidth=0.5)
    # ax[2, 0].errorbar(small_eu_values['[Eu/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_eu_values['e[Eu/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3, 0].set_xlabel('[Eu/Fe]', fontsize=12)
    ax[3, 0].set_ylabel('[$^{26}$Mg/Fe]', fontsize=12)
    ax[3, 0].set_ylim(-0.3, 0.7)
    ax[3, 0].set_xlim(-0.2, 0.5)


    # Plot mg24 vs Ba
    print(f"mg24 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], iso_mg['MgFe24'])}")
    ax[1, 1].errorbar(Ba_values['[Ba/Fe]'], iso_mg['MgFe24'], yerr=iso_mg['d_MgFe24'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    # ax[0, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg24'], yerr=small_mg['d_mg24'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[1, 1].set_ylabel('[$^{24}$Mg/Fe]', fontsize=12)
    ax[1, 1].set_xlabel('[Ba/Fe]', fontsize=12)
    ax[1, 1].set_ylim(-0.2, 0.05)


    # Plot mg25 vs Ba
    print(f"mg25 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], iso_mg['MgFe25'])}")
    ax[2,1].errorbar(Ba_values['[Ba/Fe]'], iso_mg['MgFe25'], yerr=iso_mg['d_MgFe25'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    # ax[1,1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg25'], yerr=small_mg['d_mg25'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[2,1].set_xlabel('[Ba/Fe]', fontsize=12)
    ax[2,1].set_ylabel('[$^{25}$Mg/Fe]', fontsize=12)
    # ax[2, 1].set_ylim(-0.5, 0.5)


    # Plot mg26 vs Ba
    print(f"mg26 vs Ba: {scipy.stats.pearsonr(Ba_values['[Ba/Fe]'], iso_mg['MgFe26'])}")
    ax[3, 1].errorbar(Ba_values['[Ba/Fe]'], iso_mg['MgFe26'], yerr=iso_mg['d_MgFe26'], xerr=Ba_values['e[Ba/Fe]'], fmt='o', elinewidth=0.5)
    # ax[2, 1].errorbar(small_ba_values['[Ba/H]'], small_mg['mg26'], yerr=small_mg['d_mg26'], xerr=small_ba_values['e[Ba/H]'], fmt='o', elinewidth=0.5, color='orange')
    ax[3, 1].set_xlabel('[Ba/Fe]', fontsize=12)
    ax[3, 1].set_ylabel('[$^{26}$Mg/Fe]', fontsize=12)
    ax[3, 1].set_ylim(-0.3, 0.7)


    # Save the plot
    plt.savefig(f'/home/users/qai11/Documents/Isotope-Pipeline/Paper_Figures/Results/Element_fits_X_Fe_new_calc.png', dpi=300, bbox_inches='tight')

element_plots_XFe_new(star_list)
# %%
