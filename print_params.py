"""
Title: print_params.py
Author: Quin Aicken Davies
Date: 08/08/24

Description: Prints the parameters and errors for each star in the final_params.txt and final_errors.txt files.
"""
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621']
# star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
# star = ['hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
#%%
params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_params.txt')
errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_errors.txt')

# %%
'''print the parameters and errors for each star'''
iteration = 0
for star_name in star_list:
    print('\n')
    print(star_name)
    params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star_name}_final_params.txt', sep=',',index_col=0)
    errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star_name}_final_errors.txt', sep=',',index_col=0)
    print(f"{params['teff'].values[0]:.2f}")
    print(f"{errors['teff'].values[0]:.2f}")
    # print(f'teff = {params["teff"]} pm {errors["teff"]}')
    # print(f'logg = {params["logg"]} pm {errors["logg"]}')
    # print(f'MH = {params["MH"]} pm {errors["MH"]}')
    # print(f'vmac = {params["vmac"]} pm {errors["vmac"]}')
    # print(f'vsini = {params["vsini"]} pm {errors["vsini"]}')
    iteration +=1
#%%

# %%
#plot the convergence of the parameters to the final values for each star
for i,star_name in enumerate(star_list):
    plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
    # plot_params = plot_params.append(params.loc[i],ignore_index=True)
    plt.figure()
    plt.plot(plot_params['logg'],label='logg')
    plt.xlabel('Iteration')
    plt.ylabel('logg')
    plt.legend()
    plt.figure()
    plt.plot(plot_params['MH'],label='MH')
    plt.xlabel('Iteration')
    plt.ylabel('[M/H]')
    plt.legend()
    plt.figure()
    plt.plot(plot_params['teff'],label='teff')
    plt.xlabel('Iteration')
    plt.ylabel('Teff (K)')
    plt.legend()
    plt.figure()
    plt.plot(plot_params['CHI-SQUARE'],label='$\chi^2$')
    plt.xlabel('Iteration')
    plt.ylabel('$\chi^2$')
    plt.legend()

#%%
#plot the convergence of the parameters to the final values for each star onto sublplots
# ,teff,logg,MH,alpha,Vmic,Vmac,Vsini,limb_darkening_coeff,R,CHI-SQUARE,DOF
# star = ['hd_102870']
for i,star_name in enumerate(star_list):
    plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt',header=None)
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    fig.tight_layout(pad=3)  # Add padding between subplots
    axs[0, 0].plot(plot_params[1],label='teff')
    axs[0, 0].set_ylabel('Teff (K)')
    
    axs[0, 1].plot(plot_params[2],label='logg')
    axs[0, 1].set_ylabel('logg')
    
    axs[1, 0].plot(plot_params[3],label='MH')
    axs[1, 0].set_xlabel('Iteration')
    axs[1, 0].set_ylabel('[M/H]')
    
    axs[1, 1].plot(plot_params[10],label='$\chi^2$')
    axs[1, 1].set_xlabel('Iteration')
    axs[1, 1].set_ylabel('$\chi^2$')  
    
plt.savefig(f'/home/users/qai11/Documents/Masters_Figures/{star_name}_convergence.png',dpi=150)
    
    
    
#%%
'''Plots all the convergence values onto one plot'''
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['logg'],label='logg')
#     plt.xlabel('Iteration')
#     plt.ylabel('logg')
#     plt.legend(star)    
 
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['MH'],label='MH')
#     plt.xlabel('Iteration')
#     plt.ylabel('[M/H]')
#     plt.legend(star)  
    
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['teff'],label='teff')
#     plt.xlabel('Iteration')
#     plt.ylabel('Teff (K)')
#     plt.legend(star)  
    
# plt.figure()
# for i,star_name in enumerate(star):
#     plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
#     plt.plot(plot_params['CHI-SQUARE'],label='$\chi^2$')
#     plt.xlabel('Iteration')
#     plt.ylabel('$\chi^2$')
#     plt.legend(star)     
        

# %%
#Plot the Mg triplet removal from linelist effect on the synthetic spectra
import pandas as pd
import matplotlib.pyplot as plt

spectra_names = ["synth_only_Mg_5183.txt","cut.txt","synth_NoMg.txt","synth_only_Mg_5167.322.txt","synth_only_Mg_5172"]

file1 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/_GBS_nomol_nounc/J0354027_synth_only_Mg_5183.txt'
file2 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/_GBS_nomol_nounc/J0354027_cut.txt'
file3 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/_GBS_nomol_nounc/J0354027_synth_NoMg.txt'
file4 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/_GBS_nomol_nounc/J0354027_synth_only_Mg_5167.txt'
file5 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/_GBS_nomol_nounc/J0354027_synth_only_Mg_5172.txt'

data1 = pd.read_csv(file1,delimiter='	')
data2 = pd.read_csv(file2,delimiter='	')
data3 = pd.read_csv(file3,delimiter='	')
data4 = pd.read_csv(file4,delimiter='	')
data5 = pd.read_csv(file5,delimiter='	')

plt.figure(figsize=(10,5))
plt.plot(data1['waveobs']*10,data1['flux'],label='J0354027_synth_only_Mg_5183',color='darkblue')
plt.plot(data2['waveobs']*10,data2['flux'],label='J0354027_cut',color='maroon')
plt.plot(data3['waveobs']*10,data3['flux'],label='J0354027_synth_NoMg',color='purple')
plt.plot(data4['waveobs']*10,data4['flux'],label='J0354027_synth_only_Mg_5167',color='green')
plt.plot(data5['waveobs']*10,data5['flux'],label='J0354027_synth_only_Mg_5172',color='black')
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Flux')
plt.xlim(5170,5186)
plt.ylim(0.5,1.05)
plt.legend(loc='lower right')

plt.savefig(f'/home/users/qai11/Documents/quin-masters-code/Masters_Figures/Method/All_Mg_overplot.png',dpi=150,tight_layout=True)


# %%
"""Prints the abundances for each star"""
# star = ['hd_45588', 'hd_100407', 'hd_102870', 'hd_128620', 'hd_128621', 'hd_11695', 'hd_146233', 'hd_156098', 'hd_157244', 'hd_160691', 'moon']
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
#For printing abundances
iteration = 0
for star_name in star_list:
    #Path for Mac
    try:
        abund_params = pd.read_csv(f'/Users/quin/Desktop/2024_Data/parameters_gbs_linelist_only_first_calc/{star_name}_final_indiv_abund.txt')
    except:
        #Path for Uni
        abund_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/summary_abundances_{star_name}.txt', delimiter=' ')
    print('\n')
    # print(star[iteration])
    for i in range(0,len(abund_params)):
        print(f'{abund_params["element"][i]} = {abund_params["[X/H]"][i]:.5f} \pm {abund_params["e[X/H]"][i]:.5f}')
    iteration +=1
    
# %%
import pandas as pd

table_data = []
# star = ['hd_45588', 'hd_100407', 'hd_102870', 'hd_128620', 'hd_128621', 'hd_11695', 'hd_146233', 'hd_156098', 'hd_157244', 'hd_160691', 'moon']
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# Initialize an empty dictionary to hold abundance data for each star
abundance_dict = {}

for star_name in star_list:
    # Read the abundance data for each star
    abund_params = pd.read_csv(f'/Users/quin/Desktop/2024_Data/parameters_gbs_linelist_only_first_calc/{star_name}_final_indiv_abund.txt')
    
    # Initialize a temporary dictionary to store abundance values for this star
    temp_dict = {}
    
    for i in range(len(abund_params)):
        element = abund_params["element"][i]
        abundance = f'{abund_params["[X/H]"][i]:.3f} ± {abund_params["e[X/H]"][i]:.3f}'
        temp_dict[element] = abundance
    
    # Add the abundance data for this star to the main dictionary
    abundance_dict[star_name] = temp_dict

# Convert the abundance dictionary into a DataFrame
table_df = pd.DataFrame(abundance_dict)

# Transpose the DataFrame to get the star names as columns
table_df_transposed = table_df

# Convert DataFrame to LaTeX table
latex_table = table_df_transposed.to_latex(column_format='|c|' * (len(table_df_transposed.columns)+1))

# Manually insert \hline after each row
latex_table_with_hline = latex_table.replace("\\\\", "\\\\ \\hline").replace("\\toprule", "").replace("\\bottomrule", "").replace("\\midrule", "")

# Print the LaTeX table
print(latex_table_with_hline)

# %%

"""Create latex table for the parameters"""
import pandas as pd
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# Initialize an empty dictionary to hold parameter data for each star
parameter_dict = {}
for star in star_list:
    # Read the parameter data for each star
    params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star}_final_params.txt')
    # Read the error data for each star
    errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/{star}_final_errors.txt')
    # Initialize a temporary dictionary to store parameter values for this star
    temp_dict = {}
    temp_dict['Teff'] = f'{params["teff"].values[0]:.0f} ± {errors["teff"].values[0]:.0f}'
    temp_dict['logg'] = f'{params["logg"].values[0]:.2f} ± {errors["logg"].values[0]:.2f}'
    temp_dict['[M/H]'] = f'{params["MH"].values[0]:.2f} ± {errors["MH"].values[0]:.3f}'
    # Add the parameter data for this star to the main dictionary
    parameter_dict[star] = temp_dict
    # Convert the parameter dictionary into a DataFrame
table_df = pd.DataFrame(parameter_dict)
# Transpose the DataFrame to get the star names as columns
parameter_table_df_transposed = table_df.T
# Convert DataFrame to LaTeX table
#Have the top row as StarID,Teff,GBS Teff,logg,GBS logg,MH,GBS MH
latex_table = parameter_table_df_transposed.to_latex(column_format='|c|' * (len(parameter_table_df_transposed.columns)+1))

# Manually insert \hline after each row
latex_table_with_hline = latex_table.replace("\\\\", "\\\\ \\hline").replace("\\toprule", "").replace("\\bottomrule", "").replace("\\midrule", "")
# Print the LaTeX table
print(latex_table_with_hline)
    
# %%"""Caclulate the weighted averages of isotopic abundance ratios for each star"""

import pandas as pd
import numpy as np
import ast
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
# Initialize an empty dictionary to hold abundance data for each star
abundance_dict = {}

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
    iso_abund_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/par_unc_{star_name}.csv', delimiter=',', index_col=0)
    
    #open masters stars csv
    star_info = pd.read_csv(f'/home/users/qai11/Documents/quin-masters-code/Masters_stars.csv', sep=',')
    #get the star regions
    regions = star_info[star_info['ID2'] == star_name]['regions'].apply(ast.literal_eval).values[0]
    
    #only take the average of the given regions
    iso_abund_params = iso_abund_params[iso_abund_params['region'].isin(regions)]
    # Separate the error and parameter columns
    errors = iso_abund_params.filter(like='d_')
    params = iso_abund_params.drop(columns=errors.columns).drop(columns=['region', 'pass'])
    
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
    # MgH['[X/H]'] = MgH['[X/H]'] - feh
    # print(MgH)
    # print(MgH)
    mg_24,mg_25,mg_26 = calc_ratio(abundances[1], abundances[2], abundances[3])
    # print(mg_24)
    # print(abundances[4])
    # print(errors[4], errors[5], errors[6])
    
    #calculate the log abundances
    log_solar_mg = abundances[1]+7.53
    log_solar_i_24 = log_solar_mg * (mg_24/100)
    log_solar_i_25 = log_solar_mg * (mg_25/100)
    log_solar_i_26 = log_solar_mg * (mg_26/100)
    #convert to Mg/H
    mg_24_H_is = log_solar_i_24 - (7.53*0.7899)
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
    
#save the final_df to a csv file
# final_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/weighted_avg_iso_abund_paper.csv')



#%% '''Caculate the final abundances from the line by line abundances'''
import pandas as pd

# Define the list of stars and elements
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
# star_list = ['hd_11695']
element = ["Eu", "Ba", "Mg"]

# Iterate over each star and element
for star_name in star_list:
    # Create an empty dataframe with the desired columns
    new_abund_df = pd.DataFrame(columns=['element', 'code', 'Abund', 'A(X)', '[X/H]', '[X/Fe]', 'eAbund', 'eA(X)', 'e[X/H]', 'e[X/Fe]'])

    for el in element:
        # Read the data for the current star and element
        line_abund = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/line_by_line_{el}.txt', delimiter=' ')
        
        # If the element is Eu, copy the values directly into the dataframe
        if el == "Eu":
            new_abund_df = line_abund.copy()
            #remove first and last column
            new_abund_df = new_abund_df.drop(new_abund_df.columns[[0, -1]], axis=1)
        else:
            # If there are at least 3 data points, use median and mad, otherwise use mean and mad
            if len(line_abund) >= 3:
                median_XH = line_abund['[X/H]'].median()
                mad_XH = line_abund['[X/H]'].mad()
                median_XFe = line_abund['[X/Fe]'].median()
                mad_XFe = line_abund['[X/Fe]'].mad()
                median_Abund = line_abund['Abund'].median()
                mad_Abund = line_abund['Abund'].mad()
                median_A = line_abund['A(X)'].median()
                mad_A = line_abund['A(X)'].mad()
            else:
                median_XH = line_abund['[X/H]'].mean()
                mad_XH = line_abund['[X/H]'].mad()
                median_XFe = line_abund['[X/Fe]'].mean()
                mad_XFe = line_abund['[X/Fe]'].mad()
                median_Abund = line_abund['Abund'].mean()
                mad_Abund = line_abund['Abund'].mad()
                median_A = line_abund['A(X)'].mean()
                mad_A = line_abund['A(X)'].mad()
            
            # Extract the code for the current element
            code = line_abund['code'].iloc[0]
            
            # Append the calculated values to the dataframe
            new_abund_df = new_abund_df.append({
                'element': el, 
                'code': code, 
                'Abund': median_Abund, 
                'A(X)': median_A, 
                '[X/H]': median_XH, 
                '[X/Fe]': median_XFe, 
                'eAbund': mad_Abund, 
                'eA(X)': mad_A, 
                'e[X/H]': mad_XH, 
                'e[X/Fe]': mad_XFe
            }, ignore_index=True)

    # Save the dataframe to a txt file for the current star
    new_abund_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt', sep=' ', index=False)

    print(f"Data for {star_name} saved successfully.")
    # Put the information into a latex format
    print(new_abund_df.to_latex(index=False))
    
#%%'''Print the final abundances from the line by line abundances'''

import pandas as pd

# Define the list of stars and elements
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']
element = ["Eu", "Ba", "Mg"]


#Set the index for final_abund_df to be the starname
final_abund_df = pd.DataFrame(columns=['star_name','element', '[X/H]', 'e[X/H]'])

# Loop over each star
for star_name in star_list:
    # Read the data for the current star
    file_path = f'/home/users/qai11/Documents/Fixed_fits_files/lbl_abundances/{star_name}/good_lbl/summary_abundances_{star_name}.txt'
    line_abund = pd.read_csv(file_path, delimiter=' ')

    # Add the star name as a new column
    line_abund['star_name'] = star_name  

    # Select only the necessary columns
    line_abund = line_abund[['star_name', 'element', '[X/H]', 'e[X/H]',]]

    # Append the data using concat
    final_abund_df = pd.concat([final_abund_df, line_abund], ignore_index=True)

# Print the final abundances in LaTeX format
print(final_abund_df.to_latex(index=False))
    
    

# %% """Make a table for the isotopic ratios"""

import pandas as pd
import numpy as np
star_list = ['hd_11695','hd_18884','hd_157244','hd_18907','hd_22049','hd_23249','hd_128621',
    'hd_10700','hd_100407','hd_160691','moon','hd_128620','hd_146233','hd_165499','hd_2151',
    'hd_102870','hd_45588','hd_156098']

#Make an empty df to hold the isotope information
isotope_df = pd.DataFrame(columns=['star_name','s','mg','d_mg', 'i_24', 'i_25', 'i_26','R_24','R_25',
                                   'R_26','d_i_24', 'd_i_25', 'd_i_26','d_R_24','d_R_25','d_R_26',
                                   'pass','region','ratio'])
for star_name in star_list:
    # Read all the abundance data for the star (across all regions)
    iso_abund_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/moog_tests_paper/par_unc_{star_name}.csv', delimiter=',', index_col=0)
   
    # Add the information to the isotope_df
    for i in range(len(iso_abund_params)):
        new_row = {
            'star_name': star_name,
            's': iso_abund_params['s'][i],
            'd_s': iso_abund_params['d_s'][i],
            # 'mg': iso_abund_params['mg'][i],
            # 'd_mg': iso_abund_params['d_mg'][i],
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
isotope_df.to_csv(f'/home/users/qai11/Documents/Fixed_fits_files/All_isotope_ratios_pre_avg_vpass_{iso_abund_params['pass']}.csv')

# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the data from the csv file
synth_no_mg = pd.read_csv('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_synth_NoMg.txt', delimiter='	')
synth_only_mg_5183 = pd.read_csv('/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_synth_only_Mg_5183.txt', delimiter='	')
