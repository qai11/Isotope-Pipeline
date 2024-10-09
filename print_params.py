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
star = ['hd_45588','hd_100407','hd_102870','hd_128620','hd_128621','hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
# star = ['hd_11695','hd_146233','hd_156098','hd_157244','hd_160691','moon']
#%%
params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_params.txt')
errors = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/parameters/final_errors.txt')

# %%
iteration = 0
for i in star:
    print('\n')
    print(star[iteration])
    print(f'teff = {params["teff"][iteration]:.0f} pm {errors["teff"][iteration]:.0f}')
    print(f'logg = {params["logg"][iteration]:.2f} pm {errors["logg"][iteration]:.2f}')
    print(f'MH = {params["MH"][iteration]:.2f} pm {errors["MH"][iteration]:.2f}')
    print(f'vmac = {params["vmac"][iteration]:.2f} pm {errors["vmac"][iteration]:.2f}')
    print(f'vsini = {params["vsini"][iteration]:.2f} pm {errors["vsini"][iteration]:.2f}')
    iteration +=1
    
# %%
#plot the convergence of the parameters to the final values for each star
for i,star_name in enumerate(star):
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
star = ['hd_45588']
for i,star_name in enumerate(star):
    plot_params = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    fig.tight_layout(pad=3)  # Add padding between subplots
    axs[0, 0].plot(plot_params['teff'],label='teff')
    axs[0, 0].set_ylabel('Teff (K)')
    
    axs[0, 1].plot(plot_params['logg'],label='logg')
    axs[0, 1].set_ylabel('logg')
    
    axs[1, 0].plot(plot_params['MH'],label='MH')
    axs[1, 0].set_xlabel('Iteration')
    axs[1, 0].set_ylabel('[M/H]')
    
    axs[1, 1].plot(plot_params['CHI-SQUARE'],label='$\chi^2$')
    axs[1, 1].set_xlabel('Iteration')
    axs[1, 1].set_ylabel('$\chi^2$')  
    
plt.savefig(f'/home/users/qai11/Documents/Masters_Figures/Method/{star_name}_convergence.png',dpi=150)
    
    
    
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

spectra_names = ["synth_only_Mg_5183.txt","cut.txt","synth_NoMg.txt","synth_only_Mg_5167.322.txt","synth_only_Mg_5172"]

file1 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_synth_only_Mg_5183.txt'
file2 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_cut.txt'
file3 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_synth_NoMg.txt'
file4 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_synth_only_Mg_5167.txt'
file5 = '/home/users/qai11/Documents/Fixed_fits_files/hd_102870/J0354027_synth_only_Mg_5172.txt'

data1 = pd.read_csv(file1,delimiter='	')
data2 = pd.read_csv(file2,delimiter='	')
data3 = pd.read_csv(file3,delimiter='	')
data4 = pd.read_csv(file4,delimiter='	')
data5 = pd.read_csv(file5,delimiter='	')

plt.figure(figsize=(10,5))
plt.plot(data1['waveobs'],data1['flux'],label='J0354027_synth_only_Mg_5183',color='darkblue')
plt.plot(data2['waveobs'],data2['flux'],label='J0354027_cut',color='maroon')
plt.plot(data3['waveobs'],data3['flux'],label='J0354027_synth_NoMg',color='purple')
plt.plot(data4['waveobs'],data4['flux'],label='J0354027_synth_only_Mg_5167',color='green')
plt.plot(data5['waveobs'],data5['flux'],label='J0354027_synth_only_Mg_5172',color='black')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Flux')
plt.xlim(516,519)
plt.ylim(0.1,1.05)
plt.legend(loc='lower right')

plt.savefig(f'/home/users/qai11/Documents/Masters_Figures/Method/All_Mg_overplot.png',dpi=150)
# %%

#For printing abundances
iteration = 0
for star_name in star:
    #Path for Mac
    abund_params = pd.read_csv(f'/Users/quin/Desktop/2024_Data/parameters_gbs_linelist_only_first_calc/{star_name}_final_indiv_abund.txt')
    print('\n')
    print(star[iteration])
    for i in range(0,len(abund_params)):
        print(f'{abund_params["element"][i]} = {abund_params["[X/H]"][i]:.5f} pm {abund_params["e[X/H]"][i]:.5f}')
    iteration +=1
# %%
import pandas as pd

table_data = []
star = ['hd_45588', 'hd_100407', 'hd_102870', 'hd_128620', 'hd_128621', 'hd_11695', 'hd_146233', 'hd_156098', 'hd_157244', 'hd_160691', 'moon']

# Initialize an empty dictionary to hold abundance data for each star
abundance_dict = {}

for star_name in star:
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
