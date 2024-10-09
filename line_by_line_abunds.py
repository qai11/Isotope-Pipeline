"""
Title: line_by_line_abunds.py
Author: Heather Sinclair-Wentworth, edits by Quin Aicken Davies
Date: 06/10/24

Description: Find the abundances of elements by line-by-line analysis using iSpec
"""
#%%
import os
import sys
import numpy as np
import logging
import multiprocessing
from multiprocessing import Pool
import numpy as np
from astropy.io import fits
import pandas as pd
import matplotlib.pyplot as plt
import time
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


#%%
start=time.time()


def params_abun(folder, star, para, code = "MOOG", line_by_line = False, element = None, filename=None,):
    # star_spectrum = ispec.read_spectrum("/home/users/hsi63/Masters/Spectra/" + filename)    
    # normalized_star_spectrum, name = get_final_spec(filename)
    # print(target_name)
    #rv and median stacked spectrum
    star_name = star
    # normalized_star_spectrum = ispec.read_spectrum("/home/users/qai11/Documents/Fixed_fits_files/{star_name}/rv_corrected/median_spectrum_{star_name}.txt")
    normalized_star_spectrum = ispec.read_spectrum(f'/home/users/qai11/Documents/Fixed_fits_files/{star_name}/{star_name}_adjusted.fits')
    parameters = pd.read_csv(f'/home/users/qai11/Documents/Fixed_fits_files/iteration_parameters/{star_name}_iter_param.txt')
    
    # with fits.open("/home/users/hsi63/Masters/Spectra_edited/" + star_name + ".fits") as hdul:
    #     header = hdul[0].header
    # # Use a fixed value because the spectrum is already normalized
    star_continuum_model = ispec.fit_continuum(normalized_star_spectrum, fixed_value=1.0, model="Fixed value")

    # Parameters, change as needed
#####
    teff = parameters['teff'][len(parameters)-1]
    logg = parameters['logg'][len(parameters)-1]
    MH = parameters['MH'][len(parameters)-1]
    Vmic = parameters['Vmic'][len(parameters)-1]
    Vmac = parameters['Vmac'][len(parameters)-1]
    alpha = parameters['alpha'][len(parameters)-1]
    vsini = parameters['Vsini'][len(parameters)-1]
    R=82000
    vrad = 0
    limb_darkening_coeff = 0.6
    max_iterations = 10

    # Selected model amtosphere, linelist and solar abundances
    model = ispec_dir + "/input/atmospheres/MARCS.GES/"

    #Get linelist
    atomic_linelist_file = ispec_dir + '/input/linelists/transitions/Quin_GES_LIST.420_920nm/atomic_lines.tsv'


    #Get solar abundances
    if "ATLAS" in model:
        solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.1998/stdatom.dat"
    else:
        # MARCS
        # solar_abundances_file = ispec_dir + "/input/abundances/Kurucz/stdatom.dat"
        solar_abundances_file = ispec_dir + "/input/abundances/Grevesse.2007/stdatom.dat"


    #Get isotopes file
    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"

    # Load chemical information and linelist
    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=np.min(normalized_star_spectrum['waveobs']), wave_top=np.max(normalized_star_spectrum['waveobs']))

    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01] # Select lines that have some minimal contribution in the sun

    isotopes = ispec.read_isotope_data(isotope_file)


    # Load model atmospheres
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)

    # Load SPECTRUM abundances
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)


    # Free parameters
    # free_params = ["Teff", "logg", "MH", "vmic", "vmac", "alpha", "vsini"]
    free_params = ["vmic", "vmac", "vsini", "alpha"]

    # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)

#####? not 100% sure what this is doing
    chemical_elements_file = ispec_dir + "/input/abundances/chemical_elements_symbols.dat"
#####?

    chemical_elements = ispec.read_chemical_elements(chemical_elements_file)

######################################################################################

#Regions around lines, same as plotting them in GUI, can be on or off. Need on for line by line abundances

    line_regions = None
        
    #No line regions seems to be better for teff and logg
#########################################################################################

    #If doing line by line abundances
    issue = []
    df0 = pd.DataFrame()
    
    if element is None:
        element = "all"
        line_regions = ispec.read_line_regions("/home/users/qai11/Documents/quin-masters-code/Linelists/all_lines.csv")
    else:
        if len(element) == 1:
            line_regions = ispec.read_line_regions(f"/home/users/qai11/Documents/quin-masters-code/Linelists/{element} _lines.csv")
        else:
            line_regions = ispec.read_line_regions(f"/home/users/qai11/Documents/quin-masters-code/Linelists/{element}_lines.csv")
    
    #Use fixed given params
    free_params = []
    
    #Check line regions, make sure they lie within spectral range
    line_regions = line_regions[(line_regions['wave_peak'] > normalized_star_spectrum['waveobs'][0]) & (line_regions['wave_peak'] < normalized_star_spectrum['waveobs'][-1])]

    while normalized_star_spectrum["flux"][0] == 0:
        normalized_star_spectrum = normalized_star_spectrum[1:]
    while normalized_star_spectrum["flux"][-1] == 0:
        normalized_star_spectrum = normalized_star_spectrum[:-1]
        
    print(line_regions)
    # print(star_spectrum["waveobs"][0], star_spectrum["waveobs"][-1])
    # normalised_star_spectrum = normalized_star_spectrum[(normalized_star_spectrum["waveobs"][0] > 490.5) & (normalized_star_spectrum["waveobs"][-1] < 919.5)]

    
    # filter = ispec.create_wavelength_filter(normalized_star_spectrum, wave_base=500.0, wave_top=675.0)
    # normalized_star_spectrum = normalized_star_spectrum[filter]
    for i, line in enumerate(line_regions):
        mar_val = 0.25
        # try:
        # if line['wave_peak'].between():
        # Directory and file names
        #element_name = "_".join(line['element'].split())
        element_name = "_".join(line['note'].split())
        common_filename = "example_" + code + "_individual_" + element_name + "_%.4f" % line['wave_peak']

        # Free individual element abundance (WARNING: it should be coherent with the selected line regions!)
        #free_abundances = ispec.create_free_abundances_structure([line['element'].split()[0]], chemical_elements, solar_abundances)
        free_abundances = ispec.create_free_abundances_structure([line['note'].split()[0]], chemical_elements, solar_abundances)
        free_abundances['Abund'] += MH # Scale to initial metallicity

        linelist_free_loggf = None

        # Line by line
        individual_line_regions = line_regions[i:i+1] # Keep recarray structure

        # Segment
        segments = ispec.create_segments_around_lines(individual_line_regions, margin=mar_val)

        wfilter = ispec.create_wavelength_filter(normalized_star_spectrum, regions=segments) # Only use the segment

        if len(normalized_star_spectrum[wfilter]) == 0 or np.any(normalized_star_spectrum['flux'][wfilter] == 0):
            continue
        # print(normalized_star_spectrum[wfilter])
        #free_abundances = None
        
        if float(line["wave_peak"]) > 890:
            dif = abs(pd.DataFrame(normalized_star_spectrum["waveobs"]) - float(line["wave_peak"])).idxmin().values
            line_flux = normalized_star_spectrum["flux"][dif]
            
            if line_flux > 0.975:
                continue
        
        try:
            obs_spec, modeled_synth_spectrum, params, errors, abundances_found, loggf_found, status, stats_linemasks = \
            ispec.model_spectrum(normalized_star_spectrum[wfilter], star_continuum_model, \
            modeled_layers_pack, atomic_linelist, isotopes, solar_abundances, free_abundances, linelist_free_loggf, teff, \
            logg, MH, alpha, Vmic, Vmac, vsini, \
            limb_darkening_coeff, R, vrad, free_params, segments=segments, \
            linemasks=individual_line_regions, \
            enhance_abundances=False, \
            use_errors = True, \
            vmic_from_empirical_relation = False, \
            vmac_from_empirical_relation = False, \
            max_iterations=max_iterations, \
            tmp_dir = None, \
            code=code)
            a = np.transpose(abundances_found)
            df4 = pd.DataFrame(data = a)
            df4["wave_peak"] = individual_line_regions["wave_peak"]
            df4["niter"] = status["niter"]
            df4["note"] = individual_line_regions["note"]
            df0 = pd.concat([df0, df4], ignore_index = True)
            logging.info("Saving synthetic spectrum...")  
        except:
            # print("Issue star")
            print(star)
            # issue.append(target_name)
            
        # except:
        #     print("evil line", line)

    # print("There were ", len(issue), "stars")
    # print(issue)
    return (df0, issue, normalized_star_spectrum.meta)



def run_abunds(t,element):
    # targets = ["NGC362", "NGC104", "NGC1261", "NGC1851", "NGC1904", "NGC2808", "NGC4372", "NGC4590", "NGC5927", "M12", "NGC6752", "M15", "M2", "NGC4833"]
    #Careful with single letter ones, they have a space after
    ele = element
    for el in ele:
        
        # star_list = pd.read_csv(f"/home/users/hsi63/Masters/Finals/{t}_final.csv", sep=',')
        # star_list = pd.read_csv("/home/users/hsi63/Masters/Finals/problem.csv", sep=',')
        star_list = t
        
        star_list = star_list[star_list["TEFF"] > 0]
        star_list = star_list.reset_index(drop=True)

        final_abunds = pd.DataFrame()
        final_abunds["OBJECT"] = star_list["OBJECT"]
        final_abunds["TEFF"] = star_list["TEFF"]
        final_abunds["LOGG"] = star_list["LOGG"]
        final_abunds["FEH"] = star_list["FEH"]
        final_abunds["Ra"] = star_list["RA"]
        final_abunds["Dec"] = star_list["DECLINATION"]
        final_abunds[f"Ref_{el}1"] = star_list[f"{el.upper()}1"]
        
        if len(el) == 1:
            line_list= pd.read_csv(f"/home/users/qai11/Documents/quin-masters-code/Linelists/{el} _lines.csv", sep='\t')
        else:
            line_list= pd.read_csv(f"/home/users/qai11/Documents/quin-masters-code/Linelists/{el}_lines.csv", sep='\t')
        for line in line_list["wave_peak"]:
            final_abunds[f"{line}_abund"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_abund_err"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_A(X)"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_A(X)_err"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_X/H"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_X/H_err"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_X/Fe"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan
            final_abunds[f"{line}_X/Fe_err"] = np.zeros(len(final_abunds["OBJECT"])) * np.nan

        # print(final_abunds.columns)
        al = []
        al_e = []
        bad_stars = []

        for i, star in enumerate(star_list["OBJECT"]):
            # file = ast.literal_eval(star_list["File_list"][i])
            dataframe, iss, header = params_abun("/home/users/qai11/Documents/quin-masters-code", star, [star_list["TEFF"][i], star_list["LOGG"][i], star_list["FEH"][i], 0.4, 1, 4, 2], line_by_line=True, element= el)
            # print(dataframe.columns)
            # print(dataframe["A(X)"])
            if len(dataframe) != 0:
                al.append(np.median(dataframe["A(X)"]))
                al_e.append(np.median(dataframe["eA(X)"]))
                for j in range(len(dataframe)):
                    if str(dataframe["wave_peak"][j]) + "_abund" in list(final_abunds.columns):
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_abund"] = dataframe["Abund"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_abund_err"] = dataframe["eAbund"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_A(X)"] = dataframe["A(X)"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_A(X)_err"] = dataframe["eA(X)"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_X/H"] = dataframe["[X/H]"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_X/H_err"] = dataframe["e[X/H]"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_X/Fe"] = dataframe["[X/Fe]"][j]
                        final_abunds.loc[i, f"{dataframe['wave_peak'][j]}_X/Fe_err"] = dataframe["e[X/Fe]"][j]
            else:
                al.append(np.nan)
                al_e.append(np.nan)
            bad_stars.append(iss)
        print(bad_stars)
        final_abunds[el] = al
        final_abunds[f"{el}_e"] = al_e
        # print(header)
        # final_abunds["SNR"] = header["SNR"]
        # final_abunds["RES"] = header["SPEC_RES"]
        # final_abunds["GES_FLD"] = star_list["GES_FLD"]

        final_abunds.to_csv(f"/home/users/qai11/Documents/Fixed_fits_files/{star}/{el}/{t}_{el}_abunds.csv")


#Make work :)))) 
# elements = ["Mg", "Al"]


# def make_final(data):
#     targets = ["NGC362", "NGC104", "NGC1261", "NGC1851", "NGC1904", "NGC2808", "NGC4372", "NGC4590", "NGC5927", "M12", "NGC6752", "M15", "M2", "NGC4833"]
#     for el in elements:
#         al = []
#         for t in targets:
#             single = pd.read_csv(f"/home/users/hsi63/Masters/abunds/{el}/{t}_{el}_abunds.csv", sep = ",", index_col=0)
#             al.append(single)
#             single.to_csv(f"/home/users/hsi63/Masters/abunds/{el}/{t}_{el}_abunds.csv", sep = ",")
#         combine = pd.concat(al)
#         combine.to_csv(f"/home/users/hsi63/Masters/abunds/{el}/All_{el}_abunds.csv", sep = ",", index = False)

# targets = ["NGC362", "NGC104", "NGC1261", "NGC1851", "NGC1904", "NGC2808", "NGC4372", "NGC4590", "NGC5927", "M12", "NGC6752", "M15", "M2", "NGC4833"]
# targets = ["M15"]
# targets = ["NGC4833"]
# targets = ["NGC6752"]


end = time.time()
print(f'It took {end-start} seconds to run.')


# %%
