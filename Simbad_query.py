"""
Created on Mon Mar 22 2024

@author: Quin Aicken Davies
"""
#%%
from astroquery.simbad import Simbad
import pandas as pd
from matplotlib import pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord, Angle
import pandas, numpy
import os


#%%
'''Import the Benchmark star list'''
star_list = pd.read_csv('Benchmark_Star_Names.csv')

'''Put the column of HD_Number into a list'''
HD_Number = star_list['HD_Number'].to_list()

'''Create file to query Simbad'''
bench_simbad = Simbad()

bench_simbad.TIMEOUT = 240

Simbad.list_votable_fields() #Lists what can you call

#%%
'''Input values you want to be queried'''
bench_simbad.add_votable_fields('flux(V)','fe_h','distance','pm','rot','rv_value','sptype')

# %%
'''Query Simbad'''
benchmark_stars = bench_simbad.query_objects(HD_Number)

# %%
'''Print the keys from the Query for Review'''
print(benchmark_stars.keys())

# %%

'''Plot all the positions of the data'''
# Load data
dataframe = benchmark_stars.to_pandas()
#%%
#prints the dataframe
#dataframe

#%% Processing coordinates
# Converting the units in RA & DEC (both string values, RA is also not in degrees) collums into degrees values
c = SkyCoord(
    dataframe["RA"].to_numpy(), 
    dataframe["DEC"].to_numpy(), 
    unit=(u.hourangle, u.deg)
)

# Converting to radians & fitting the RA units to 180 to -180
ra_radians = Angle(c.ra).wrap_at(180 * u.degree).radian
dec_radians = numpy.deg2rad(c.dec)

# Adding the radians back to the origianl dataframe
dataframe["dec_rad"] = dec_radians.value
dataframe["ra_rad"] = ra_radians

# Adding to the dataframe a filtered version
dataframe = dataframe[dataframe["dec_rad"] <= numpy.deg2rad(10)]


#%% Plotting Observed Targets
plt.figure(figsize=(10,10))
plt.subplot(111, projection="aitoff")
plt.scatter(dataframe["ra_rad"], dataframe["dec_rad"], s = 5, label=f"Not Observed: {len(dataframe)}")
plt.legend(loc="upper center")
plt.ylabel("DEC")
plt.xlabel("RA")
plt.grid()
plt.title("Isotope Targets", pad=20)
plt.savefig(os.path.dirname(os.path.abspath(__file__))+'/Targets')
plt.show()

#%% Plotting magnitudes
_, edges, _ = plt.hist(dataframe["FLUX_V"], bins=8, alpha=0.5, label="Total stars")

plt.legend()
plt.ylabel("total amount of stars per magnitude")
plt.xlabel("Visual magnitude")
plt.savefig("plots/Magnitude.pdf")
plt.show()

# %%

