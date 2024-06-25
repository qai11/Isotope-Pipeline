"""
Title: Exp_calc.py
Author: Quin Aicken Davies
Date: 19/06/2024

Description: calculates an exposure time for 1 metre fibre 1
"""

# %%
import numpy as np

def magnitude_to_photon_flux(magnitude, zero_point_flux, wavelength):
    # Convert zero point flux from W/m^2/nm to photons/s/m^2/nm
    h = 6.626e-34  # Planck's constant (Joule*seconds)
    c = 3e8        # Speed of light (meters/second)
    energy_per_photon = h * c / wavelength  # Energy per photon (Joules)
    zero_point_flux_photons = zero_point_flux / energy_per_photon  # Photons/s/m^2/nm
    return zero_point_flux_photons * 10**(-magnitude / 2.5)

def exposure_time(SNR, wavelength, area, quantum_efficiency, spectrograph_efficiency, photon_flux):
    # Calculate energy per photon
    h = 6.626e-34  # Planck's constant (Joule*seconds)
    c = 3e8        # Speed of light (meters/second)
    energy_per_photon = h * c / wavelength
    
    # Calculate exposure time
    t = (SNR**2 * energy_per_photon) / (area * quantum_efficiency * spectrograph_efficiency * photon_flux)
    return t

# Inputs
SNR = 100
wavelength = 500e-9  # 500 nm in meters
area = np.pi * (0.5)**2  # Effective area of a 1 meter telescope (m^2)
quantum_efficiency = 0.73
spectrograph_efficiency = 0.188
magnitude = 3.4  # Example magnitude of the source

# Zero-point flux (W/m^2/nm) for magnitude 0 star in the V-band
zero_point_flux = 3.6e-8

# Convert magnitude to photon flux
photon_flux = magnitude_to_photon_flux(magnitude, zero_point_flux, wavelength)

# Calculate exposure time
exposure_time_seconds = exposure_time(SNR, wavelength, area, quantum_efficiency, spectrograph_efficiency, photon_flux)

print(f"Exposure time: {exposure_time_seconds} seconds")

# %%
