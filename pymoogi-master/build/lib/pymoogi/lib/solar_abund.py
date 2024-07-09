# This is a list of solar abundaces copied from moog/Batom.f
# Zero in added for first index so indexes represent atomic numbers

  #  xabu = the set of current solar (when available) or meteorite       
  #  abundances, scaled to log(h) = 12.00 .  The data are from Asplund
  #  et al. (2009, Ann. Rev. Ast. Ap., 47, 481).

def get_solar_abund():

    solar = [0, 
        12.00,10.93, 1.05, 1.38, 2.70, 8.43, 7.83, 8.69, 4.56, 7.93,    
        6.24, 7.60, 6.45, 7.51, 5.41, 7.12, 5.50, 6.40, 5.03, 6.34,     
        3.15, 4.95, 3.93, 5.64, 5.43, 7.50, 4.99, 6.22, 4.19, 4.56,     
        3.04, 3.65, 2.30, 3.34, 2.54, 3.25, 2.52, 2.87, 2.21, 2.58,    
        1.46, 1.88,-5.00, 1.75, 0.91, 1.57, 0.94, 1.71, 0.80, 2.04,     
        1.01, 2.18, 1.55, 2.24, 1.08, 2.18, 1.10, 1.58, 0.72, 1.42,     
        -5.00, 0.96, 0.52, 1.07, 0.30, 1.10, 0.48, 0.92, 0.10, 0.84,     
        0.10, 0.85,-0.12, 0.85, 0.26, 1.40, 1.38, 1.62, 0.92, 1.17,     
        0.90, 1.75, 0.65,-5.00,-5.00,-5.00,-5.00,-5.00,-5.00, 0.02,     
        -5.00,-0.54,-5.00,-5.00,-5.00]

    return solar


