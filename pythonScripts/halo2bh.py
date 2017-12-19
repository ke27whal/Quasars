### Purpose: To calculate the black hole masses for a given set of DM halo masses
### res1, res2 = halo2bhmass(mhalo)

### input:
### mhalo - halo masses (in M_sun/h)
### hmf - the output of halo_mass_function

### output:
### bhmass -  black hole masses (in M_sun/h)
### bhmf - black hole mass function

### Written by: Kelly Whalen (Dartmouth College)

def halo2bh(mhalo, hmf):
    import numpy as np
    
    # Halo Mass to Galaxy mass

    c = 0.129
    M0 = 10.**11.4
    alpha = 0.926
    beta = 0.261
    gamma = 2.440

    galaxySample = c*((mhalo/M0)**(-1*alpha) + (mhalo/M0)**(beta))**(-1*gamma)

    stellarmassSample = galaxySample*mhalo

    logSM = np.log10(stellarmassSample)
    logHaloMass = np.log10(stellarmassSample/galaxySample)

    convert = logHaloMass/logSM

    GMF = hmf * convert

    ### Galaxy Mass to Black Hole Mass

    bhmass = 10**(8.20 + 1.12*np.log10(stellarmassSample/1.e11))

    logMbh = 8.2 + 1.12*(np.log10(stellarmassSample) - 11)
    logMstellarSample = ((logMbh - 8.2)/1.12) + 11

    startobhSample = logMstellarSample/logMbh

    bhmf = GMF * startobhSample

    return bhmass, bhmf
