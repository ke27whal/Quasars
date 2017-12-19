### Based on sigma_.pro written by Mike DiPompeo
###  Calculate RMS of mass density field for given halo mass/redshift.

### Usage: sigm = sigma_m(mhalos, z)

### input:
### mhalos - halo masses in (in M_sun/h)
### z - redshift

### output:
### sigma_M

### Tranferred to Python by: Kelly Whalen (Dartmouth College)

def sigma_m(mhalo, z):
    
    import sys, platform, os
    sys.path.insert(0, '/Users/kellywhalen/Downloads/CAMB-0.1.6.1/pycamb') # Your current path
    import camb
    from camb import model, initialpower
    import numpy as np
    
    #Setting up a few parameters
    h0 = 0.702
    omega_m = 0.275
    omega_l = 0.725
    omega_b = 0.046
    omega_bh2 = omega_b*h0**2
    omega_ch2 = (omega_m - omega_b)*(h0**2.)
    spec_index = ns=0.96
    rho_crit=2.7745e11
    rho_mean=rho_crit*omega_m
    
    

    #Get power spectrum using CAMB
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=h0*100., ombh2=omega_bh2, omch2=omega_ch2)
    pars.set_dark_energy() #re-set defaults
    pars.InitPower.set_params(ns=spec_index, As = 2.46e-9)
    
    
    #Not non-linear corrections couples to smaller scales than you want
    pars.set_matter_power(redshifts= [z], kmax=10.0)
    
   #Non-Linear spectra (Halofit)
    pars.NonLinear = model.NonLinear_both
    results = camb.get_results(pars)
    results.calc_power_spectra(pars)
    wavenum, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=50, npoints = 200)

    pk = pk[0]
    Rf=((3./4.)*(mhalo/(np.pi*rho_mean)))**(1./3)
    sigm = np.zeros(len(mhalo))
 
    for i in range(len(mhalo)):
        WHat=(3./((wavenum*Rf[i])**3.))*(np.sin(wavenum*Rf[i])-((wavenum*Rf[i])*np.cos(wavenum*Rf[i])))
        integ=pk*(WHat**2.)*(wavenum**2.)
        res = np.trapz(integ, wavenum)
        sigm[i] = ((1./(2.*(np.pi**2.)))*res)**(1./2.) 
 
    return sigm
