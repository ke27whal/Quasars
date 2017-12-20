### Purpose: To calculate the bias for a given set of halo masses
### bias = mhalo2bias(mhalo, z)

### input:
### mhalo - halo masses (in M_sun/h)
### z - redshift

### output:
### bias

### Tranferred to Python by: Kelly Whalen (Dartmouth College)

### Note: This code is explicitly for Tinker et al. 2010 cosmology- this can be changed in sigma_m.py script

def mhalo2bias(mhalo, z):

    # importing necessary packages
    import numpy as np
    from sigma_m import sigma_m

    # set up parameters
    h0 = 0.702
    omega_m = 0.275
    omega_b = 0.046
    omega_l = 0.725
    spec_ind = 0.96
    Delta=200.

    # Set Hubble parameter scaling E(z)
    E_z=np.sqrt((omega_m*(1.+z)**3.) + omega_l)

    # Scale omega_m with z
    omega_m_z=omega_m*((1.+z)**3.)/(E_z**2.)

    # Get sigma_m
    sig_m = sigma_m(mhalo, z)

    #Use approximation of NFW 97 for delta (valid in universe with
    #Lambda, while Tinker delta_c=1.69 only for Omega_m=1)

    delta_c=0.15*((12.*np.pi)**(2./3.))*((omega_m_z)**(0.0055))

    nu=delta_c/sig_m

    y=np.log10(Delta)
    Abig=1.+0.24*y*np.exp(-(4/y)**4.)
    asmall=0.44*y-0.88
    Bbig=0.183
    bsmall=1.5
    Cbig=0.019+0.107*y+0.19*np.exp(-(4/y)**4.)
    csmall=2.4

    term1 = Abig*((nu**asmall)/((nu**asmall)+(delta_c**asmall)))
    term2 = Bbig*(nu**bsmall)
    term3 = Cbig*(nu**csmall)
   
    bias = 1. - term1 + term2 + term3

    return bias
