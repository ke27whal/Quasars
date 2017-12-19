### Based on halo_mass_function_.pro written by Mike DiPompeo
###  Get N(M) for a halo mass (or masses) using Tinker et al. 2010 cosmology

### Usage: res=halo_mass_function(log_mhalo, z)

### input:
### mhalos - halo masses in (in M_sun/h)
### z - redshift

### output:
### n_mhalo

### Tranferred to Python by: Kelly Whalen (Dartmouth College)


def halo_mass_function(mhalo, z):
    import numpy as np
    from scipy.interpolate import interp1d
    from sigma_m import sigma_m
    
    delta = 200.
    
    h0 = 0.702
    omega_m = 0.275
    omega_l = 0.725
    omega_b = 0.046
    omega_bh2 = omega_b*h0**2
    omega_ch2 = (omega_m - omega_b)*(h0**2.)
    spec_index = ns=0.96
    rho_crit=2.7745e11
    rho_mean=rho_crit*omega_m
    
    # Set Hubble parameter scaling E(z)
    E_z = np.sqrt((omega_m*(1. + z)**3) + omega_l)
    
    # Scale omega_m, omega_l with z
    omega_m_z=omega_m*((1.+z)**3.)/(E_z**2.)
    
    # Get sigma_m
    sig_m = sigma_m(mhalo, z)
    
    # Use approximation of NFW 97 for delta (valid in universe with Lambda, while Tinker delta_c=1.69 only for 
    # Omega_m=1), 
    delta_c = 0.15*((12*np.pi)**(2./3.))*((omega_m_z)**(0.0055))
    nu = delta_c/sig_m
    
    # Tinker 2010 model
    Deltas=[200.,300.,400.,600.,800.,1200.,1600.,2400.,3200.]
    alphas=[0.368,0.363,0.385,0.389,0.393,0.365,0.379,0.355,0.327]
    betas=[0.589,0.585,0.544,0.543,0.564,0.623,0.637,0.673,0.702]
    gams=[0.864,0.922,0.987,1.09,1.20,1.34,1.50,1.68,1.81]
    phis=[-0.729,-0.789,-0.910,-1.05,-1.20,-1.26,-1.45,-1.50,-1.49]
    etas=[-0.243,-0.261,-0.261,-0.273,-0.278,-0.301,-0.301,-0.319,-0.336]
    
    f1 = interp1d(Deltas, alphas, kind = 'quadratic')
    alpha = f1(delta)
    
    f2 = interp1d(Deltas, betas, kind = 'quadratic')
    beta = f2(delta)
    
    f3 = interp1d(Deltas, gams, kind = 'quadratic')
    gam = f3(delta)
    
    f4 = interp1d(Deltas, phis, kind = 'quadratic')
    phi = f4(delta)
    
    f5 = interp1d(Deltas, etas, kind = 'quadratic')
    eta = f5(delta)
    
    term1=1+((beta*nu)**((-2)*phi))
    term2=nu**(2.*eta)
    term3=np.exp(((-1.)*gam*(nu**2.))/2.)
    f=alpha*term1*term2*term3

    # Get N(M) in physical units
    rho_crit=2.7745e11
    rho_matter=rho_crit*omega_m

    log_mhalo = np.log10(mhalo)
   
    d_sig_M=sigma_m(10**(log_mhalo+0.001), z)
    d_lognu=np.log10(delta_c/d_sig_M)-np.log10(nu)
    dlognu_dlogM=d_lognu/0.001

    n_mhalo=(((nu*f*rho_matter*dlognu_dlogM)/(mhalo**2.))/0.7**4)*mhalo*np.log(10)

    return n_mhalo
