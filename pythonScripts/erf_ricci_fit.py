### Purpose: The goal of this program is to fit an error function to the empirical relationship between obscured fraction and Eddington ratio as seen in Ricci et al. 2017. This function also has the option of calculating the obscured fraction for a simulated population of quasars, as well as the median halo masses for the obscured and unobscured populations of quasars generated using the error function fit to the Ricci et al. 2017 relationship.               


### erf, ledd_plot, obsc_med, unobsc_med, obscFraction = erf_ricci_fit(min_fobs, w, x_shift, mhalo_lum, ledd_lum, randos_lum = None, weights_lum = None)

### input:
### min_fobs -The minimum obscured fraction out at high Eddington ratios
### w - ''width parameter,'' adjusts the range over which the erf is fit, larger number results in a steeper slope at mid Eddington
### x_shift - A parameter that moves the inflection point of the erf fit, more positive value results in a shift to the right
###  mhalo_lum - A distribution of DM halo masses after a luminosity cut has been made
###  ledd_lum - A distribution of Eddington ratios after a luminosity cut has been made
### randos_lum - An optional parameter; an array of uniform random numbers 0<n<1, of len(mhalo_lum), used to identify obscured and unobscured quasars
### weights_lum - An optional parameter; an array of weights of len(mhalo_lum) to be used for the weighted medians

### output:
### erf -  error function fit, an array of len(mhalo_lum)
### ledd_plot - array of log(Eddington ratio) of len(erf) over which erf was fit
### obsc_med - weighted median DM halo mass for obscured sample [M_solar]
### unobsc_med - weighted median DM halo mass for unobscured sample [M_solar]
### obscFraction - obscured fraction of generated quasar population

### Written by: Kelly Whalen (Dartmouth College)

def erf_ricci_fit(min_fobs, w, x_shift, mhalo_lum, ledd_lum, randos_lum = None, weights_lum = None):

    # Import libraries and called functions
    import numpy as np
    from scipy.interpolate import interp1d
    from wmedian import wmedian
    import scipy

    # Set the bounds for the Ricci relationship to be fit
    x = min_fobs
    rlogledd = [-4.0, -3.8, -3.0, -2.2, -1.7, -1.2, -0.8, 0.0, 1.0]
    rfobs = [0.85, 0.85, 0.85, 0.8, 0.70, 0.6, x, x, x]

    # Initializes Eddington ratio range over which erf will fit Ricci relationship 
    interp_ledd = np.linspace(-4., 1.0, len(mhalo_lum))

    # Fitting the error function
    erf_xaxis_loop_x = np.linspace(-1, 1, len(mhalo_lum)) * w - w*x_shift
    erf_ycompress_loop_x = -((max(rfobs) - min(rfobs))/2.)*scipy.special.erf(erf_xaxis_loop_x)
    erf_test_loop_x = erf_ycompress_loop_x + (max(rfobs) - max(erf_ycompress_loop_x))

    # Interpolating to the post-Lcut simulated data
    erf_test_interpF_loop_x = interp1d(sorted(interp_ledd), erf_test_loop_x, bounds_error = False)
    erf_test_interpolated_loop_x = erf_test_interpF_loop_x(interp_ledd)
    
    lum_fobs_interp_f = interp1d(sorted(interp_ledd), erf_test_loop_x, bounds_error = False)
    lum_fobs_interp = lum_fobs_interp_f(np.log10(ledd_lum))
    
    # Calculate median halo masses and obscured fraction if input parameters are present
    if len(randos_lum) != 0:
        obscFlag_erf_loop_x = (randos_lum < lum_fobs_interp) 
        unobscFlag_erf_loop_x = (randos_lum >= lum_fobs_interp)
    
        avg_obsc_x = wmedian(mhalo_lum[obscFlag_erf_loop_x], weights = weights_lum[obscFlag_erf_loop_x])
        avg_unobsc_x = wmedian(mhalo_lum[unobscFlag_erf_loop_x], weights = weights_lum[unobscFlag_erf_loop_x])

        obscuredFraction = np.float(sum(obscFlag_erf_loop_x)) / np.float(sum(obscFlag_erf_loop_x) + sum(unobscFlag_erf_loop_x))
        
        return erf_test_interpolated_loop_x, interp_ledd, avg_obsc_x, avg_unobsc_x, obscuredFraction
    
    else: 
        return erf_test_interpolated_loop_x, interp_ledd
