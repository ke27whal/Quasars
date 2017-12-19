### Purpose: To calculate the probability distribution for a range of Eddington ratios
### res = schechterLedd(Lambda)

### input:
### Lambda - sample of Eddington ratios

### output:
### schechter - probability distribution


def schechterLedd(Lambda):
    import numpy as np
    
    eta = 1.0
    alpha = 0.5
    t0 = 1.e8

    
    schech = t0 * (Lambda/eta)**(-alpha)*np.exp(-1*(Lambda))

    schechter = schech/np.sum(schech)

    return schechter
