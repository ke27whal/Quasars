### This code will calculate the weighted median for a given array and weights

def wmedian(data, weights):

    import numpy as np

    sweights = weights[np.argsort(data)]
    sdata = data[np.argsort(data)]

    totw = np.sum(sweights)

    cumulw = np.cumsum(sweights)

    medianindex = np.argwhere(cumulw < totw/2.)[-1]

    if len(sdata)%2 == 0:
        wmedian = (sdata[medianindex] + sdata[medianindex + 1])/2.
    else:
        wmedian = sdata[medianindex]

    return wmedian
