"""
Analyzing data in the HDF5 database
"""

eV = 0.03674932534

import numpy
from scipy.stats.kde import gaussian_kde
import sys

try:
    import tables
except ImportError:
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()



def Histogram(data, numbins = None, scale = 72):
    data = numpy.asanyarray(data)
    numdata = len(data)
    data_min = min(data)
    data_max = max(data)
    if numbins == None:
        numbins = int(numdata ** 0.5) #sqrt(n)
        #Scott's choice
        #spacing = 3.5 * numpy.var(data) ** 0.5 / numdata ** (1.0 / 3)
        #numbins = numpy.ceil((data_max - data_min) / spacing)

    print 'Histogram'
    print '*********'
    histogram, samplepoints = numpy.histogram(data, numbins)
    scale_factor = float(scale) / max(histogram)
    print 'N =', numdata, 'h =', numbins
    print '  x    %age'
    for bin_id in range(numbins):
        bin_ord = data_min + (bin_id + 0.5) / numbins * (data_max - data_min)
        percentage = float(histogram[bin_id]) / numdata * 100
        print '%7.3f' % bin_ord, '%7.4f%%' % percentage, \
            '*'*int(histogram[bin_id] * scale_factor)

    print
    print 'Gaussian kernel estimate'
    print '************************'
    kernel = gaussian_kde(data)
    kernel_sample = kernel.evaluate(samplepoints)
    for bin_id in range(numbins):
        bin_ord = samplepoints[bin_id]
        percentage = float(kernel_sample[bin_id]) / numdata * 100
        print '%7.3f' % bin_ord, '%7.4f%%' % percentage, \
            '*'*int(kernel_sample[bin_id] * numdata ** 0.5 * scale_factor)



def histogram_energies():
    filename = 'h2pc-data.h5'
    h5data = tables.openFile(filename, mode = 'r')
    excitation_energies = []
    for node in h5data.walkNodes():
        if node._v_name == 'Energy':
            if len(node) >= 2:
                energy = (node[1] - node[0]) / eV
                if 0 < energy < 10:
                    excitation_energies.append(energy)

    h5data.close()
    Histogram(excitation_energies)



if __name__ == '__main__':
    histogram_energies()
