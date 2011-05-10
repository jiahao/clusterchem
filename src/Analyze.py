#!/usr/bin/env python
"""
Analyze
~~~~~~~

Analyze data contained in a HDF5 database.

This is a standalone executable.

:copyright: 2011 Jiahao Chen
:license: Unknown

.. versionadded:: 0.1
"""

import numpy
from scipy.stats.kde import gaussian_kde
import sys
from Units import eV

try:
    import tables
except ImportError:
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()



def Histogram(data, numbins = None, scale = 72, barchar = '*'):
    """
    Prints a histogram to the console as well as a fit to a Gaussian kernel.

    :argument data: data to be histogrammed.
    :type data:
        A one-dimensional :py:mod:`numpy.ndarray`
        or iterable acceptable to :py:func:`numpy.asanyarray`

    :keyword numbins:
        Number of bins in the histogram to be created.
        (Default: None)
        If None, will be automatically estimated.

    :type numbins: integer or None

    :keyword integer scale: Longest bar in histogram to print. (Default: 72)

    :keyword string barchar:
        The character to print to console
        to render the histogram. (Default: '*')

    :returns: None

    This is a fancy wrapper around :py:func:`numpy.histogram`
    and :py:func:`scipy.stats.gaussian_kde`.

    .. versionadded:: 0.1
    """
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
            barchar * int(kernel_sample[bin_id] * numdata ** 0.5 * scale_factor)



def histogram_energies(h5filename = 'h2pc-data.h5'):
    """
    Queries a HDF5 file to extract energies and generates a histogram
    using :py:func:`Histogram`.

    :argument string h5filename:
        Name of a :term:`PyTables`-compatible HDF5 file containing energy data.

    :returns: None

    .. versionadded:: 0.1
    """
    h5data = tables.openFile(h5filename, mode = 'r')
    excitation_energies = []
    dump = []
    for node in h5data.walkNodes():
        if node._v_name == 'Energy':
            if len(node) >= 2:
                energy = (node[1] - node[0]) / eV
                if 0 < energy < 10:
                    excitation_energies.append(energy)
                    dump.append(str(energy) + '   ' + str(node))
    h5data.close()

    #Dump energies
    dump.sort()
    f = open('energies.dat', 'w')
    f.write('\n'.join(dump))
    f.close()

    Histogram(excitation_energies)



if __name__ == '__main__':
    histogram_energies()
