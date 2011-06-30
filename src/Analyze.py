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

import logging, numpy, sys
from scipy.stats.kde import gaussian_kde
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
        if node._v_name == 'Energy' and 'Fixed' in node._v_pathname:
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



def analyze_h2pc_tdip(h5filename = 'h2pc-data.h5'):
    """Calculate the dot product of the transition dipole
    with the H-H cordinate vector"""

    from numpy import dot, arccos, pi
    from numpy.linalg import norm


    buf = []
    h5data = tables.openFile(h5filename, mode = 'r')
    for node in h5data.walkNodes():
        if node._v_name == 'Dipole':
            #Extract transition dipole from matrix
            dipole = node[:,0,1]
            
            #Look up geometry
            path = node._v_pathname.split('/')
            resid  = path[-2]
            snapshot = path[1]

            try:
                geometry = h5data.getNode('/'+snapshot+'/CHARMM_CARD')
            except tables.exceptions.NoSuchNodeError: continue

            #Extract basic hydrogens
            atoms = [ x['Coord'] for x in geometry.iterrows() \
                      if x['ResID'] == resid and 'NQ' in x['Type'] ]

            if len(atoms) != 2: continue #silently fail
            assert len(atoms) == 2, """\
Wrong number of free base nitrogens in %r:%r
Expected 2 but found %d
Coordinates:
""" % (snapshot, resid, len(atoms)) + '\n'.join([str(x) for x in atoms])
            
            atomvec = atoms[1] - atoms[0]

            #Extract energies
            try:
                energies = h5data.getNode('/'.join(path[:-1]+['Energy']))
                excite = (energies[1] - energies[0]) /eV
            except (IndexError, tables.exceptions.NoSuchNodeError):
                excite = 0.0

            if norm(dipole) != 0.0:
                angle = arccos(dot(atomvec, dipole)/(norm(atomvec)*norm(dipole)))
                angle = min(angle, pi - angle)

                if len(energies) > 2: #Try next higher state
                    newexcite = (energies[2] - energies[0]) /eV
                    newdipole = node[:,0,2]
                    newangle = arccos(dot(atomvec, newdipole)/(norm(atomvec)*norm(newdipole)))
                    newangle = min(newangle, pi - newangle)
                    if newangle < angle:
                        #print 'Root flip at', snapshot, resid
                        angle = newangle
                        excite = newexcite

                if len(energies) > 2:
                    assert min(energies[2:]) > energies[1]

                #buf.append((snapshot, resid, excite, angle))
                if True or angle > 0.9 or excite > 2.1:
                    buf.append((snapshot, resid, excite, angle))
                    idx = len(buf)
                    #print idx, snapshot, resid, 'tddft-nonpol'
                    print snapshot, resid, '%10.6f eV' % excite, '%10.6f rad.' % angle

    #for idx, (snapshot, resid, excite, angle) in enumerate(buf):
    #    print snapshot, resid, '%10.6f eV' % excite, '%10.6f rad.' % angle
        #if angle < 0.9 or excite > 2.1:
        #    print idx+1, snapshot, resid, 'tddft-nonpol'

    h5data.close()


def LocateMissingData(h5filename = 'h2pc-data.h5', mollist = '/home/cjh/rmt/h2pc-data/sub-0/1/mol.list'):
    """
    Scans the HDF5 file for missing data.

    :param string h5filename: Name of HDF5 file
    :param string mollist: Name of file containing a list of active resids

    This scans the HDF5 file for missing or invalid energies and dipoles.

    Things that could go wrong:
    - The ground state energy is zero.
    - The excited state energy is zero or missing.
    - The transition dipole moment is zero.

    Outputs to console a list of serialized job information suitable for
    running with :py:mod:`SGEInterface`.
    """

    from numpy import zeros

    resids = []
    for l in open(mollist):
        try:
            resids.append(str(int(l)))
        except ValueError: pass

    #Names of data nodes and expected dimensions
    data = [('Energy', 2), ('Dipole', (3, 2, 2))]

    #Types of calculations
    calctypes = ['gs', 'dscf', 'td', 'gs-nonpol', 'dscf-nonpol', 'td-nonpol']

    #Enumerate coordinates
    
    h5data = tables.openFile(h5filename, mode = 'r')

    coords = []
    for node in h5data:
        if node._v_name == 'CHARMM_CARD':
            parentname = node._v_parent._v_name
            if '_' not in parentname:
                coords.append(parentname)

    isok = {}
    for coord in coords:
        for resid in resids:
            for calctype in calctypes:
                isok[(coord, resid, calctype)] = False

    for node in h5data:
        path = node._v_pathname.split('/')
        try:
            coord = path[1]
            env   = path[2] #Fixed or WithDrude
            resid = path[-2]
        except IndexError: continue
        
        for name, dim in data:
            if node._v_name == name:                
                if name == 'Energy':
                    okenergy1 = (node[0] != 0.0)
                    if env == 'Fixed': calctype = 'gs-nonpol'
                    else: calctype = 'gs'
                    isok[(coord, resid, calctype)] = okenergy1
                    if len(node) >= 2:
                        okenergy2 = (node[1] != 0.0) 
                        if env == 'Fixed': calctype = 'dscf-nonpol'
                        else: calctype = 'dscf'
                        isok[(coord, resid, calctype)] = okenergy2

                elif name == 'Dipole':                    
                    okdipole = (node.shape == dim) and not (node[:,:,:] == zeros(dim)).all()
                    if env == 'Fixed': calctype = 'td-nonpol'
                    else: calctype = 'td'
                    isok[(coord, resid, calctype)] = okdipole

    h5data.close()


    #Collates job data
    n = 1
    buf = []
    for coord in coords:
        for resid in resids:  
            #for calctype in [x for x in calctypes if 'nonpol' in x]:
            for calctype in ['gs-nonpol', 'dscf-nonpol', 'td-nonpol']:
                if not isok[(coord, resid, calctype)]:
                    buf.append('\t'.join((str(n), coord, resid, calctype)))
                    n += 1

    print '\n'.join(buf)
    

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Archives data into HDF5 database')
    parser.add_argument('--h5file', action = 'store', default = 'h2pc-data.h5', help = 'Name of HDF5 database file')
    parser.add_argument('--loglevel', action = 'store', default = logging.INFO, type = int, help = 'Logging level')
    args = parser.parse_args()

    logging.basicConfig(level = args.loglevel)

    #histogram_energies(args.h5file)
    analyze_h2pc_tdip(args.h5file)
    #LocateMissingData(args.h5file)
