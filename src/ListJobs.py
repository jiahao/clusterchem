#!/usr/bin/env python

"""
Lists jobs with missing data.

If run as a standalone script, executes :py:func:`ListJobsWithMissingData`.

.. versionadded:: 0.1
"""

import tables

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Lists missing jobs in \
HDF5 file')
    parser.add_argument('-i', '--h5file', action = 'store',
			default = 'h2pc-data.h5',
			help = 'Name of HDF5 database file')

    parser.add_argument('-m', '--mollist', action = 'store',
			default = '../mol.list',
			help = 'Name of residue list file')
    args = parser.parse_args()


def ListJobsWithMissingData(filename = 'h2pc-data.h5',
    mol_list = '../mol.list'):
    """
    Generates a list of molecule that lack nonpolarizable
    :term:`QM/MM` transition dipole data.

    Prints to console: index, snapshot_id, site_id, and type of
    data missing (currently, only td-nonpol)

    :param string filename: Name of HDF5 data file.
    :param string mol_list: File containing list of relevant site ids.

    :returns: :rtype: None
    """

    resids = set()
    for l in open(mol_list):
        try:
            resid = int(l)
            resids.add(str(resid))
        except ValueError:
            pass


    badnodes = set()
    h5data = tables.openFile(filename, mode = 'r')
    for topnode in h5data.iterNodes('/'):
        if 'Model' in topnode._v_name:
            continue

        testnameroot = '/'+topnode._v_name + '/Fixed/Sites/'
        for resid in resids:
            nodename = testnameroot + resid
            try:
                node = h5data.getNode(nodename, 'Energy')
                if node[0] == 0.0:
                    badnodes.add(nodename)
                    print node._v_pathname, 'Ground state energy is 0'
                if len(node) < 2 or node[1] == 0.0:
                    badnodes.add(nodename)
                    print node._v_pathname, 'Excited state energy is 0 or missing'
            except tables.NoSuchNodeError:
                badnodes.add(nodename)
                print nodename, 'No energy data'

            try:
                node = h5data.getNode(nodename, 'Dipole')
            except tables.NoSuchNodeError:
                badnodes.add(nodename)
                print nodename, 'No dipole data'
                
    for n, node in enumerate(sorted(badnodes)):
        t = node.split('/')
        print n+1, t[1], t[4], 'tddft-nonpol'



if __name__ == '__main__':
    ListJobsWithMissingData(args.h5file, args.mollist)

