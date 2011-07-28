#!/usr/bin/env python

"""
Lists jobs with missing data.

If run as a standalone script, executes :py:func:`ListJobsWithMissingData`.

.. versionadded:: 0.1
"""

import tables

def ListJobsWithMissingData(filename = 'h2pc-data.h5',
    mol_list = '/home/cjh/rmt/inputdata/mol.list'):
    """
    Generates a list of molecule that lack nonpolarizable
    :term:`QM/MM` transition dipole data.

    Prints to console: index, snapshot_id, site_id, and type of
    data missing (currently, only td-nonpol)

    :param string filename: Name of HDF5 data file.
    :param string mol_list: File containing list of relevant site ids.

    :returns: :rtype: None
    """

    badnodes = []

    resids = []
    for l in open(mol_list):
        resids.append(l.strip())

    h5data = tables.openFile(filename, mode = 'r')

    for topnode in h5data.iterNodes('/'):
        testnameroot = '/'+topnode._v_name + '/Fixed/Sites/'
        for resid in resids:
            nodename = testnameroot + resid
            """
            try:
                node = h5data.getNode(nodename, 'Energy')
                if node[0] == 0.0:
                    badnodes.append((node._v_pathname, 0))
                if len(node) < 2 or node[1] == 0.0:
                    badnodes.append((node._v_pathname, 1))
            except tables.NoSuchNodeError:
                badnodes.append((nodename+'/Energy', 0))
                badnodes.append((nodename+'/Energy', 1))
            """
            try:
                node = h5data.getNode(nodename, 'Dipole')
            except tables.NoSuchNodeError:
                badnodes.append((nodename+'/Dipole', 0))
                
    for n, (node, _) in enumerate(badnodes):
        t = node.split('/')
        print n+1, t[1], t[4], 'td-nonpol'



if __name__ == '__main__':
    ListJobsWithMissingData()

