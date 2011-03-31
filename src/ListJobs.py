#!/usr/bin/env python

import tables

def main():
    badnodes = []

    resids = []
    for l in open('/home/cjh/rmt/inputdata/mol.list'):
        resids.append(l.strip())

    filename = 'h2pc-data.h5'
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
                
    for n, (node, idx) in enumerate(badnodes):
        t = node.split('/')
        print n+1, t[1], t[4], 'td-nonpol'


if __name__ == '__main__':
    main()
    #cleanup()
