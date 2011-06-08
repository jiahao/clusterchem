#!/usr/bin/env python

"""
Miscellaneous utilities

If run as a standalone script, will execute :func:`dedrude`, passing
command-line arguments to it.

.. versionadded:: 0.1
"""

import logging

def dedrude(infilename, outfilename):
    """
    Strips Drude particles from :term:`CHARMM card file`

    :param string infilename:
        Name of :term:`CHARMM card file` to read in that
        describes a geometry that contains Drude particles.

    :param string outfilename:
        Name of :term:`CHARMM card file` to write out
        that contains the geometry without Drude particles.

    This function removes Drude particles from a CHARMM CARD file.
    Drude particles have an element name that begins with D.

    .. versionadded:: 0.1
    """
    buf = []    
    state = 'header'

    log = logging.getLogger()
    log.info('Stripping Drude particles from CHARMM CARD '+infilename)

    for l in open(infilename):
        t = l.split()
        if l[0] == '*':
            buf.append(l)
        elif state == 'header':
            if t[1] == 'EXT':
                numatomslinenum = len(buf) - 1
                buf.append(l)
                atom = 0
                state = 'atom'
        elif state == 'atom':
            if t[3][0] != 'D': #If not Drude
                atom += 1
                buf.append('%10d' % atom)
                buf.append(l[10:])
    buf[numatomslinenum] = '%10d' % atom + '  EXT'

    f = open(outfilename, 'w')
    f.write(''.join(buf))
    f.close()
    log.info('Wrote new CHARMM CARD '+outfilename)


if __name__ == '__main__':
    import sys

    logging.getLogger('CHARMMUtil')
    if len(sys.argv) < 3:
        logging.error('Specify input and output CHARMM CARD filenames')
        exit()

    logging.basicConfig(level = logging.INFO)    
    dedrude(sys.argv[1], sys.argv[2])
