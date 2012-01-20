#!/usr/bin/env python

"""
Emergency recovery utility

This program is designed to recreate data from the log file
created by the Archive utility. Based on the textual representations,
it will try to recreate a HDF5 containing the same data.
"""

import bz2, logging, numpy, os, tables
from numpy import array

from Units import kcal_mol
from HDF5Interface import RecordIntoHDF5CArray, RecordIntoHDF5EArray


def RecoverFromLog(h5filename, bzlogfilename):
    mode = 'seek'
    doflush = False
    buf = []
    stacklevel = 0

    compress = tables.Filters(complevel = 9, complib = 'zlib')
    h5data = tables.openFile(h5filename, mode = 'a', filters = compress,
                           title = 'Recovered')
    
    for line in bz2.BZ2File(bzlogfilename):
        data = None
        if 'Archiving HDF5 CArray' in line:
                location = line.split()[-1][:-1]
                mode = 'carray'
                doflush = True
        elif 'Archiving HDF5 EArray' in line:
                location = line.split()[-1][:-1]
                mode = 'earray'
                doflush = True
        if mode == 'earray':
            if 'position' in line:
                t = line.split()
                position = int(t[-1])
                data = float(t[-4])
        elif mode == 'carray':
            linebuf = line.split()
            if len(linebuf)>0 and linebuf[0].count('[') == len(linebuf[0]):
                token = linebuf.pop(0)
                linebuf[0] = token+linebuf[0]
            if len(linebuf)>0:
                buf.append(','.join(linebuf).replace('NaN','numpy.nan'))
            stacklevel += line.count('[')-line.count(']')
            if stacklevel == 0:
                parseme = ','.join(buf)
                if parseme != '':
                    try:
                        data = numpy.array(eval(parseme))
                    except (NameError, SyntaxError):
                        pass
                    buf = []
              
        if data is not None:
            basepath, name = os.path.split(location)
            if mode == 'earray':
                print location, mode, position, repr(data)
                RecordIntoHDF5EArray(data * kcal_mol, h5data, basepath, name, position)
            else:
                print location, mode, repr(data)
                RecordIntoHDF5CArray(data, h5data, basepath, name)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Recovers data from archival log')
    parser.add_argument('--h5file', action = 'store', default = 'recovered.h5',
        help = 'Name of HDF5 database file to write')
    parser.add_argument('--logfile', action = 'store', default = 'archival.log.bz2',
        help = 'Name of archival log to read')
    args = parser.parse_args()

    RecoverFromLog(args.h5file, args.logfile)
    

