#!/usr/bin/env python
"""
Extracts CHARMM Cards from a HDF5 database.
"""

import datetime, os, tables, warnings

#Turn off NaturalNameWarning from PyTables
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

DefaultCHARMMCARDSuffix = 'cor'

def ExtractCHARMMCARDs(h5filename):
    """
    Extract all CHARMM CARD entries in HDF5 file and write them out
    as separate files on disk.
    """
    #pylint: disable=W0212
    
    h5file = tables.openFile(h5filename, mode = 'r')
    for node in h5file.walkNodes('/', classname='Table'):
        if 'CHARMM CARD coordinates'.upper() not in node.title.upper():
            continue

        CHARMMCardFilename = '.'.join((node._v_pathname[1:].\
                replace('/CHARMM_CARD','').replace('/','--'), 
                DefaultCHARMMCARDSuffix))
        print "Extracting CHARMM CARD coordinates to", CHARMMCardFilename
        try:
            WriteCHARMM_CARD(node, CHARMMCardFilename)
        except AssertionError:
            print "Skipping existing file"

def WriteCHARMM_CARD(table, CARDfilename, isExtended=False, Overwrite=False):
    """
    Write HDF5 CHARMM CARD table to a text file on disk.
    """

    FileBuffer = ["* AUTOMATICALLY EXTRACTED FROM HDF FILE" + CARDfilename,
                "* MACHINE: " + os.uname()[1],
                "* TIME: " + str(datetime.datetime.now()),
                str(table.nrows)]
    for row in table:
        data = dict()
        for field in ('AtomNo', 'ResidNo', 'Res', 'Type', 'Coord', 'SegID',
                'ResID', 'Weighting'):
            data[field] = row[field]
        data['x'] = data['Coord'][0]
        data['y'] = data['Coord'][1]
        data['z'] = data['Coord'][2]
        if isExtended:
            #This is the extended CRD format as specified in
            #the CHARMM documentation, io.doc
            FileBuffer.append("%(AtomNo)10d%(ResidNo)10d  %(Res)8s  %(Type)8s\
%(x)20.10f%(y)20.10f%(z)20.10f  %(SegID)8s  %(ResID)8s%(Weighting)20.10f" \
% data)
        else:
            #This is the normal CRD format
            FileBuffer.append("%(AtomNo)5d%(ResidNo)5d %(Res)4s %(Type)4s\
%(x)10.5f%(y)10.5f%(z)10.5f %(SegID)4s%(ResID)4s%(Weighting)10.5f" % data)

    if not Overwrite:
        assert not(os.path.exists(CARDfilename)), "File already exists: " \
            + CARDfilename

    with open(CARDfilename, 'w') as f:
        f.write('\n'.join(FileBuffer))

if __name__ == '__main__':
    ExtractCHARMMCARDs('h2pc-data.h5')
