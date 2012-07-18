#!/usr/bin/env python
"""
Extracts CHARMM Cards from a HDF5 database.
"""

from glob import glob
import logging, os, shutil, tables, warnings

#Turn off NaturalNameWarning from PyTables
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

#Internal modules
from ParseOutput import ParseOutput
from Units import kcal_mol
import SGE
from HDF5Interface import CHARMM_CARD, LoadCHARMM_CARD, RecordIntoHDF5EArray, \
    RecordIntoHDF5CArray, OpenHDF5Table

from QChemTDDFTParser import QChemTDDFTParser
from QChemIO import QChemInput

DefaultCHARMMCARDSuffix='cor'

import os
import datetime

def ExtractCHARMMCARDs(h5filename):
    h5file = tables.openFile(h5filename, mode = 'r')
    for node in h5file.walkNodes('/', classname='Table'):
        if node.title == 'CHARMM CARD Coordinates':
            CHARMMCardFilename='.'.join((node._v_pathname[1:].\
                    replace('/CHARMM_CARD','').replace('/','--'), 
                    DefaultCHARMMCARDSuffix))
            print "Extracting CHARMM CARD coordinates to", CHARMMCardFilename
            try:
                WriteCHARMM_CARD(node, CHARMMCardFilename)
            except AssertionError:
                print "Skipping existing file"

def WriteCHARMM_CARD(table, CARDfilename):
    FileBuffer=["* AUTOMATICALLY EXTRACTED FROM HDF FILE"+CARDfilename,
                "* MACHINE: "+os.uname()[1],
                "* TIME: "+str(datetime.datetime.now()),
                str(table.nrows)]
    for row in table:
        mydata = row[:]
        data=dict()
        for field in ('AtomNo', 'ResidNo', 'Res', 'Type', 'Coord', 'SegID', 'ResID', 'Weighting'):
            data[field]=row[field]
        data['x']=data['Coord'][0]
        data['y']=data['Coord'][1]
        data['z']=data['Coord'][2]
        #This is the extended CRD format as specified in the CHARMM documentation, io.doc
        FileBuffer.append("%(AtomNo)10d%(ResidNo)10d  %(Res)8s  %(Type)8s\
%(x)20.10f%(y)20.10f%(z)20.10f  %(SegID)8s  %(ResID)8s%(Weighting)20.10f" % data)
    
    #assert not(os.path.exists(CARDfilename)), "File already exists: "+CARDfilename
    with open(CARDfilename, 'w') as f:
        f.write('\n'.join(FileBuffer))

if __name__ == '__main__':
    ExtractCHARMMCARDs('h2pc-tddft.h5')
