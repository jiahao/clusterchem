#!/usr/bin/env python
"""
Extracts the list of molecules to iterate over from an existing HDF5 database
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
    h5data = tables.openFile(h5filename, mode = 'r')
    try:
        MolList = h5data.getNode('/Model', 'MolList', 'Table')
        return map(int, MolList)

    except tables.exceptions.NoSuchNodeError:
        #Plan B: regenerate it based on the nodes traversed
        ResIDs = dict()
        for node in h5data.walkNodes('/', classname='Group'):
            try:
                ThisResID = int(node._v_name)
                ResIDs[ThisResID]=True
            except ValueError:
                continue
        return sorted(ResIDs.keys())

if __name__ == '__main__':
    data = ExtractCHARMMCARDs('h2pc-tddft.h5')
    for resid in data:
        print resid

