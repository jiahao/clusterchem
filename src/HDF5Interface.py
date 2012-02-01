#!/usr/bin/env python

"""
HDF5 Interface classes and functions
"""

import logging, numpy, os, tables

class CHARMM_RTF_Short(tables.IsDescription):
    """
    PyTables data structure containing some CHARMM RTF information
    pertaining to atomic charges
    """
    #: Atom Type
    Type = tables.StringCol(4)
    #: Atomic charge
    Charge = tables.Float64Col()

    def __init__(self):
        tables.IsDescription.__init__(self)


class ResidList(tables.IsDescription):
    """
    PyTables data structure containing a list of molecule residue lists
    """
    #: Residue ID
    ResID = tables.UInt32Col()

    def __init__(self):
        tables.IsDescription.__init__(self)


class CHARMM_CARD(tables.IsDescription):
    """
    PyTables data structure for storing a CHARMM card file in a HDF5 database.
    """
    #: Atomic number
    AtomNo = tables.UInt64Col()
    #: Residue number
    ResidNo = tables.UInt64Col()
    #: Residue label of up to 4 characters
    Res = tables.StringCol(4)
    #: Type
    Type = tables.StringCol(4)
    #: Coordinates
    Coord = tables.Float32Col(shape = 3,)
    #: SegID
    SegID = tables.StringCol(4)
    #: ResID
    ResID = tables.StringCol(4)
    #: Weighting
    Weighting = tables.Float32Col()

    def __init__(self):
        tables.IsDescription.__init__(self)


def LoadCHARMM_CARD(h5table, CARDfile):
    """
    Reads a :term:`CHARMM card file` and saves it into a :term:`HDF5` database.
    
    :argument h5table: Name of :term:`HDF5` table
    :type h5table: :term:`PyTables` :py:class:`tables.File` object
    
    :argument string CARDfile: Name of :term:`CHARMM card file`

    :returns: None

    .. versionadded:: 0.1
    """

    logger = logging.getLogger('Archive.LoadCHARMM_CARD')

    state = 'title'
    for l in open(CARDfile):
        if state == 'title':
            if not (len(l) > 0 and l[0] == '*'):
                state = 'numatoms'
        if state == 'numatoms':
            numatoms = int(l.split()[0])
            thisnumatoms = 0
            state = 'read'
        elif state == 'read':
            t = l.split()
            data = h5table.row
            data['AtomNo'] = int(t[0])
            data['ResidNo'] = int(t[1])
            data['Res'] = t[2]
            data['Type'] = t[3]
            data['Coord'] = float(t[4]), float(t[5]), float(t[6])
            data['SegID'] = t[7]
            data['ResID'] = t[8]
            data['Weighting'] = float(t[9])
            data.append()
            thisnumatoms += 1

    if thisnumatoms != numatoms:
        logger.error("""Did not read in expected number of atoms.
Expected %d but only %d atoms were received.""", numatoms, thisnumatoms)

    logger.info('Loaded CHARMM file %s into HDF5 table %s', CARDfile,
                str(h5table))
    h5table.flush()


def OpenHDF5Table(h5data, where, name, description, title,
    DoOverwrite = False, DoAppend = True):
    """
    @param h5data HDF5 file handler
    @param where Parent group
    @param name Name of new table
    @param description Object describing table
    @param title String description
    @param DoOverwrite
    @param DoAppend
    """

    logger = logging.getLogger('OpenHDF5Table')
    try:
        h5table = h5data.getNode(where, name)
        if DoOverwrite:
           logger.info('Overwriting existing %s in HDF5: %s',
                      name, os.path.join(where, name))
           h5data.removeNode(where, name)
           
           h5table = h5data.createTable(where, name, description, title,
               createparents = True)

    except tables.exceptions.NodeError: #Table does not exist
        h5table = h5data.createTable(where, name, description, title,
            createparents = True)

    return h5table


def RecordIntoHDF5CArray(data, h5data, location, name,
    atom = tables.atom.FloatAtom(), DoOverwrite = False):
    """
    Puts data into HDF5 CArray.

    If destination does not exist, will be automatically created.

    @param data Data to record
    @param h5data HDF5 file
    @param location Location in HDF5
    @param name Name to make
    @param atom Pytables data type to record
    @param DoOverwrite whether to overwrite existing data
    """

    logger = logging.getLogger('Archive.RecordIntoHDF5CArray')

    try:
        myarray = h5data.createCArray(where = location, name = name,
           atom = atom, shape = data.shape, createparents = True)
    except tables.exceptions.NodeError:
        #Already in there
        if DoOverwrite:
            logger.info('Overwriting existing HDF5 CArray: %s', 
                        os.path.join(location, name))
            myarray = h5data.getNode(location, name)
            #Same shape? If not, delete and recreate
            if myarray.shape != data.shape:
                h5data.removeNode(where = location, name = name)
                myarray = h5data.createCArray(where = location, title = name,
                    name = name, atom = atom, shape = data.shape,
                    createparents = True)
        else:
            logger.info('HDF5 CArray already exists: %s\nSkipping,', 
                        os.path.join(location, name))
            return

    myarray[:] = data
 
    logger.info('Archiving HDF5 CArray %s:\n%s', os.path.join(location, name),
                str(myarray[:]))


def RecordIntoHDF5EArray(data, h5data, location, name, position = None,
    atom = tables.atom.FloatAtom(), DoOverwrite = False):

    """
    Puts data into HDF5 CArray.

    If destination does not exist, will be automatically created.

    @param data Data to record
    Data is either a scalar or an entire array.
    If scalar, position MUST be specified.
    @param h5data HDF5 file
    @param location Location in HDF5
    @param name Name to make
    @param position EArray index to insert at. (Optional)
    @param atom Pytables data type to record
    @param DoOverwrite whether to overwrite existing data

    @note Currently hard coded for one-dimensional arrays
    """

    logger = logging.getLogger('Archive.RecordIntoHDF5EArray')

    try:
        myarray = h5data.createEArray(where = location, name = name,
                      atom = atom, shape = (0,), title = name,
                      expectedrows = 2, createparents = True)
    except tables.exceptions.NodeError:
        #Already in there
        if DoOverwrite:
            myarray = h5data.getNode(location, name)
        else:
            logger.info('HDF5 EArray already exists: %s\nSkipping.',
                        os.path.join(location, name))
            return

    try:
        N = len(data)
        while N > len(myarray):
            myarray.append([numpy.NaN])

        myarray[:] = data
        logger.info('Archiving HDF5 EArray %s:', os.path.join(location, name))
        for idx, datum in enumerate(data):
            logger.info('Archived '+str(datum)+' in position '+str(idx))
    except TypeError:
        if position == None:
            logger.error("""Requested archival of scalar datum into EArray but \
position in EArray was not given.
Datum was NOT archived.""")
            return

        while position >= len(myarray):
            myarray.append([numpy.NaN])
        
        myarray[position] = data

        logger.info('Archiving HDF5 EArray %s:', os.path.join(location, name))
        logger.info('Archived %s in position %s', str(data), str(position))
