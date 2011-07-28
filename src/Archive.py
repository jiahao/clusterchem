#!/usr/bin/env python
"""
Archives data into a :term:`HDF5` database.

This is a standalone executable requiring :term:`PyTables`
and the homebrew :py:mod:`QChemIO`.

.. versionadded:: 0.1

"""

#TODO SHOULD WE ARCHIVE ELECTRONIC STRUCTURE CONVERGENCE INFORMATION?

from glob import glob
import logging, numpy, os, os.path, shutil, sys, warnings

try:
    import tables
except ImportError:
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()

#Turn off NaturalNameWarning from PyTables
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

#My custom module
try:
    import QChemIO
except ImportError:
    sys.stderr.write("""Could not find QChemIO module.
Please check that PYTHONPATH is set correctly.
""")
    exit()

#Internal modules
from ParseOutput import ParseOutput
from Units import kcal_mol
import SGE



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = 'Archives data into HDF5 database')
    parser.add_argument('workingdir', nargs = '?', default = '.', help = 'Root to recurse from')
    parser.add_argument('--h5file', action = 'store', default = 'h2pc-data.h5', help = 'Name of HDF5 database file')
    parser.add_argument('--taskfile', action = 'store', default = 'sgearraytasklist.txt', help = 'Name of SGE array task file')
    parser.add_argument('--nuke', action = 'store_true', default = False, help = 'Deletes files after parsing')
    parser.add_argument('--nukeinc', action = 'store_true', default = False, help = 'Deletes incomplete output files')
    parser.add_argument('--checkpoint', action = 'store_true', default = False, help = 'Create new Pytables checkpoint in HDF5 file')
    parser.add_argument('--force', action = 'store_true', default = False, help = 'Forces reloading of existing data')
    parser.add_argument('--loglevel', action = 'store', default = logging.INFO, type = int, help = 'Logging level')
    args = parser.parse_args()



class CHARMM_CARD(tables.IsDescription):
    """
    :term:`PyTables` data structure for storing a :term:`CHARMM card file` in a :term:`HDF5` database.

    .. versionadded:: 0.1
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
Expected %d but only %d atoms were received.""" % (numatoms, thisnumatoms))

    logger.info('Loaded CHARMM file '+ CARDfile + ' into HDF5 table ' + str(h5table))
    h5table.flush()



def CheckIfRunning(filename = None, taskfile = args.taskfile):
    """
    Is the process responsible for producing filename still running?

    Need to infer this from qstat

    Get list of running SGE processes
    then check working directories
    """

    logger = logging.getLogger('Archive.CheckIfRunning')

    #Break down filename into job info
    filename_token = filename.split(os.sep)
    myjobtype = filename_token[-2]
    mymol = filename_token[-3]
    mysnapshot = filename_token[-4]

    #Retrieve mapping between job array id and job info
    tasklist = []
    for line in open('sgearraytasklist.txt'):
        tasklist.append(line.split())
    myjobids = [job[0] for job in tasklist if job[1] == myjobtype and job[2] == mymol and job[3] == mysnapshot]

    if len(myjobids) == 0: 
        logger.warning('Cannot find job array corresponding to %s, assuming not running anymore', filename)
        return False

    assert len(myjobids) == 1, 'More than one job array matching '+filename
    
    myjobid = myjobids[0]

    #Poll queue
    for job in SGE.SGE().getuserjobs():
        if myjobid in job._ja_tasklist: #It's still there!
            logger.info('Job is in SGE queue with status %s', job.state)
            return True

    #No matches
    logger.info('Job is not listed in SGE queue')
    return False



def RunQChemAgain(filename):
    """
    Rewrite Q-Chem input deck with different convergence algorithm.
    DIIS --> 
    If funny behavior initially, RCA_DIIS
    If oscillations toward the end, DIIS_GDM
    or MOM - set MOM_START
    If MOM is to be used to aid
    convergence, an SCF without MOM should be run to determine when the SCF
    starts oscillating. MOM should be set to start just before the oscillation
    MOM_START
    If MOM is to be used to aid
    convergence, an SCF without MOM should be run to determine when the SCF
    starts oscillating. MOM should be set to start just before the oscillation
    """

    #TODO
    logger = logging.getLogger('Archive.RunQChemAgain')
    logger.critical('Run again manually: %s', filename)
    return False



def RecordIntoHDF5CArray(data, h5data, location, name, atom = tables.atom.FloatAtom(), DoOverwrite = False):
    """Data to record
       HDF5 file
       Location in HDF5
       Name to make"""

    logger = logging.getLogger('Archive.RecordIntoHDF5CArray')

    try:
        myarray = h5data.createCArray(where = location, name = name, atom = atom,
           shape = data.shape, createparents = True)
    except tables.exceptions.NodeError:
        #Already in there
        if DoOverwrite:
            logger.info('Overwriting existing HDF5 CArray '+os.path.join(location, name))
            myarray= h5data.getNode(location, name)
            #Same shape? If not, delete and recreate
            if myarray.shape != data.shape:
                h5data.removeNode(where = location, name = name)
                myarray = h5data.createCArray(where = location, title = name,
                    name = name, atom = atom, shape = data.shape, createparents = True)
        else:
            logger.info('HDF5 CArray already exists: '+os.path.join(location, name)+'\nSkipping.')
            return

    myarray[:] = data
 
    logger.info('Archiving HDF5 CArray ' + os.path.join(location, name) + ':\n' + str(myarray[:]))



def RecordIntoHDF5EArray(data, h5data, location, name, position = None, atom = tables.atom.FloatAtom(), DoOverwrite = False):

    logger = logging.getLogger('Archive.RecordIntoHDF5EArray')

    #Currently hard coded for one-dimensional arrays
    try:
        myarray = h5data.createEArray(where = location, name = name,
                      atom = atom, shape = (0,), title = name,
                      expectedrows = 2, createparents = True)
    except tables.exceptions.NodeError:
        #Already in there
        if DoOverwrite:
            myarray = h5data.getNode(location, name)
        else:
            logger.info('HDF5 EArray already exists: '+os.path.join(location, name)+'\nSkipping.')
            return

    #Data is either a scalar or an entire array
    #If scalar, position MUST be specified

    try:
        N = len(data)
        while N > len(myarray):
            myarray.append([numpy.NaN])

        myarray[:] = data
        logger.info('Archiving HDF5 EArray ' + os.path.join(location, name) + ':')
        for idx, datum in enumerate(data):
            logger.info('Archived '+str(datum)+' in position '+str(idx))
    except TypeError:
        if position == None:
            logger.error("""Requested archival of scalar datum into EArray but position in EArray was not given.
Datum was NOT archived.""")
            return

        while position >= len(myarray):
            myarray.append([numpy.NaN])
        
        myarray[position] = data

        logger.info('Archiving HDF5 EArray ' + os.path.join(location, name) + ':')
        logger.info('Archived '+str(data)+' in position '+str(position))



def LoadEmUp(path = '.', h5filename = 'h2pc-data.h5', Nuke = False, NukeInc = False, DoCheckPoint = False, DoOverwrite = False):
    """
    Recursively reads data from all files in a given path and loads it into a
    :term:`HDF5` database.
    
    :argument string path: Path to read data from recursively. (Default: '.')
    
    :argument string h5filename: Name of :term:`HDF5` database to create. \
        or update (Default: 'h2pc-data.h5')

    :argument Boolean Nuke: Delete files after successful archival (Default: False)
    :argument Boolean NukeInc: Delete incomplete output files detected (Default: False)
    :argument Boolean DoCheckPoint: Create new check point after successful
         archival (Default: False)

    :argument Boolean DoOverwrite: Force overwrite of existing data with data in existing files (Default: False)

    :returns: None

    .. TODO :: replace with a real way to test if node exists in Pytables
    .. TODO :: Failed to detect symlink

    .. versionadded:: 0.1
    """
    logger = logging.getLogger('Archive.LoadEmUp')

    compress = tables.Filters(complevel = 9, complib = 'zlib')
    h5data = tables.openFile(h5filename, mode = 'a', filters = compress,
        title = 'H2Pc')
    if not h5data.isUndoEnabled(): h5data.enableUndo(filters = compress)
    for root, subdirs, files in os.walk(path):
        for cwdfile in files:

            filename = os.path.join(root, cwdfile) #Relative path-qualified filename



            #Skip symlinks
            if os.path.islink(filename): continue

            ##############################
            # Parse CHARMM coordinate file
            ##############################

            if len(cwdfile) > 4 \
                and cwdfile[-4:] == '.cor' \
                and cwdfile[-6:-4] != '.d':

                # CARDs with Drude particles are skipped since they
                # will need to be re-equilibrated as part of the SCF process

                cardtitle = 'CHARMM CARD Coordinates'
                cardname = name = 'CHARMM_CARD'
                
                cwdh5name = cwdfile[:-4]#.replace('-', '_')

                location = os.path.join('/' + root, cwdh5name)
                try:
                    thischarmmcard = h5data.createTable(
                        name = cardname, \
                        where = location, \
                        description = CHARMM_CARD, \
                        title = cardtitle, createparents = True)
                    LoadCHARMM_CARD(thischarmmcard, cwdfile)
                except tables.exceptions.NodeError:
                    if DoOverwrite:
                        logger.info('Overwriting existing CHARMM CARD in HDF5 table '+os.path.join(location, name))
                        h5data.removeNode(location, name)
                        thischarmmcard = h5data.createTable(
                            where = location, \
                            name = cardname, \
                            description = CHARMM_CARD, \
                            title = cardtitle, createparents = True)
                        LoadCHARMM_CARD(thischarmmcard, cwdfile)
                    else:
                        logger.info('HDF5 table already exists: '+os.path.join(location, name)+ ' skipping archival.')
                        continue


            ####################################
            # Parse CHARMM or Q-Chem output file
            ####################################

            if len(cwdfile) > 4 and cwdfile[-4:] == '.out':

                #Assume path has the structure [.../]GeomName/Site/jobtype
                Site = root.split(os.sep)[-2]
                GeometryName = root.split(os.sep)[-3]
                calctype = os.path.split(root)[1]

                State = Energy = Energies = Dipole = Environment = None
                data = ParseOutput(filename)
                isDone = data[-1]
                EnergyList = data[1]
                Dipole = data[2]

                if not isDone:
                    #TODO Check if job is still running on SGE
                    isRunning = CheckIfRunning(filename)
                    if isRunning:
                        #If it is running, do nothing
                        logger.info('Skipping incomplete calculation in file '+ filename)
                        continue
                    
                    #TODO If it isn't, check for errors
                    
                    #If there are convergence errors:
                    RunQChemAgain(filename)

                    if NukeInc:
                        logger.info('Deleting incomplete calculation in directory '+ root)
                        shutil.rmtree(root)
                        continue

                if '-nonpol' in calctype:
                    Environment = 'Fixed'
                else:
                    Environment = 'WithDrude'
                    #TODO Record Drude particle positions
                    logger.warning("Here I should be recording the Drude particle positions but I'm not!")
 
                if 'gs' in calctype:
                    State = 0
                    if len(EnergyList) > 0:
                        Energy = EnergyList[-1]
                    else: continue
                
                elif 'dscf' in calctype:
                    State = 1
                    if len(EnergyList) > 0:
                        Energy = EnergyList[-1]

                elif cwdfile == 'td.out' and 'td' in calctype:
                    Dipole = data[2]

                elif 'tddft' in calctype:
                    from QChemTDDFTParser import QChemTDDFTParser
                    QCtddft = QChemTDDFTParser(filename)
                    TDAEnergies, TDADipole, Energies, Dipole = QCtddft.GetEnergiesAndDipole()
                    TDAMultiplicities, TDAOscillatorStrengths, Multiplicities, OscillatorStrengths = QCtddft.GetOtherData()

                    if TDAEnergies != None: TDAEnergies /= kcal_mol
                    if Energies != None: Energies /= kcal_mol
               
                elif cwdfile == 'qchem.out':
                    """Assume all energy calculations go through the
                    CHARMM/Q-Chem interface.

                    This avoids double processing of the polarizable simulations"""
                    continue


                ################################################
                # Record completed calculations into HDF5 file #
                ################################################

                #Update only if job is done
                if not isDone: continue

                location = os.path.join('/' + GeometryName, Environment, 'Sites', Site)

                #Update a single state's energy
                if Energy != None:
                    RecordIntoHDF5EArray(Energy * kcal_mol, h5data, location, 'Energy', State,  DoOverwrite = DoOverwrite)
                #Update all energies
                if Energies != None:
                    RecordIntoHDF5EArray(Energies * kcal_mol, h5data, location, 'Energy', DoOverwrite = DoOverwrite)

                #Update a transition dipole
                if Dipole != None:
                    RecordIntoHDF5CArray(Dipole, h5data, location, 'Dipole', DoOverwrite = DoOverwrite)

                if 'tddft' in calctype: #Store also TDA energies and dipoles
                    RecordIntoHDF5EArray(TDAEnergies * kcal_mol, h5data, location, 'TDAEnergy', DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(TDADipole, h5data, location, 'TDADipole', DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(TDAMultiplicities, h5data, location, 'TDAMultiplicities', atom = tables.atom.IntAtom(), DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(Multiplicities, h5data, location, 'Multiplicities', atom = tables.atom.IntAtom(), DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(TDAOscillatorStrengths, h5data, location, 'TDAOscillatorStrengths', DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(OscillatorStrengths, h5data, location, 'OscillatorStrengths', DoOverwrite = DoOverwrite)


                ############
                # Clean up #
                ############

                if Nuke:
                    logger.info('Deleting files in directory ' + root)
                    shutil.rmtree(root)

        #Clean up dangling subdirectories
        #XXX HORRIBLE HACK
        #This needs to check if the job is still running first...
        if False and Nuke and len(subdirs) == 0 and 'h2zn' in root and \
            (len(root.split(os.path.sep)) - len(path.split(os.path.sep))) <= 2:
            #We only want to go up to path/h2zn-???{,?}/??? to check for lack
            #of further subdirectories
            logger.info('Cleaning up directory ' + root)
            shutil.rmtree(root)

            #Naive recursion all the way up to original path, but
            #DO NOT DELETE ORIGINAL
            rootup = root
            while (len(rootup.split(os.path.sep)) - len(path.split(os.path.sep))) > 1:
                #Go up - there MUST be a better way to do this!
                rootup = os.path.sep.join(rootup.split(os.path.sep)[:-1])
                subdirs = os.walk(rootup).next()[1]
                if len(subdirs) == 0:
                    logger.info('Cleaning up directory ' + rootup)
                    shutil.rmtree(rootup)
        
        #Clean up redirected SGE stdout/stderr streams
        if root == path:
            filenames = glob('*.stdout.log')
            if len(filenames) > 0:
                logger.info('Cleaning up SGE stdout/stderr streams in ' + root)

            for filename in filenames:
                os.unlink(filename)


    if h5data.isUndoEnabled() and DoCheckPoint:
        mark = h5data.mark()
        logger.info('Checkpointing HDF5 database at mark point '+str(mark))

    h5data.close()



if __name__ == '__main__':
    logging.getLogger('CHARMMUtil')

    logging.basicConfig(level = args.loglevel)

    LoadEmUp(args.workingdir, args.h5file, args.nuke, args.nukeinc, args.checkpoint, args.force)
