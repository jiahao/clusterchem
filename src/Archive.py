#!/usr/bin/env python
"""
Archives data into a :term:`HDF5` database.

This is a standalone executable requiring :term:`PyTables`
and the homebrew :py:mod:`QChemIO`.

.. versionadded:: 0.1

"""

#TODO SHOULD WE ARCHIVE ELECTRONIC STRUCTURE CONVERGENCE INFORMATION?

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

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description = 'Archive into HDF5 file')
    parser.add_argument('workingdir', nargs = '?', default = '.',
                        help = 'Root to recurse from')
    parser.add_argument('--h5file', action = 'store', default = 'h2pc-data.h5',
                        help = 'Name of HDF5 database file')
    parser.add_argument('--taskfile', action = 'store',
                        default = 'sgearraytasklist.txt',
                        help = 'Name of SGE array task file')
    parser.add_argument('--nuke', action = 'store_true', default = False,
                        help = 'Deletes files after parsing')
    parser.add_argument('--nukeinc', action = 'store_true', default = False,
                        help = 'Deletes incomplete output files')
    parser.add_argument('--checkpoint', action = 'store_true', default = False,
                        help = 'Create new Pytables checkpoint in HDF5 file')
    parser.add_argument('--force', action = 'store_true', default = False,
                        help = 'Forces reloading of existing data')
    parser.add_argument('--loglevel', action = 'store', default = logging.INFO,
                        type = int, help = 'Logging level')
    args = parser.parse_args()



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
    for line in open(taskfile):
        tasklist.append(line.split())
    myjobids = [job[0] for job in tasklist if job[1] == myjobtype and 
                job[2] == mymol and job[3] == mysnapshot]

    if len(myjobids) == 0: 
        logger.warning("""Cannot find job array corresponding to %s
Assuming not running anymore""", filename)
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



def LoadEmUp(path = '.', h5filename = 'h2pc-data.h5', Nuke = False,
             NukeInc = False, DoCheckPoint = False, DoOverwrite = False):
    """
    Recursively reads data from all files in a given path and loads it into a
    :term:`HDF5` database.
    
    :argument string path: Path to read data from recursively. (Default: '.')
    
    :argument string h5filename: Name of :term:`HDF5` database to create. \
        or update (Default: 'h2pc-data.h5')

    :argument Boolean Nuke: Delete files after successful archival
         (Default: False)
    :argument Boolean NukeInc: Delete incomplete output files detected
         (Default: False)
    :argument Boolean DoCheckPoint: Create new check point after successful
         archival (Default: False)

    :argument Boolean DoOverwrite: Force overwrite of existing data with data
         in existing files (Default: False)

    :returns: None

    .. TODO :: Failed to detect symlink
    """
    logger = logging.getLogger('Archive.LoadEmUp')

    compress = tables.Filters(complevel = 9, complib = 'zlib')
    h5data = tables.openFile(h5filename, mode = 'a', filters = compress,
        title = 'H2Pc')
    if not h5data.isUndoEnabled():
        h5data.enableUndo(filters = compress)

    for root, subdirs, files in os.walk(path):
        for cwdfile in files:
            #Relative path-qualified filename
            filename = os.path.join(root, cwdfile)

            #Skip symlinks
            if os.path.islink(filename):
                continue

            ##############################
            # Parse CHARMM coordinate file
            ##############################

            if len(cwdfile) > 4 \
                and cwdfile[-4:] == '.cor' \
                and cwdfile[-6:-4] != '.d':

                # CARDs with Drude particles are skipped since they
                # will need to be re-equilibrated as part of the SCF process

                cwdh5name = cwdfile[:-4]
                location = os.path.join('/' + root, cwdh5name)

                Coords = OpenHDF5Table(h5data, location, 'CHARMM_CARD',
                    CHARMM_CARD, 'CHARMM CARD coordinates', DoOverwrite)

                if Coords is not None:
                    LoadCHARMM_CARD(Coords, cwdfile)

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
                    isRunning = CheckIfRunning(filename)
                    if isRunning:
                        #If it is running, do nothing
                        logger.info('Skipping incomplete calculation in %s',
                                    filename)
                        continue
                    
                    #TODO If it isn't, check for errors
                    
                    #If there are convergence errors:
                    RunQChemAgain(filename)

                    if NukeInc:
                        logger.info('Deleting incomplete calculation in \
directory: %s', root)
                        shutil.rmtree(root)
                        continue

                if '-nonpol' in calctype:
                    Environment = 'Fixed'
                else:
                    Environment = 'WithDrude'
                    #TODO Record Drude particle positions
                    logger.warning("Here I should be recording the Drude \
particle positions but I'm not!")
 
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
                    QCtddft = QChemTDDFTParser(filename)
                    TDAEnergies, TDADipole, Energies, Dipole = \
                        QCtddft.GetEnergiesAndDipole()
                    TDAMultiplicities, TDAOscillatorStrengths, \
                        Multiplicities, \
                        OscillatorStrengths = QCtddft.GetOtherData()

                    if TDAEnergies != None:
                        TDAEnergies /= kcal_mol
                    if Energies != None:
                        Energies /= kcal_mol
               
                elif cwdfile == 'qchem.out':
                    #Assume all energy calculations go through the
                    #CHARMM/Q-Chem interface.
                    #
                    #This avoids double processing of the polarizable
                    #simulations
                    continue


                ################################################
                # Record completed calculations into HDF5 file #
                ################################################

                #Update only if job is done
                if not isDone:
                    continue

                location = os.path.join('/' + GeometryName, Environment,
                                        'Sites', Site)

                #Update a single state's energy
                if Energy != None:
                    RecordIntoHDF5EArray(Energy * kcal_mol, h5data, location,
                                         'Energy', State,
                                         DoOverwrite = DoOverwrite)
                #Update all energies
                if Energies != None:
                    RecordIntoHDF5EArray(Energies * kcal_mol, h5data, location,
                                         'Energy', DoOverwrite = DoOverwrite)

                #Update a transition dipole
                if Dipole != None:
                    RecordIntoHDF5CArray(Dipole, h5data, location, 'Dipole',
                                         DoOverwrite = DoOverwrite)

                if 'tddft' in calctype: #Store also TDA energies and dipoles
                    RecordIntoHDF5EArray(TDAEnergies * kcal_mol, h5data,
                        location, 'TDAEnergy', DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(TDADipole, h5data, location,
                        'TDADipole', DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(TDAMultiplicities, h5data, location,
                        'TDAMultiplicities', atom = tables.atom.IntAtom(),
                        DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(Multiplicities, h5data, location,
                        'Multiplicities', atom = tables.atom.IntAtom(),
                        DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(TDAOscillatorStrengths, h5data,
                        location, 'TDAOscillatorStrengths',
                        DoOverwrite = DoOverwrite)
                    RecordIntoHDF5CArray(OscillatorStrengths, h5data, location,
                        'OscillatorStrengths', DoOverwrite = DoOverwrite)

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
            while (len(rootup.split(os.path.sep)) - \
                   len(path.split(os.path.sep))) > 1:
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

    LoadEmUp(args.workingdir, args.h5file, args.nuke, args.nukeinc,
             args.checkpoint, args.force)
