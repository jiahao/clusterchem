#!/usr/bin/env python
"""
Archives data into a HDF5 database.

This is a standalone executable requiring PyTables
and the homebrew QChemIO.

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
    parser.add_argument('-n', '--nuke', action = 'store_true', default = False,
                        help = 'Deletes files after parsing')
    parser.add_argument('--nukeinc', action = 'store_true', default = False,
                        help = 'Deletes incomplete output files')
    parser.add_argument('--checkpoint', action = 'store_true', default = True,
                        help = 'Create new Pytables checkpoint in HDF5 file')
    parser.add_argument('-f', '--force', action = 'store_true', default = False,
                        help = 'Forces reloading of existing data')
    parser.add_argument('--loglevel', action = 'store', default = 'info',
                help = 'Logging level (debug, info, warning, error, critical)')
    parser.add_argument('-l', '--logfile', action = 'store',
                        default = 'archival.log',
                        help = 'Name of log file to write (Default: archival.log; specify None to disable')
    args = parser.parse_args()



def CheckIfRunning(filename = None, taskfile = args.taskfile):
    """
    Is the process responsible for producing filename still running?

    Need to infer this from qstat

    Get list of running SGE processes
    then check working directories
    """

    logger = logging.getLogger('Archive.CheckIfRunning')

    #If taskfile doesn't exist, assume it's not running
    try:
        open(taskfile).close()
    except IOError:
        logger.warning("""Cannot open taskfile %s
Assuming not running anymore""", taskfile)
        return False

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



def RunQChemAgain(filename, inputfile = None):
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


    If input file is None, then assume same directory as filename + qchem.working.in

    returns False if nothing was done, True if a new input was written.
    """

    #TODO
    logger = logging.getLogger('Archive.RunQChemAgain')

    #What was the error?
    mode = 'seek'
    skiplines = 0
    error = ''
    for line in open(filename):
        if skiplines > 0:
            skiplines -= 1
            continue

        if mode == 'seek':
            if 'Q-Chem fatal error occurred' in line:
                skiplines = 1
                mode = 'readerror'
                error = ''
        elif mode == 'readerror':
            if line.strip() == '':
                break #Assume no more error to read if blank line found; done
            else:
                error += line

    if inputfile is None:
        inputfile = os.path.join(os.path.dirname(filename), 'qchem.working.in')

    try:
        with open(inputfile) as f:
            pass
    except IOError:
        logger.error('Error detected in output file %s but input file %s \
was not found\nRun again manually: %s', filename, inputfile, filename)
        return False

    error = error.strip()

    if error == 'MaxIt Reached in CisIt0':
        #Need to increase MAX_CIS_CYCLES
        newinput = QChemInput(inputfile)
        current_value = newinput.GetCurrentJob().rem_get('max_cis_cycles')
        if current_value is None: #Assume it's 30
            current_value = 30

        #Double it and try again
        newinput.GetCurrentJob().rem_set('max_cis_cycles', str(2*current_value))
        newinput.write()
        return True
    elif error == '':
        #Execution terminated early
        #Is there a coredump?
        try:
            open('core', 'rb').close()
            #There is a coredump
            #Now what?
        except IOError: #No file
            #Terminated early with no good reason?
            pass
            
    else:
        logger.critical("""Error in %s was: %s
I don't know how to deal with this error.
Run again manually: %s""", filename, error, inputfile)


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
    logger.info('Opening HDF5 file: %s', h5filename)

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
            if len(cwdfile) > 6 and cwdfile[-6:] == '.d.cor':
                #Here is a CARD with Drude particles
                #does an equivalent card without Drudes exist?
                newfile = filename[:-5]+'cor'
                if os.path.exists(newfile):
                    continue

                from CHARMMUtil import dedrude
                targetfile = os.path.join(root, newfile)
                dedrude(filename, newfile)
                logging.info('Stripped Drude particles from %s', cwdfile)
                
                cwdh5name = cwdfile[:-6]
                location = os.path.join('/' + root, cwdh5name)

                Coords = OpenHDF5Table(h5data, location, 'CHARMM_CARD',
                    CHARMM_CARD, 'CHARMM CARD coordinates', DoOverwrite)

                if Coords.nrows < 1:
                    LoadCHARMM_CARD(Coords, targetfile)
                    
            elif len(cwdfile) > 4 \
                and cwdfile[-4:] == '.cor' \
                and cwdfile[-6:-4] != '.d':
                # CARDs with Drude particles are skipped since they
                # will need to be re-equilibrated as part of the SCF process

                cwdh5name = cwdfile[:-4]
                location = os.path.join('/' + root, cwdh5name)

                Coords = OpenHDF5Table(h5data, location, 'CHARMM_CARD',
                    CHARMM_CARD, 'CHARMM CARD coordinates', DoOverwrite)

                if Coords.nrows < 1:
                    LoadCHARMM_CARD(Coords, filename)
              
            ####################################
            # Parse CHARMM or Q-Chem output file
            ####################################

            if len(cwdfile) > 4 and cwdfile[-4:] == '.out':
                #Assume path has the structure [.../]GeomName/Site/jobtype
                Site = root.split(os.sep)[-2]
                GeometryName = root.split(os.sep)[-3]
                calctype = os.path.split(root)[1]

                State = Energy = Energies = Dipole = Environment = None
                try:
                    data = ParseOutput(filename)
                except IOError: #file gone
                    logger.warning('File %s disappeared before it could be \
parsed.')
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
                    if QCtddft is None: continue #Incomplete file
                    TDAEnergies, TDADipole, Energies, Dipole = \
                        QCtddft.GetEnergiesAndDipole()
                    TDAMultiplicities, TDAOscillatorStrengths, \
                        Multiplicities, \
                        OscillatorStrengths = QCtddft.GetOtherData()

                    if TDAEnergies != None:
                        TDAEnergies /= kcal_mol
                    if Energies != None:
                        Energies /= kcal_mol
                        print Energies               
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
                    try:  
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
                    except TypeError: #No data returns None - something here will complain
                                      #usually the unit conversion from kcal_mol
                        pass
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

            num_errors = 0
            error_jobs = set()
            for filename in filenames:
                haserror = False
                for line in open(filename):
                    if 'Traceback (most recent call last):' in line:
                        #Uh-oh, Python crashed
                        haserror = True
                        error_jobs.add(filename.split('.')[1])
                        break
                if not haserror:
                    os.unlink(filename)
                else:
                    num_errors += 1

            if num_errors > 0:
                logger.warning("""Found %d jobs that crashed!
Holding rest of job array: jobids %s.
Use qalter -h U -u $USER to resume.""", num_errors, ' '.join(error_jobs))
                from subprocess import Popen, PIPE, STDOUT
                Popen('qhold '+' '.join(error_jobs), shell=True)

    if h5data.isUndoEnabled() and DoCheckPoint:
        mark = h5data.mark()
        logger.info('Checkpointing HDF5 database at mark point %d', mark)

    h5data.close()



if __name__ == '__main__':
    import sys
    logger = logging.getLogger('Archive')
    logger.setLevel(getattr(logging, args.loglevel.upper()))
    if args.logfile != 'None':
        #Add another log handler to duplicate output to file and stdout
        logfile = logging.FileHandler(args.logfile)
        logger.addHandler(logfile)

    LoadEmUp(args.workingdir, args.h5file, args.nuke, args.nukeinc,
             args.checkpoint, args.force)
