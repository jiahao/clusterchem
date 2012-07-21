#!/usr/bin/env python

"""
Generates QM/MM electronic structure jobs and submits them to a Sun Grid
Engine queue.
"""
import logging, os, sys, warnings
from glob import glob

from CHARMMUtil import dedrude
import SGE

#Compatibility kludge
try:
    from HDF5Interface import CHARMM_RTF_Short, ResidList, OpenHDF5Table, \
        CHARMM_CARD, LoadCHARMM_CARD
    import tables
    #Turn off NaturalNameWarning from PyTables
    warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

except ImportError:
    pass

from tables.nodes import filenode

from JobHandler import SGEArrayTaskHandler
            
def PrepareJobs(runpack = '/home/cjh/rmt/runpack/*', 
    joblist = 'sgearraytasklist.txt', sge_jobname = 'qm-mm',
    DoOverwrite = False, h5filename = 'qm-mm.h5'):

    """
    Prepares jobs to be executed in a term:`SGE` array job environment.

    This expects the following to be in the current working directory OR already
    present in the HDF5 data file:

    - At least one CHARMM CARD coordinate file or its corresponding CHARMM_CARD table
    - At least one CHARMM RTF file containing data on atomic charges.
    - At least one Q-Chem input template named QCHEM-*.IN.
    - A text file called 'mol.list' containing a list of all the CHARMM resids
      for which electronic structure calculations should be generated.
      The format is simply one number per line.

    This does the following:
        1. Scans all files in the current directory for files matching
           @c \*.d.cor, which will be assumed to be :term:`CHARMM card file` s
           containing coordinates for the molecular system, and also Drude
           particles. For each such file, it will generate the corresponding
           :term:`CHARMM card file` s without Drude particles using dedrude().

        1. Copies files specified by runpack using linkwild() to a
           directory named after the root of the :term:`CHARMM card file` name.
           If this directory does not exist, it will be created.

        1. Creates a symlink to the :term:`CHARMM card file`. If necessary, the
           filename will be converted into upper case due to a limitation of
           CHARMM.

        1. Writes data needed by :py:func:`domyjob` to a file named in joblist.

        1. Generates the command line needed to submit the array job to a Sun
           Grid Engine environment.

    This function requires the homebrew :py:mod:`SGE` module.

    .. versionadded:: 0.1
    """
    logger = logging.getLogger('SubmitJob.PrepareJobs')

    ###################################################################
    # 1. Check whether data file has all the necessary information to #
    #    run calculations.                                            #
    ###################################################################

    if not os.path.exists(h5filename):
        logging.info('Creating new HDF5 file: %s', h5filename)
        h5data = tables.openFile(h5filename,'a')
    else:
        h5data = tables.openFile(h5filename, 'r')
    
    logging.debug('HDF5 file opened: %s', h5filename)
    ToLoad = list()

    # A. Does file contain a non-empty MolList?
    #    MolList is the list of molecule (CHARMM residue IDs) to generate
    #    QM/MM calculations for.
    try:
        numMols = h5data.getNode('/Model', 'MolList').nrows
        assert numMols > 0
        logger.info('QM/MM calculations will be generated for %d molecules',
                numMols)
    except (AssertionError, tables.NoSuchNodeError):
        logger.info('No MolList found.')
        ToLoad.append('MolList')

    # B. Are there CHARMMM CARD coordinates in the current directory that
    #    are not already in the HDF5 file?

    #This code will fail if there is this file.
    assert not os.path.exists("Model.cor"), "Name collision with internal data table /Model! I don't know how to deal with this."
    CHARMM_CARDFiles=glob('*.cor')
    InternalCHARMM_CARDs = [entry._v_name+'.cor' for entry in h5data.listNodes('/')]
    ToLoadCHARMM_CARDs = set(CHARMM_CARDFiles) - set(InternalCHARMM_CARDs)

    for CHARMM_CARDSkipped in set(InternalCHARMM_CARDs) & set(CHARMM_CARDFiles):
        logger.info('Skipped existing CHARMM CARD: %s', CHARMM_CARDSkipped)
    if len(ToLoadCHARMM_CARDs) > 0:
        logger.info('Found CHARMM CARDs to load.')
        for crd in ToLoadCHARMM_CARDs:
            logger.debug('Will load CHARMM CARD: %s', crd)
    else:
        logger.info('Did not find any CHARMM CARDs to load.')

    # C. Does the file already contain charge data from CHARMM RTF files?
    try:
        h5rtf = h5data.getNode('/Model', 'CHARMM_RTF')
        assert h5rtf.nrows > 0
        logger.info('CHARMM RTF data found in HDF5 file, skipping any RTFs present in directory')
    except (AssertionError, tables.NoSuchNodeError):
        logger.info('No RTF found.')
        ToLoad.append('RTF')
    
    # D. Do we have Q-Chem input templates in the current director that?
    #    are not already in the HDF5 file?
    QChemInputTemplates=glob('QCHEM-*.IN')
    try:
        InternalTemplates = ['QCHEM-'+entry._v_name+'.IN' for entry in h5data.listNodes('/Model/QChemTemplates')]
    except tables.NoSuchNodeError:
        InternalTemplates = []
    ToLoadTemplates = set(QChemInputTemplates) - set(InternalTemplates)

    for TemplateSkipped in set(InternalTemplates) & set(QChemInputTemplates) :
        logger.info('Skipped existing Q-Chem input template: %s', TemplateSkipped)
    if len(ToLoadTemplates) > 0:
        logger.info('Found Q-Chem input templates to load.')
        for crd in ToLoadTemplates:
            logger.debug('Will load Q-Chem template: %s', crd)
    else:
        logger.info('Did not find any Q-Chem input templates to load.')

    # E. Does the file already contain a list of jobs to be run?
    try:
        numJobs = h5data.getNode('/Model', 'JobsToRun').nrows
        assert numJobs > 0
        logger.info('Existing job schedule found.')
    except (AssertionError, tables.NoSuchNodeError):
        logger.info('No job schedule found. Will generate a new one.')
        ToLoad.append('Jobs')
    
    h5data.close()
    
    ###################################################################
    # 2. Load up necessary data                                       #
    ###################################################################
    
    #XXX HOW TO LOCK FILE??? WHAT IF THERE ARE JOBS RUNNING?
    h5data = tables.openFile(h5filename, 'a')

    # A. Load up MolList
    if 'MolList' in ToLoad:
        assert os.path.exists('mol.list'), 'No mol.list found'

        MolList = OpenHDF5Table(h5data, '/Model', 'MolList',
            ResidList, 'List of molecule ids for which QM data exist',
            DoOverwrite)

        logger.info('Adding mol.list')
        for res_id in open('mol.list'):
            try:
                data = MolList.row
                data['ResID'] = int(res_id)
                data.append()
            except ValueError:
                pass
        

    # B. Load up CHARMM CRDs
    alljobs = []
    for coordfile in ToLoadCHARMM_CARDs:
        logger.info('Adding CHARMM CARD: %s', coordfile)
        coordfileroot = coordfile.split('.')[0]

        Coords = OpenHDF5Table(h5data, '/'+coordfileroot, 'CHARMM_CARD',
            CHARMM_CARD, 'CHARMM CARD File', DoOverwrite)

        if Coords is not None:
            LoadCHARMM_CARD(Coords, coordfile)

        alljobs += iteratejobs(coordfileroot)

    # C. Load up charge data in CHARMM RTF
    if 'RTF' in ToLoad:
        for CHARMM_RTFile in glob('*.rtf'):
            logger.info('Adding CHARMM RTF: %s', CHARMM_RTFile)
            RTF = OpenHDF5Table(h5data, '/Model', 'CHARMM_RTF',
                CHARMM_RTF_Short, 'CHARMM Topology File', DoOverwrite, DoAppend = True)

            assert RTF is not None
            for line in open(CHARMM_RTFile):
                t = line.split()
                if len(t) > 0 and t[0].upper() != 'ATOM': continue 
                try:
                    data = RTF.row
                    data['Type'] = t[1]
                    data['Charge'] = float(t[3])
                    logger.debug('Atom type %s has charge %12.8f', t[1], float(t[3]))
                    data.append()
                except (ValueError, IndexError):
                    pass
            break #REMOVE ME
        if len(glob('*.rtf'))==0:
            raise ValueError, """No RTF files found in current directory.

Cannot continue without atomic charges for QM/MM electrostatic embedding."""


    # D. Load up Q-Chem input templates
    for QChemTemplate in ToLoadTemplates:
        try:
            h5data.getNode(where='/', name='Model')
        except tables.NoSuchNodeError:
            h5data.createGroup(where='/', name='Model')
        try:
            h5data.getNode(where='/Model', name='QChemTemplates')
        except tables.NoSuchNodeError:
            h5data.createGroup(where='/Model', name='QChemTemplates')

        fnode = filenode.newNode(h5data, where='/Model/QChemTemplates', 
                name=QChemTemplate[6:-3])
                            #0123456..-321
                            #QCHEM-XXXX.IN
        fnode.write(open(QChemTemplate).read())
        fnode.close()

    # E. Regenerate jobs list
    #TODO clean up - this uses alljobs which was populated in 2B
    if 'Jobs' in ToLoad:
        #Populate text file
        if os.path.exists(joblist) and not DoOverwrite:
            logger.info('Skipping population of %s', joblist)
            #count numjobs
            numjobs = 0
            for line in open(joblist):
                len_line = len(line.split())
                if len_line == 4:
                    numjobs += 1
                elif len_line != 0:
                    raise ValueError, 'Invalid format of line:' + line
        else:
            with open(joblist, 'w') as f:
                for i, j in enumerate(alljobs):
                    f.write(str(i+1)+' \t')
                    f.write(' '.join(j) + '\n')
    
            numjobs = len(alljobs)
            logger.info('Indexed %d jobs', numjobs)

        #Populate HDF5 file
        class SGEArrayJobDispatch(tables.IsDescription):
            Geometry = tables.StringCol(80)
            ResID = tables.UInt32Col()
            JobType = tables.StringCol(20) #TODO this really should be an enumerated type
            def __init__(self):
                tables.IsDescription.__init__(self)

        h5job = OpenHDF5Table(h5data, '/Model', 'JobList',
            SGEArrayJobDispatch, 'SGE Array job dispatching information',
            DoOverwrite)
        assert h5job is not None
        
        for idx, (crd, resid, jobtype) in enumerate(alljobs):
            data = h5job.row
            data['Geometry'] = crd
            data['ResID'] = int(resid)
            data['JobType'] = jobtype
            data.append()
            logger.debug('Job taskID %7d: %s %s %s', idx+1, crd, resid, jobtype)
        logger.info('Enumerated %d jobs for SGE job array', len(alljobs))

    h5data.close()
   
    ###################################################################
    # 3. Verify input                                                 #
    ###################################################################
    #TODO  
    # AB. Check that MolList actually iterates over molecules that exist
    #     in each CHARMM CRD
    # C.  Check that CHARMM RTF specifies charges for ALL atom types
    # TODO Import from JobHandler.py
    # D. Run a sample Q-Chem input to make sure it is not insane
    # E. Regenerate jobs list

    ###################################################################
    # 4. Submit job to the queue                                      #
    ###################################################################

    #Now create SGE array job and submit self as an array job
    x = SGE.qsubOptions()
    x.args.N = sge_jobname
    x.args.t = '1-'+str(numjobs)
    x.args.p = '-500'
    x.args.o = '$JOB_NAME.$JOB_ID.$TASK_ID.stdout.log'.replace('$', '\$')
    FullPathToJobHandler = os.path.dirname(sys.argv[0])+os.sep+'JobHandler.py'
    x.args.command = 'python '+FullPathToJobHandler+' --taskfile '+joblist
    x.args.b = 'y'
    print '#To submit the array job, run this command:'
    x.execute(mode = 'echo')


def iteratejobs(corfileroot, isPolarizable = False, jobtypes = ('tddft',),
        mollist = 'mol.list'):
    """
    Creates a hierarchy of subdirectories in corefileroot
    of the form 'resid/jobtype' and return information about the jobs created.

    Iterates over residues in mollist, then over a hardcoded set of jobtypes.

    The code for populating individual directories is currently commented out.

    If a file ending in .out is found at any time during the traversal, an
    empty list will be returned immediately.

    jobtypes is an iterable with the following strings as valid elements:
    
    tddft - Time-dependent density functional theory
    Additional experimental job types:
    gs - Ground state calculation only
    dscf - Delta SCF calculation for excited state
    td   - CDFT-CI transition dipole calculation using a gs and a dscf output

    For each of these there is an additional variant - tddft-nonpol,
    td-nonpol, etc.  each of which specifies a nonpolarizable environment.
    tddft, td, etc. refers to the default polarizable calculation. The -nonpol
    suffix will be automatically added if isPolarizable = False.

    @returns Parameters for py:func:`domyjob`    
    :rtype: list of tuples (_, 3)

    .. versionadded:: 0.1
    """

    logger = logging.getLogger('iteratejobs')

    resids = []
    for l in open(mollist):
        try:
            resid = int(l)
            resids.append(resid)
        except ValueError:
            pass
    jobdata = []

    #Initialize calculations

    for job in jobtypes:
        if not isPolarizable:
            job += '-nonpol'
        for resid in resids:
            #If no output files exist
            if len(glob(os.path.join(corfileroot, str(resid), job,
                '*.out'))) == 0:
                jobdata.append((corfileroot, str(resid), job))
            else:
                logger.warning('Output file already exists: %s. Aborting', 
                    ' '.join(glob(os.path.join(corfileroot, str(resid), job,
                        '*.out'))))
                return [] #Abort if somewhere, an output file exists.

    return jobdata

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'SGE Array Job manager')
    parser.add_argument('--taskfile', action = 'store', default = \
        'sgearraytasklist.txt', help = 'Name of SGE array task info file')
    parser.add_argument('--h5file', action = 'store', default = \
        'qm-mm.h5', help = 'Name of HDF5 database file to create or open')
    parser.add_argument('--overwrite', action = 'store_true', default = \
        False, help = 'Continue calculation even output file exists')

    parser.add_argument('--qcprog', action = 'store', default = \
        'qchem40.beta', help = 'Name of Q-Chem binary to execute')

    parser.add_argument('--loglevel', action = 'store', default='INFO',
            type = str, help = 'Logging level')

    args = parser.parse_args()

    logging.basicConfig()
    level = getattr(logging, args.loglevel.upper(), None)
    logging.getLogger().setLevel(level)
    if 'SGE_TASK_ID' not in os.environ:
        PrepareJobs(joblist = args.taskfile, DoOverwrite = args.overwrite, h5filename = args.h5file)
    else: #In SGE array task
        SGEArrayTaskHandler(ge_task_id = os.environ['SGE_TASK_ID'],
            taskinfofile = args.taskfile, overwrite = args.overwrite,
            qchemcmd = args.qcprog)


