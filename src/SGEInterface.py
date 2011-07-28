#!/usr/bin/env python

"""
Generates QM/MM electronic structure jobs and submits them to a Sun Grid
Engine queue.
"""
import logging, os
from glob import glob

import SGE

#Compatibility kludge
try:
    from HDF5Interface import CHARMM_RTF_Short, ResidList, OpenHDF5Table, \
        CHARMM_CARD, LoadCHARMM_CARD
    import tables
except ImportError:
    pass

from JobHandler import SGEArrayTaskHandler
            

def PrepareJobs(runpack = '/home/cjh/rmt/runpack/*', 
    joblist = 'sgearraytasklist.txt', sge_jobname = 'h2pc-qmmm-rmt',
    DoOverwrite = False, h5filename = 'h2pc-data.h5'):

    """
    Prepares jobs to be executed in a term:`SGE` array job environment.

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

    logger = logging.getLogger('SGEInterface.PrepareJobs')

    if os.path.exists(joblist) and not DoOverwrite:
        print 'Skipping population of', joblist
        #count numjobs
        numjobs = 0
        for _ in open(joblist):
            numjobs += 1
    else:
        alljobs = []

        #Collect all data into HDF5 file
        h5data = tables.openFile(h5filename, 'w')

        #############################
        # Load in residues to do QM #
        #############################
        MolList = OpenHDF5Table(h5data, '/Model', 'MolList',
            ResidList, DoOverwrite)

        if MolList is not None:
            for res_id in open('../mol.list'):
                logger.info('Adding mol.list')
                try:
                    data = MolList.row
                    data['ResID'] = int(res_id)
                    data.append()
                except ValueError:
                    pass

        ###############################################
        # Load in atomic charges from CHARMM topology #
        ###############################################
        for CHARMM_RTFile in glob(runpack+'.rtf'):
            RTF = OpenHDF5Table(h5data, '/Model', 'CHARMM_RTF',
                CHARMM_RTF_Short, DoOverwrite, DoAppend = True)

            if RTF is not None:
                for line in open(CHARMM_RTFile):
                    logger.info('Adding CHARMM RTF: %s', CHARMM_RTFile)
                    t = line.split()
                    try:
                        data = RTF.row
                        data['Type'] = t[1]
                        data['Charge'] = float(t[3])
                        data.append()
                    except ValueError:
                        raise ValueError, 'Parse error in RTF:', CHARMM_RTFile

        #####################################
        # Load coordinates from CHARMM CARD #
        #####################################
        for coordfile in glob('*.cor'):
            logger.info('Adding CHARMM CARD: %s', coordfile)
            coordfileroot = coordfile.split('.')[0]

            Coords = OpenHDF5Table(h5data, coordfileroot, 'CHARMM_RTF',
                CHARMM_CARD, DoOverwrite)

            if Coords is not None:
                LoadCHARMM_CARD(Coords, coordfile)

            alljobs += iteratejobs(coordfileroot)

        f = open(joblist, 'w')
        for i, j in enumerate(alljobs):
            f.write(str(i+1)+' \t')
            f.write(' '.join(j) + '\n')
        f.close()

        numjobs = len(alljobs)
        logger.info('Indexed %d jobs', numjobs)

    #Now create SGE array job and submit self as an array job
    x = SGE.qsubOptions()
    x.args.N = sge_jobname
    x.args.t = '1-'+str(numjobs)
    x.args.o = '$JOB_NAME.$JOB_ID.$TASK_ID.stdout.log'.replace('$', '\$')
    x.args.command = 'python JobHandler.py --taskfile '+joblist
    x.args.b = 'y'
    print '#To submit the array job, run this command:'
    x.execute(mode = 'echo')


def iteratejobs(corfileroot, isPolarizable = False, jobtypes = ('tddft',),
        mollist = '../mol.list'):
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
    parser.add_argument('--overwrite', action = 'store_true', default = \
        False, help = 'Continue calculation even output file exists')

    parser.add_argument('--qcprog', action = 'store', default = \
        'qchem40', help = 'Name of Q-Chem binary to execute')

    parser.add_argument('--loglevel', action = 'store', default=logging.INFO,
            type = int, help = 'Logging level')

    args = parser.parse_args()

    logging.basicConfig(level = args.loglevel)

    if 'SGE_TASK_ID' not in os.environ:
        PrepareJobs(joblist = args.taskfile, DoOverwrite = args.overwrite)
    else: #In SGE array task
        SGEArrayTaskHandler(ge_task_id = os.environ['SGE_TASK_ID'],
            taskinfofile = args.taskfile, overwrite = args.overwrite,
            qchemcmd = args.qcprog)


