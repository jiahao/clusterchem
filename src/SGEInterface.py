#!/usr/bin/env python

"""
Generates QM/MM electronic structure jobs and submits them to a Sun Grid
Engine queue.
"""
import logging, os, sys

try:
    import tables
    HavePyTables = True
except ImportError:
    HavePyTables = False

from glob import glob
from OSUtils import chdirn
from subprocess import Popen, PIPE, STDOUT

import SGE
from OSUtils import chdirn
from QChemInterface import QChemInputForElectrostaticEmbedding, \
         QChemInputForTransitionDipole
from QChemIO import QChemInput

if HavePyTables:
    class CHARMM_RTF_Short(tables.IsDescription):
        """
        PyTables data structure containing some CHARMM RTF information
        pertaining to atomic charges
        """
        #: Atom Type
        Type = tables.StringCol(4)
        #: Atomic charge
        Charge = tables.Float64Col()

    class ResidList(tables.IsDescription):
        """
        PyTables data structure containing a list of molecule residue lists
        """
        #: Residue ID
        ResID = tables.UInt32Col()


def OpenHDF5Table(h5data, location, name, description, DoOverwrite = False, DoAppend = True):

    logger = logging.getLogger('OpenHDF5Table')
    try:
        h5table = h5data.createTable(name = name, where = location,
            title = name, description = description, createparents = True)
    except tables.exceptions.NodeError:
        if not DoOverwrite:
            logger.info('HDF5 table already exists: %s: skipping.',
                os.path.join(location, name))

            if not DoAppend:
                return None
            else:
                h5table = h5data.getNode(name = name, where = location)

        logger.info('Overwriting existing %s in HDF5: %s',
                   name, os.path.join(location, name))

        h5data.removeNode(location, name)

        h5table = h5data.createTable(name = name, where = location,
            title = name, description = description, createparents = True)

    return h5table



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
        from Archive import CHARMM_CARD, LoadCHARMM_CARD

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
            coordfileroot = coordfile.split('.')[0].replace('-', '_')

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
    x.args.command = 'python '+sys.argv[0]+' --taskfile '+joblist
    x.args.b = 'y'
    print '#To submit the array job, run this command:'
    x.execute(mode = 'echo')



def linkwild(src, dest, overwrite=False):
    """
    Symbolic link creation with wildcard support.
    
    @note Symlinks matched by src will *not* be copied to avoid possible
    circular references.

    @param src A specification of source files parseable by glob.glob()

    @param dest A destination (directory).

    @param overwrite: Whether to overwrite files that exist in :py:var:dest.
    """

    logger = logging.getLogger('linkwild')

    for fname in glob(src):
        destfname = os.path.join(dest, os.path.basename(fname))
        if overwrite:
            doit = True
        else:
            doit = not os.path.exists(destfname)
        if (os.path.isfile(fname) or os.path.islink(fname)) and doit:
            if overwrite and os.path.exists(destfname):
                os.unlink(destfname)
            try:
                os.symlink(fname, destfname)
                logger.info('Made symlink %s --> %s', os.path.abspath(fname),
                        os.path.abspath(destfname))
            except OSError:
                logger.warning('Warning: could not make symlink: %s -X-> %s',
                        os.path.abspath(fname), os.path.abspath(destfname))


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



def RunQMMM(infile, corfile, qchemcnt, res1, np = 1):
    """
    Executes the QM/MM wrapper script charmmpol

    If the outputfile exists, execution will be skipped.

    @param string infile: Name of :term:`CHARMM input file`
    @param string corfile: Name of :term:`CHARMM card file`
    @param string qchemcnt: Name of Q-Chem template file 
    @param integer np: Number of CPU cores to request (default: 1) \
    @param integer res1: Index of molecule in the QM region
    
    .. NOTE:: This currently overwrites the :envvar:`NP` environment variable

    .. @todo:: DOCUMENT TEMPLATES
    .. @todo:: This is currently a fancy wrapper around Lee-Ping, Shane and
        Tim's charmmpol script. REPLACE WITH INTERNAL CALL TO CHARMMPOL.PY
    """
    outfile = infile[:-2]+'out'
    corfile = corfile.upper()
    qchemcnt = qchemcnt.upper()
    
    logger = logging.getLogger('SGEInterface.RunQMMM')

    if not os.path.exists(outfile):
        os.environ['NP'] = str(np)
        p = Popen('charmmpol '+infile+' CORFILE='+ corfile +\
                ' QCHEMCNT='+qchemcnt+' np='+str(np)+ ' res1='+str(res1),
                stdout=PIPE, stderr=STDOUT, shell=True)
        stdout = p.stdout.read().strip()
        if stdout:
            logger.warning(stdout)
    else:
        logger.warning('Output file already exists: %s. Skipping calculation',
                os.path.abspath(outfile))



def domyjob(corfileroot, resid, jobtype, overwrite = False, qchemcmd = 'qchem'):
    """
    @todo DOCUMENT

    @param qchem Q-Chem binary or QCPROG for transition dipole calculation

    .. versionadded:: 0.1
    """

    #Historical note: prior to the installation of qchem40.beta,
    #the executable used was
    #qchem= '/home/tkowalcz/qchem/qchem-2ecoupling/linux64/exe/qcprog.exe'
    olddir = os.getcwd()

    #mkdir -p $corfileroot/$resid/$jobtype && cd [...]
    chdirn(corfileroot)
    linkwild('/home/cjh/rmt/runpack/*', '.')

    chdirn(str(resid))
    chdirn(jobtype)

    #Initialize
    linkwild(os.path.join('..', '..', '*'), '.')

    if '-nonpol' in jobtype:
        charmmcardfile = corfileroot+'.COR'
    else:
        charmmcardfile = corfileroot+'.D.COR'
    charmmcardfile = charmmcardfile.upper()

    outfiles = glob('*.out')
    if overwrite:
        for filename in outfiles:
            os.unlink(filename)
            print 'Deleted existing output file', filename

    rtfiles = ('h2pc.rtf', 'znpc.rtf')

    if jobtype == 'gs-nonpol':
        #RunQMMM('h2pc.in', charmmcardfile, 'QCHEM-GROUND.IN', resid)
        Q = QChemInput('QCHEM-GROUND.IN')
        Q = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                rtfiles, Q)
        Q.write('qchem.working.in')
        Q.execute('qchem.out', qchemcmd = qchemcmd, parse = False)
    elif jobtype == 'dscf-nonpol':
        #RunQMMM('h2pc.in', charmmcardfile, 'QCHEM-DSCF.IN', resid)
        Q = QChemInput('QCHEM-DSCF.IN')
        Q = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                rtfiles, Q)
        Q.write('qchem.working.in')
        Q.execute('qchem.out', qchemcmd = qchemcmd, parse = False)
    elif jobtype == 'tddft-nonpol':
        Q = QChemInput('QCHEM-TDDFT.IN')
        Q = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                rtfiles, Q)
        Q.write('qchem.working.in')
        Q.execute('qchem.out', qchemcmd = qchemcmd, parse = False)


    elif jobtype == 'gs':
        RunQMMM('h2pc-pol.in', charmmcardfile, 'QCHEM-GROUND.IN', resid)
    elif jobtype == 'dscf':
        RunQMMM('h2pc-pol.in', charmmcardfile, 'QCHEM-DSCF.IN', resid)
    elif jobtype == 'td':
        filenames = glob('../gs/*/qchem*working.in')\
                  + glob('../dscf/*/qchem*working.in')
        if len(filenames) == 2:
            print 'Running polarizable QM/MM transition dipole moment \
                    for resid =', resid
        Q = QChemInputForTransitionDipole(filenames[0], filenames[1])

        Q.filename = 'td.in'
        Q.jobs[1].rem_delete('purify')
        Q.execute(qchemcmd = qchemcmd, parse = False)
        if len(filenames) > 2:
            print 'Too many possible inputs:', ' '.join(filenames)

    elif jobtype == 'td-nonpol':
        print 'Running nonpolarizable QM/MM transition dipole moment for \
                resid =', resid

        filenames = glob('../gs-nonpol/*/qchem*working.in') \
                  + glob('../dscf-nonpol/*/qchem*working.in')
        if len(filenames) == 2:
            print 'Reusing input from', ', '.join([os.path.abspath(f) \
                    for f in filenames])
            Q = QChemInputForTransitionDipole(filenames[0], filenames[1])
        else:
            print 'Too few or too many files found:', \
                ', '.join([os.path.abspath(f) for f in filenames])
            print 'Regenerating state descriptions using default rem blocks'

            Q1 = QChemInput('QCHEM-GROUND.IN')
            Q1 = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                    'h2pc.rtf', Q1)

            Q2 = QChemInput('QCHEM-DSCF.IN')
            Q2 = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                    'h2pc.rtf', Q2)

            Q = QChemInputForTransitionDipole(Q1, Q2)

        Q.filename = 'td.in'
        Q.jobs[1].rem_delete('purify')
            
        Q.execute(qchemcmd = qchemcmd, parse = False)

    else:
        assert False, 'Unknown jobtype '+str(jobtype)

    #Return to old directory
    os.chdir(olddir)



def SGEArrayTaskHandler(ge_task_id, taskinfofile = 'sgearraytasklist.txt', 
        overwrite = False, qchemcmd = 'qchem'):
    """
    This gets called if the main program is run in a Sun Grid Engine array
    job environment.

    Looks up job assigned to array id ge_task_id that is specified in
    taskinfofile, loads the parameters saved in taskinfofile and passes them
    to domyjob().

    @param ge_task_id: the array job id provided by Sun Grid Engine. 
        This is actually an integer but the string representation is more
        convenient to to work with.

    @param taskinfofile: Name of file containing task information.
        (Default: 'sgearraytasklist.txt')
    
    The file named in taskinfofile is a text file containing data in the
    following format:
        Col 1: job_id (integer)
        Col 2: path_to_CHARMM_CARD_file (string)
        Col 3: residue_id (integer)
        Col 4: job_type (string)
    

    See domyjob() for valid values of job_type.

    @returns :rtype: None

    .. versionadded:: 0.1
    """
    for line in open(taskinfofile):
        t = line.split()
        if len(t) >= 4 and t[0] == ge_task_id:
            thiscorfileroot, thisresid, thisjobtype = t[1], t[2], t[3]
            break
    domyjob(thiscorfileroot, thisresid, thisjobtype, overwrite, qchemcmd)


            
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


