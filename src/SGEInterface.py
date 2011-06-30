#!/usr/bin/env python

"""
Generates QM/MM electronic structure jobs and submits them to a Sun Grid Engine queue.

.. versionadded:: 0.1
"""
import logging, os, sys
from glob import glob

def PrepareJobs(runpack = '/home/cjh/rmt/runpack/*', joblist = 'sgearraytasklist.txt', 
    sge_jobname = 'h2pc-qmmm-rmt', overwrite = False):

    """
    Prepares jobs to be executed in a term:`SGE` array job environment.

    This does the following:
        1. Scans all files in the current directory for files matching \*.d.cor,
           which will be assumed to be :term:`CHARMM card file` s containing coordinates
           for the molecular system, and also Drude particles. For each such file,
           it will generate the corresponding :term:`CHARMM card file` s without Drude particles
           using :py:func:`dedrude`.

        1. Copies files specified by runpack using :py:func:`copywild` to a
           directory named after the root of the :term:`CHARMM card file` name. If this
           directory does not exist, it will be created.

        1. Creates a symlink to the :term:`CHARMM card file`. If necessary, the
           filename will be converted into upper case due to a limitation of
           CHARMM.

        1. Writes data needed by :py:func:`domyjob` to a file named in joblist.

        1. Generates the command line needed to submit the array job to a Sun
           Grid Engine environment.

    This function requires the homebrew :py:mod:`SGE` module.

    .. versionadded:: 0.1
    """

    from CHARMMUtil import dedrude
    from OSUtils import chdirn

    if os.path.exists(joblist) and not overwrite:
        print 'Skipping population of', joblist

        #count numjobs
        numjobs = 0
        for _ in open(joblist):
            numjobs += 1
    else:
        alljobs = []

        for coordfile in glob('*.cor'):
            coordfileroot = coordfile[:-4]

            isPolarizable = False
            if coordfileroot[-2:] == '.d': #Has Drude particles
                isPolarizable = True
                #Generates coordinate file without Drude particles
                nonpolcoordfile = coordfileroot[:-2] +'.cor'
                if not os.path.exists(nonpolcoordfile):
                    dedrude(coordfile, nonpolcoordfile)

                #XXX Monkey patch to hijack polarizable job generation with non-polarizable ones
                isPolarizable = False
                coordfileroot = coordfileroot[:-2]

            if True or not os.path.exists(coordfileroot):
                print 'Preparing', coordfileroot
                chdirn(coordfileroot)

                #Copy runpack to the coordfileroot
                copywild(runpack, '.')

                #Copy coordinate files
                #NB Due to a limitation of CHARMM's environment variable
                # parsing mechanism, the coordinate filenames
                # must be in upper case.
                if not os.path.exists(coordfile.upper()): 
                    os.symlink(os.path.join('..', coordfile), coordfile.upper())
                if not os.path.exists(nonpolcoordfile.upper()):
                    os.symlink(os.path.join('..', nonpolcoordfile), \
                               nonpolcoordfile.upper())

                #makejobs(coordfileroot)
                os.chdir('..')

            alljobs += iteratejobs(coordfileroot, isPolarizable)

        f = open(joblist, 'w')
        for i, j in enumerate(alljobs):
            f.write(str(i+1)+' \t')
            f.write(' '.join(j) + '\n')
        f.close()

        numjobs = len(alljobs)
        print 'Indexed', numjobs, 'jobs'


    #Now create SGE array job and submit self as an array job
    import SGE
    x = SGE.qsubOptions()
    x.args.N = sge_jobname
    x.args.t = '1-'+str(numjobs)
    x.args.o = '$JOB_NAME.$JOB_ID.$TASK_ID.stdout.log'.replace('$', '\$')
    x.args.command = 'python '+sys.argv[0]+' --taskfile '+joblist
    x.args.b = 'y'
    #x.args.v = 'PYTHONPATH=$PYTHONPATH:'+os.path.dirname(sys.argv[0])
    print '#To submit the array job, run this command:'
    x.execute(mode = 'echo')



def copywild(src, dest, overwrite=False, Verbose = True):
    """
    Copying subroutine that supports wildcards.

    NOTE. This routine does not actually copy the files physically, but
    instead creates symlinks.

    NOTE. Symlinks matched by src will *not* be copied to avoid possible
    circular references.

    :keyword string src: A specification of source files parseable by \
        :py:mod:`glob`.
    :keyword string dest: A destination (directory).
    :keyword Boolean overwrite: Whether to overwrite files that exist in \
        :py:var:dest.
    :keyword Boolean Verbose: Whether to print additional output

    .. versionadded:: 0.1
    """
    for fname in glob(src):
        destfname = os.path.join(dest, os.path.basename(fname))
        if overwrite:
            doit = True
        else:
            doit = not os.path.exists(destfname)
        if (os.path.isfile(fname) or os.path.islink(fname)) and doit:
            #if not os.path.exists(destfname) or (overwrite and os.path.exists(destfname)):
            if overwrite and os.path.exists(destfname): os.unlink(destfname)
            try:
                os.symlink(fname, destfname)
                if Verbose: print os.path.abspath(fname), '-->', os.path.abspath(destfname)
            except OSError:
                print 'Warning: could not make symlink. ', os.path.abspath(fname), '-X->', os.path.abspath(destfname)


def iteratejobs(corfileroot, isPolarizable = False, jobtypes = ('tddft',), mollist = '../mol.list'):
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

    For each of these there is an additional variant - tddft-nonpol, td-nonpol, etc.
    each of which specifies a nonpolarizable environment. tddft, td, etc. refers to
    the default polarizable calculation. The -nonpol suffix will be automatically
    added if isPolarizable = False.

    :returns: Parameters for py:func:`domyjob`    
    :rtype: list of tuples (_, 3)

    .. versionadded:: 0.1
    """
    resids = []
    for l in open(mollist):
        try:
            resid = int(l)
            resids.append(resid)
        except ValueError: pass
    jobdata = []

    #Initialize calculations

    for job in jobtypes:
        if not isPolarizable: job += '-nonpol'
        for resid in resids:
            #If no output files exist
            if len(glob(os.path.join(corfileroot, str(resid), job, '*.out'))) == 0:
                #copywild(os.path.join(corfileroot,'*'), os.path.join(corfileroot, str(resid), job))
                jobdata.append((corfileroot, str(resid), job))
            else:
                print 'Exists:', ' '.join(glob(os.path.join(corfileroot, str(resid), job, '*.out')))
                return [] #Abort if somewhere, an output file exists.

    return jobdata



def RunQMMM(infile, corfile, qchemcnt, res1, np = 1):
    """
    Executes the QM/MM wrapper script charmmpol

    If the outputfile exists, execution will be skipped.

    :param string infile: Name of :term:`CHARMM input file`
    :param string corfile: Name of :term:`CHARMM card file`
    :param string qchemcnt: Name of Q-Chem template file 
    :param integer np: Number of CPU cores to request (default: 1) \
    :param integer res1: Index of molecule in the QM region
    
    .. NOTE:: This currently overwrites the :envvar:`NP` environment variable

    .. TODO:: DOCUMENT TEMPLATES
    .. TODO:: This is currently a fancy wrapper around Lee-Ping, Shane and Tim's charmmpol script. REPLACE WITH INTERNAL CALL TO CHARMMPOL.PY
    """
    from subprocess import Popen, PIPE
    outfile = infile[:-2]+'out'
    corfile = corfile.upper()
    qchemcnt = qchemcnt.upper()
    
    if not os.path.exists(outfile):
        os.environ['NP'] = str(np)
        p = Popen('charmmpol '+infile+' CORFILE='+ corfile +' QCHEMCNT='+qchemcnt+' np='+str(np)+ ' res1='+str(res1), stdout=PIPE, shell=True)
        print p.stdout.read()
    else:
        print 'Output file already exists:', os.path.abspath(outfile)
        print 'Skipping calculation'



def domyjob(corfileroot, resid, jobtype, overwrite = False, qchemcmd = 'qchem'):
    """
    .. TODO:: DOCUMENT

    :param string qchem: Q-Chem binary or QCPROG for transition dipole calculation
    Calls :py:func:`RunQMMM`

    .. versionadded:: 0.1
    """

    #Historical note: prior to the installation of qchem40.beta, the executable used was
    #qchem= '/home/tkowalcz/qchem/qchem-2ecoupling/linux64/exe/qcprog.exe'
    from OSUtils import chdirn
    try:
        from QChemInterface import QChemInputForElectrostaticEmbedding, QChemInputForTransitionDipole
    except ImportError:pass

    from QChemIO import QChemInput

    olddir = os.getcwd()

    #mkdir -p $corfileroot/$resid/$jobtype && cd [...]
    chdirn(corfileroot)
    copywild('/home/cjh/rmt/runpack/*', '.')

    chdirn(str(resid))
    chdirn(jobtype)

    #Initialize
    copywild(os.path.join('..','..','*'), '.')

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
        Q = QChemInputForElectrostaticEmbedding(resid, charmmcardfile, rtfiles, Q)
        Q.write('qchem.working.in')
        Q.execute('qchem.out', qchemcmd = qchemcmd, parse = False)
    elif jobtype == 'dscf-nonpol':
        #RunQMMM('h2pc.in', charmmcardfile, 'QCHEM-DSCF.IN', resid)
        Q = QChemInput('QCHEM-DSCF.IN')
        Q = QChemInputForElectrostaticEmbedding(resid, charmmcardfile, rtfiles, Q)
        Q.write('qchem.working.in')
        Q.execute('qchem.out', qchemcmd = qchemcmd, parse = False)
    elif jobtype == 'tddft-nonpol':
        #Do not use CHARMM interface
        Q = QChemInput('QCHEM-TDDFT.IN')
        Q = QChemInputForElectrostaticEmbedding(resid, charmmcardfile, rtfiles, Q)
        Q.write('qchem.working.in')
        Q.execute('qchem.out', qchemcmd = qchemcmd, parse = False)


    elif jobtype == 'gs':
        RunQMMM('h2pc-pol.in', charmmcardfile, 'QCHEM-GROUND.IN', resid)
    elif jobtype == 'dscf':
        RunQMMM('h2pc-pol.in', charmmcardfile, 'QCHEM-DSCF.IN', resid)
    elif jobtype == 'td':
        try:
            os.environ['QCPROG'] = _qchem_td_qcprog
        except NameError: #_qchem_td_qcprog not defined
            pass

        filenames = glob('../gs/*/qchem*working.in')\
                  + glob('../dscf/*/qchem*working.in')
        if len(filenames) == 2:
            print 'Running polarizable QM/MM transition dipole moment for resid =', resid
        Q = QChemInputForTransitionDipole(filenames[0], filenames[1])

        Q.filename = 'td.in'
        Q.jobs[1].rem_delete('purify')
        Q.execute(qchemcmd = qchemcmd, parse = False)
        if len(filenames) > 2:
            print 'Too many possible inputs:', ' '.join(filenames)

    elif jobtype == 'td-nonpol':
        try:
            os.environ['QCPROG'] = _qchem_td_qcprog
        except NameError:
            pass

        print 'Running nonpolarizable QM/MM transition dipole moment for resid =', resid

        filenames = glob('../gs-nonpol/*/qchem*working.in') \
                  + glob('../dscf-nonpol/*/qchem*working.in')
        if len(filenames) == 2:
            print 'Reusing input from', ', '.join([os.path.abspath(f) for f in filenames])
            Q = QChemInputForTransitionDipole(filenames[0], filenames[1])
        else:
            print 'Too few or too many files found:', ', '.join([os.path.abspath(f) for f in filenames])
            print 'Regenerating state descriptions using default rem blocks'

            Q1 = QChemInput('QCHEM-GROUND.IN')
            Q1 = QChemInputForElectrostaticEmbedding(resid, charmmcardfile, 'h2pc.rtf', Q1)

            Q2 = QChemInput('QCHEM-DSCF.IN')
            Q2 = QChemInputForElectrostaticEmbedding(resid, charmmcardfile, 'h2pc.rtf', Q2)

            Q = QChemInputForTransitionDipole(Q1, Q2)

        Q.filename = 'td.in'
        Q.jobs[1].rem_delete('purify')
            
        Q.execute(qchemcmd = qchemcmd, parse = False)


    else:
        assert False, 'Unknown jobtype '+str(jobtype)

    #Return to old directory
    os.chdir(olddir)



def SGEArrayTaskHandler(ge_task_id, taskinfofile = 'sgearraytasklist.txt', overwrite = False, qchemcmd = 'qchem'):
    """
    This gets called if the main program is run in a Sun Grid Engine array job environment.

    Looks up job assigned to array id ge_task_id that is specified in taskinfofile, loads the parameters
    saved in taskinfofile and passes them to :py:func:`domyjob`.

    :param string ge_task_id: the array job id provided by Sun Grid Engine. \
    This is actually an integer but the string representation is more convenient to to work with.

    :param string taskinfofile: Name of file containing task information. (Default: 'sgearraytasklist.txt')
    
    The file named in taskinfofile is a text file containing data in the following format:
        job_id (integer)   path_to_CHARMM_CARD_file (string)   residue_id (integer)    job_type (string)
        ...

    See :py:func:`domyjob` for valid values of job_type.

    :returns: :rtype: None

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
    parser = argparse.ArgumentParser(description = 'SGE Array Job generator and handler')
    parser.add_argument('--taskfile', action = 'store', default = \
        'sgearraytasklist.txt', help = 'Name of SGE array task information file')
    parser.add_argument('--overwrite', action = 'store_true', default = \
        False, help = 'Continue calculation even output file exists')

    parser.add_argument('--qcprog', action = 'store', default = \
        'qchem40', help = 'Name of Q-Chem binary to execute')

    parser.add_argument('--loglevel', action = 'store', default = logging.INFO, type = int, help = 'Logging level')

    args = parser.parse_args()

    logging.basicConfig(level = args.loglevel)

    if 'SGE_TASK_ID' not in os.environ:
        PrepareJobs(joblist = args.taskfile, overwrite = args.overwrite)
    else: #In SGE array task
        SGEArrayTaskHandler(ge_task_id = os.environ['SGE_TASK_ID'], \
            taskinfofile = args.taskfile, overwrite = args.overwrite, qchemcmd = args.qcprog)


