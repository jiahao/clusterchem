#!/usr/bin/env python

"""
Generates QM/MM electronic structure jobs and submits them to a Sun Grid
Engine queue.
"""
import logging, os

from glob import glob
from subprocess import Popen, PIPE, STDOUT

from OSUtils import chdirn, linkwild
from QChemInterface import QChemInputForElectrostaticEmbedding, \
         QChemInputForTransitionDipole
from QChemIO import QChemInput

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

    logger = logging.getLogger('domyjob')
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

    #Link coordinate file
    try:
        os.symlink(os.path.join('..', '..', '..', corfileroot+'.cor'),
            corfileroot.upper()+'.COR')
    except OSError, e:
        logger.warning('Python reported error making symlink to %s: %s. Ignoring'
                , corfileroot.upper()+'.COR', e)

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
        print 'Running nonpolarizable QM/MM TDDFT for \
                resid =', resid
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

            Q1 = QChemInput('QCHEM-GROUND')
            Q1 = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                    ['h2pc.rtf', 'znpc.rtf'], Q1)

            Q2 = QChemInput('QCHEM-DSCF')
            Q2 = QChemInputForElectrostaticEmbedding(resid, charmmcardfile,
                    ['h2pc.rtf', 'znpc.rtf'], Q1)

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
        logging.error('Need to run this as a Grid Engine array job')
    else: #In SGE array task
        SGEArrayTaskHandler(ge_task_id = os.environ['SGE_TASK_ID'],
            taskinfofile = args.taskfile, overwrite = args.overwrite,
            qchemcmd = args.qcprog)


