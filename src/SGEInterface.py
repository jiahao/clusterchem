#!/usr/bin/env python

"""
An experimental script

Compute transition dipole moment

Jiahao Chen <jiahao@mit.edu> 2011-02-18

TODO use SGE interface
"""

#Number of processors
np=1

#QCPROG for transition dipole calculation
_qchem_td_qcprog = '/home/tkowalcz/qchem/qchem-2ecoupling/linux64/exe/qcprog.exe'

try:
    from CHARMMUtil import dedrude
except ImportError: pass

from QChemInterface import QChemInputForTransitionDipole
from glob import glob
import os

def main():
    alljobs = []

    for coordfile in glob('*.d.cor'):
        coordfileroot = coordfile[:-6]

        #Generates coordinate file without Drude particles
        nonpolcoordfile = coordfileroot +'.cor'
        if not os.path.exists(nonpolcoordfile):
            dedrude(coordfile, nonpolcoordfile)

        if True or not os.path.exists(coordfileroot):
            print 'Preparing', coordfileroot
            if not os.path.exists(coordfileroot): os.mkdir(coordfileroot)
        
            os.chdir(coordfileroot)
            #Copy runpack
            copywild('/home/cjh/rmt/runpack/*', '.')

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

        alljobs += iteratejobs(coordfileroot)

        #exit()

    f = open('sgearraytasklist.txt', 'w')
    for i, j in enumerate(alljobs):
        f.write(str(i+1)+' \t')
        f.write(' '.join(j) + '\n')
    f.close()

    numjobs = len(alljobs)
    print 'Indexed', numjobs, 'jobs'

    #Now create SGE array job and submit self as an array job
    import SGE
    import sys
    x = SGE.qsubOptions()
    x.args.N = 'h2pc-qmmm-rmt'
    x.args.t = '1-'+str(numjobs)
    x.args.o = '$JOB_NAME.$JOB_ID.$TASK_ID.stdout.log'.replace('$', '\$')
    x.args.command = sys.argv[0]
    print '#To submit the array job, run this command:'
    x.execute(mode = 'echo')



def copywild(src, dest, overwrite=False, Verbose = True):
    for fname in glob(src):
        destfname = os.path.join(dest, os.path.basename(fname))
        if overwrite:
            doit = True
        else:
            doit = not os.path.exists(destfname)
        if (os.path.isfile(fname) or os.path.islink(fname)) and doit:
            #if not os.path.exists(destfname) or (overwrite and os.path.exists(destfname)):
            if overwrite and os.path.exists(destfname): os.unlink(destfname)
            if not os.path.exists(destfname): os.symlink(fname, destfname)
            if Verbose: print os.path.abspath(fname), '-->', os.path.abspath(destfname)



def iteratejobs(corfileroot):
    'populates jobdata'
    resids = []
    for l in open('/home/cjh/rmt/qmmm-test/trdip/mol.list'):
        resid = int(l)
        resids.append(resid)

    jobdata = []

    #Initialize calculations
    for resid in resids:
        #for job in ['gs-nonpol', 'dscf-nonpol', 'gs', 'dscf']:
        #for job in ['td', 'td-nonpol']:
        for job in ['dscf-nonpol']:
            #If no output files exist
            if len(glob(os.path.join(corfileroot, str(resid), job, '*.out'))) == 0:
                #copywild(os.path.join(corfileroot,'*'), os.path.join(corfileroot, str(resid), job))
                jobdata.append((corfileroot, str(resid), job))
            else:
                print 'Exists:', ' '.join(glob(os.path.join(corfileroot, str(resid), job, '*.out')))
                return [] #Abort if somewhere, an output file exists.

    return jobdata



def RunQMMM(infile, corfile, qchemcnt, res1):
    from subprocess import Popen, PIPE
    outfile = infile[:-2]+'out'
    if not os.path.exists(outfile):
        os.environ['NP'] = str(np)
        p = Popen('charmmpol '+infile+' CORFILE='+ corfile +' QCHEMCNT='+qchemcnt+' np='+str(np)+ ' res1='+str(res1), stdout=PIPE, shell=True)
        print p.stdout.read()
    else:
        print 'Output file already exists:', os.path.abspath(outfile)
        print 'Skipping calculation'



def chdirn(thedir):
    'Changes into directory. If directory does not exist, make it.'
    if not os.path.exists(thedir):
        os.mkdir(thedir)
    os.chdir(thedir)



def domyjob(corfileroot, resid, jobtype):
    olddir = os.getcwd()

    #mkdir -p $corfileroot/$resid/$jobtype && cd [...]
    chdirn(corfileroot)
    copywild('/home/cjh/rmt/runpack/*', '.')

    chdirn(str(resid))
    chdirn(jobtype)

    #Initialize
    print os.getcwd()
    copywild(os.path.join('..','..','*'), '.')

    if jobtype == 'gs-nonpol':
        RunQMMM('h2pc.in', corfileroot+'.cor', 'QCHEM-GROUND.IN', resid)
    elif jobtype == 'dscf-nonpol':
        RunQMMM('h2pc.in', corfileroot+'.cor', 'QCHEM-DSCF.IN', resid)
    elif jobtype == 'gs':
        RunQMMM('h2pc-pol.in', corfileroot+'.d.cor', 'QCHEM-GROUND.IN', resid)
    elif jobtype == 'dscf':
        RunQMMM('h2pc-pol.in', corfileroot+'.d.cor', 'QCHEM-DSCF.IN', resid)
    elif jobtype == 'td':
        os.environ['QCPROG'] = _qchem_td_qcprog

        filenames = glob('../gs/*/qchem*working.in')\
                  + glob('../dscf/*/qchem*working.in')
        if len(filenames) == 2:
            print 'Running polarizable QM/MM transition dipole moment for resid =', resid
        Q = QChemInputForTransitionDipole(filenames[0], filenames[1])
        Q.filename = 'td.in'
        Q.jobs[1].rem_delete('purify')
        Q.execute(parse = False)
        if len(filenames) > 2:
            print 'Too many possible inputs:', ' '.join(filenames)

    elif jobtype == 'td-nonpol':
        os.environ['QCPROG'] = _qchem_td_qcprog

        filenames = glob('../gs-nonpol/*/qchem*working.in') \
                  + glob('../dscf-nonpol/*/qchem*working.in')
        if len(filenames) == 2:
            print 'Running nonpolarizable QM/MM transition dipole moment for resid =', resid
            Q = QChemInputForTransitionDipole(filenames[0], filenames[1])
            Q.filename = 'td.in'
            Q.jobs[1].rem_delete('purify')
            Q.execute(parse = False)
            if len(filenames) > 2:
                print 'Too many possible inputs:', ' '.join(filenames)

    else:
        assert False, 'Unknown jobtype '+str(jobtype)

    #Return to old directory
    os.chdir(olddir)


            
if __name__ == '__main__':
    if 'SGE_TASK_ID' not in os.environ:
        main()
    else: #In SGE array task 
        ge_task_id = os.environ['SGE_TASK_ID']
        #Look up data about job
        for line in open('sgearraytasklist.txt'):
            t = line.split()
            if len(t) >= 4 and t[0] == ge_task_id:
                thiscorfileroot, thisresid, thisjobtype = t[1], t[2], t[3]
                break
        #print l,
        domyjob(thiscorfileroot, thisresid, thisjobtype)



