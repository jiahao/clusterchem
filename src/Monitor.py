#!/usr/bin/env python

#Internal
from OSUtils import find_files
from Units import eV, kcal_mol

EnergyFlag = '\033[1;33m!Energy\033[1;m'
GRMSFlag = '\033[1;33m!GRMS\033[1;m'
DoneFlag = '\033[1;30mDone\033[1;m'
RunningFlag = '\033[1;30mRunning\033[1;m'
AbortFlag = '\033[1;41mABORTED\033[1;m'

SmallAbortFlag = '*'
SmallRunningFlag = '.'
NotAvailableFlag = '   N/A    '

_qstat_bin = 'qstat'

__doc__ = """Polarizable QM/MM (QChem/CHARMM) Jobs monitor
Jiahao Chen <jiahao@mit.edu> 2011-02-24

Key:
"""+EnergyFlag +"""\tEnergy is not approximately monotonically convergent
"""+GRMSFlag   +"""\tRMS Gradient has pronounced spike
"""+DoneFlag   +"""\tJob has completed successfully
"""+RunningFlag+"""\tJob is still running
"""+AbortFlag  +"""\tJob aborted
"""

###

import glob, os



def ParseOutStatus(filename):
    """
    Checks status of QM/MM job based on output file

    :param string filename:
        :term:`CHARMM` or :term:`Q-Chem` output file to monitor

    :returns: the same data as :py:func:`ParseOutput`.

    .. TODO:: There is some overlap with :py:func:`ParseOutput`; consolidate.

    .. versionadded:: 0.1
    """

    grms_history = []
    e_history = []
    lastmm_e_history = []
    isdone = False
    charmmiter = 0
    qm_iterthresh = 0
    tdip = None
    mode = 'default'

    for l in open(filename):
        if mode == 'default':
            #Parsing for CHARMM output
            if l[:5] == 'SCF >' or l[:5] == 'ENER>':
                charmmiter = int(l[5:14])
                e = float(l[14:26])
                grms = float(l[40:53])
                #e, _, grms = [float(x) for x in l[14:].split()]
                e_history.append(e)
                lastmm_e_history.append(e)
                grms_history.append(grms)

            if 'SCF_CONVERGENCE' in l:
                lastmm_e_history = []
                qm_iterthresh = int(l.split()[-1])            

            if 'The kill switch has been activated' in l:
                isdone = True
                break

            if 'NORMAL TERMINATION BY END OF FILE' in l:
                isdone = True
                break


            #Parsing for Q-Chem output
            if 'Convergence criterion met' in l:
                e = float(l.split()[1]) / kcal_mol
                e_history.append(e)

            if 'PURIFY final SCF energy' in l:
                e = float(l.split()[-1]) / kcal_mol
                e_history.append(e)

            if '*** MISSION COMPLETED -- STARFLEET OUT ***' in l:
                isdone = True
                break

            #Read transition dipole moment
            if 'dipole x component in diabatic basis' in l:
                tdip = None
                import numpy, QChemIO
                q = QChemIO.QChemOutput(filename)
                x = q.ReadMatrix('dipole x component in diabatic basis')
                y = q.ReadMatrix('dipole y component in diabatic basis')
                z = q.ReadMatrix('dipole z component in diabatic basis')
                tdip = numpy.asarray(x, y, z)
                print tdip; exit()
                break
            
        """
        if mode == 'ReadTdip':
            t = l.split()
            if t[0] == 'dipole' and 'component in diabatic basis' in l:
                direction = t[1] - ord('x') #x, y, z --> 0, 1, 2

            #Test integer
            if str(int(t[-1])) == t[-1]:
                numstates = int(t[-1])
                if tdip == None:
                    tdip = zeros((numstates, numstates, 3))
            else:
                thisstate = int(t[0])
                data = [float(x) for x in t[1:]]
                tdip[thisstate, :, direction] = data
        """

    #Check for large fluctuations in rms gradient
    #The heuristic here is to detect spikes of at least twice the
    #height of the average of the grms immediately before and after
    isgrmsweird = False
    for i in range(1,charmmiter-2):
        try:
            if grms_history[i] > grms_history[i+1] + grms_history[i-1]:
                isgrmsweird = True
                break
        except IndexError:
            pass

    #Check for nonmonotonic convergence in energy
    isenergyweird = False
    tol = 1e-2
    for i in range(len(lastmm_e_history)-2):
        if e_history[i] + tol < e_history[i+1]:
            isenergyweird = True
            break
    
    return charmmiter, e_history, tdip, qm_iterthresh, isenergyweird, isgrmsweird, isdone



def main(path = '.', extension = '.out'):
    """
    main

    .. TODO:: DOCUMENT

    .. versionadded:: 0.1
    """
    from SGE import SGE
    sge_jobs = SGE()
    jobs = sge_jobs.getuserjobs()

    Aborted = []
    Data = {}
    tdip = None
    for root, filename in find_files(path, '*'+extension):
        jobname = filename#[:-len(extension)]
        charmmiter, e_history, tdip, qm_iterthresh, isenergyweird, isgrmsweird, isdone \
            = ParseOutStatus(filename)
        if not args.summary_only:
            print jobname, 
            if isenergyweird: print EnergyFlag,
            if isgrmsweird:   print GRMSFlag,
            if isdone:        print DoneFlag,
        thisjob = None
        isAborted = False
        if not isdone: #Not done, check if it's still in the queue
            #The best we can do is match directories
            thispath = os.path.realpath(root)
            for job in jobs:
                if 'r' in job.state and job.cwd == thispath:
                    thisjob = job
                    break
            if thisjob != None:
                if not args.summary_only: print RunningFlag, '(SGE %d)' % thisjob.job_number,
            else:
                isAborted = True
                if not args.summary_only: print AbortFlag,
                Aborted.append(root)
        
        #Print some info about job
        if len(e_history):
            energy = e_history[-1]
            if not args.summary_only:
                print 'Energy = %.4f kcal/mol' % energy,
        else:
            energy = 0.

        if not args.summary_only and not isdone:
            if charmmiter > 1:
                print '(%+.4f) at iteration %d, SCF threshold = 1e-%d' % \
                (e_history[-1]-e_history[-2], charmmiter, qm_iterthresh),
            else:
                print 'iteration = %d' % charmmiter,
 
        if args.kill and thisjob != None and qm_iterthresh == 8 and abs(e_history[-1]-e_history[-2]) < 5e-2:
            #Heuristic for being done - energy converged to 0.05 kcal/mol
            #ssh $hostname 'echo 1 > /scratch/cjh/$jobid.qmmmpol/killswitch'
            os.system("ssh %s \'echo 1 > /scratch/$USER/%d.qmmmpol/killswitch\'" % (thisjob.host, thisjob.job_number))
            print '\n"Dr. Strangelove, the kill switch has been deployed."'

        if args.ehist and len(e_history):
            print '\nEnergy history'
            for n, e in enumerate(e_history):
                print 'Iteration = %4d Energy = %.4f kcal/mol' % (n,e)

        if not args.summary_only: print

        #Save data for calculations
        trunc = 0


        if root.endswith('/gs'): trunc = 3
        if root.endswith('/gs-nonpol'): trunc = 10
        if root.endswith('/dscf'): trunc = 5
        if root.endswith('/dscf-nonpol'): trunc = 12
        if trunc > 0:
            jobroot = root[:-trunc]
            jobstem = root[1-trunc:]
            if jobroot not in Data: Data[jobroot] = {}
            if energy != 0.0: Data[jobroot][jobstem] = energy, not isdone, isAborted

        #This fixes a bug in Q-Chem wrapper script that forgot to parse nonpolarizable
        #delta-SCF correctly
        if 'qchem.out' in filename and 'dscf' in root:
            if len(e_history) >= 2:
                e_correction = e_history[-1] - e_history[-2]
            root = os.path.dirname(root)
            if root.endswith('/dscf'): trunc = 5
            if root.endswith('/dscf-nonpol'): trunc = 12
            jobroot = root[:-trunc]
            jobstem = root[1-trunc:]
            try:
                Data[jobroot][jobstem] = ((Data[jobroot][jobstem][0] + e_correction), not isdone, isAborted)
            except (KeyError, UnboundLocalError):
                pass

    if tdip != None:
        print 'Hey, a transition dipole moment!'
        print tdip
        exit()

    if len(Aborted) > 0:
        print 'List of aborted jobs:', ' '.join(Aborted)
        if args.nuke:
            for job in Aborted:
                for filename in glob.glob(os.path.join(job,'*.out')):
                    ext = 0
                    while True:
                        newfilename = filename+'.'+str(ext)
                        if not os.path.exists(newfilename):
                            os.rename(filename, newfilename) 
                            print filename, '-->', newfilename
                            break
                        else:
                            if ext == '': ext = 0
                            else: ext += 1

    #Calculate random stuff
    pq = Data.items()
    from heapq import heapify, heappop

    heapify(pq)
    print '\n\nSummary of excitation energies computed:'
    
    while True:
        try:
            jobroot, jobdata = heappop(pq)
        except IndexError:
            break
        print jobroot, 
        try: 
            print 'non-pol:',
            print ' %.3f eV' % ((jobdata['dscf-nonpol'][0] - jobdata['gs-nonpol'][0])*kcal_mol/eV),
            if jobdata['gs-nonpol'][2] or jobdata['dscf-nonpol'][2]: print SmallAbortFlag,
            elif jobdata['gs-nonpol'][1] or jobdata['dscf-nonpol'][1]: print SmallRunningFlag,
            else: print ' ',
        except KeyError:
            print NotAvailableFlag,
        try: 
            print 'pol:',
            print '%.3f eV' % ((jobdata['dscf'][0] - jobdata['gs'][0])*kcal_mol/eV),
            if jobdata['gs'][2] or jobdata['dscf'][2]: print SmallAbortFlag,
            elif jobdata['gs'][1] or jobdata['dscf'][1]: print SmallRunningFlag,
            else: print ' ',
        except KeyError:
            print NotAvailableFlag,
        print



if __name__ == '__main__':
    # Process arguments

    import argparse
    parser = argparse.ArgumentParser(description='Monitors QChem/CHARMM QM/MM jobs monitor')
    parser.add_argument('workingdir', help='Working directory', nargs = '?', default = '.')
    parser.add_argument('--kill', action='store_true', default = False, help='Force termination of very slowly converging jobs')
    parser.add_argument('--nuke', action='store_true', default = False, help='Automatically rename output files of aborted jobs')
    parser.add_argument('--ehist', action='store_true', default = False, help='Print convergence history of QM/MM energy')
    parser.add_argument('--summary-only', action='store_true', default = False, help='Print summary only')
    args = parser.parse_args()

    if not args.summary_only: print __doc__
    print 'Running in directory',args.workingdir
    main(args.workingdir)
        
