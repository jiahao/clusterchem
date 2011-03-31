#!/usr/bin/env python

Hartree_to_eV = 27.2113839
Hartree_to_kcal_mol = 627.5095
kcal_mol_to_ev = 23.060549

EnergyFlag = '\033[1;33m!Energy\033[1;m'
GRMSFlag = '\033[1;33m!GRMS\033[1;m'
DoneFlag = '\033[1;30mDone\033[1;m'
RunningFlag = '\033[1;30mRunning\033[1;m'
AbortFlag = '\033[1;41mABORTED\033[1;m'

SmallAbortFlag = '*'
SmallRunningFlag = '.'
NotAvailableFlag = '   N/A    '

_qstat_bin = '/opt/gridengine/bin/lx24-amd64/qstat'

__doc__ = """Polarizable QM/MM (QChem/CHARMM) Jobs monitor
Jiahao Chen <jiahao@mit.edu> 2011-02-24

Key:
"""+EnergyFlag +"""\tEnergy is not approximately monotonically convergent
"""+GRMSFlag   +"""\tRMS Gradient has pronounced spike
"""+DoneFlag   +"""\tJob has completed successfully
"""+RunningFlag+"""\tJob is still running
"""+AbortFlag  +"""\tJob aborted
"""

# Here we process arguments

import argparse

parser = argparse.ArgumentParser(description='Monitors QChem/CHARMMM QM/MM jobs monitor')
parser.add_argument('workingdir', help='Working directory', nargs = '?', default = '.')
parser.add_argument('--kill', action='store_true', default = False, help='Force termination of very slowly converging jobs')
parser.add_argument('--nuke', action='store_true', default = False, help='Automatically rename output files of aborted jobs')
parser.add_argument('--ehist', action='store_true', default = False, help='Print convergence history of QM/MM energy')
parser.add_argument('--summary-only', action='store_true', default = False, help='Print summary only')
args = parser.parse_args()
###

import fnmatch
import glob
import os
from subprocess import Popen, PIPE
from heapq import heapify, heappop, heappush

def find_files(directory, pattern):
    """Recursively traverses directory tree IN ORDER for files matching pattern
    Note: does NOT recurse into symbolic links.
    """


    stack = [directory] #A priority queue
    while stack:
        nowdir = heappop(stack)
        for base in os.listdir(nowdir):
            name = os.path.normpath(os.path.join(nowdir, base))
            if os.path.isdir(name):
                if not os.path.islink(name): #Avoid infinite loop with symlinks
                    heappush(stack, name)
                    #Use default lexicographical order for priority (sort)
            elif fnmatch.fnmatch(name, pattern):
                yield nowdir, name



def ParseOutStatus(filename):
    'Checks status of QM/MM job based on output file'

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
                e = float(l.split()[1]) * Hartree_to_kcal_mol
                e_history.append(e)

            if 'PURIFY final SCF energy' in l:
                e = float(l.split()[-1]) * Hartree_to_kcal_mol
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
          if  grms_history[i] > grms_history[i+1] + grms_history[i-1]:
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



class SGEJobData:
    def __init__(self, job_number = None, qstat = _qstat_bin):
        if job_number == None:
            self.job_number = None
        else:
            self.job_number = int(job_number)
        self.qstatbin = qstat
        #if job_number != None:
        #    self.CallQstatJ()



    def CallQstatJ(self):
        """Populates data from SGE using qstat -j

        Note. This will NOT work as expected for scheduling info and parallel environment
              but I don't expect to use these data, so I don't care.

        """

        P = Popen([self.qstatbin, '-j', str(self.job_number)], stdout = PIPE)
        for l in P.stdout:
            t = l.split()
            if t[0][-1] == ':':
                attributename = t[0][:-1]
                attributevalue = t[1]

                #Convert data type
                if attributename in ['job_number']:
                    attributevalue = int(attributevalue)
             
                if attributename == 'job_number' and self.job_number != None:
                    assert attributevalue == self.job_number
                setattr(self, attributename, attributevalue)



    def CallQstat(self):
        """Populates data from SGE using qstat

        Not all data reported by qstat can be obtained from qstat -j

        Note. If instantiating many SGEJobData, it is preferable to poll qstat directly.
        and use ParseQstatLine(). See SGEGetRunningList()
        """
        P = Popen(self.qstatbin, stdout = PIPE)
        for l in P.stdout:
            self.ParseQstatLine(l)

    def ParseQstatLine(self, l):
        """Populates data from SGE using a line from qstat

        Not all data reported by qstat can be obtained from qstat -j
        """
        t = l.split()
        jobid = int(t[0])
        if self.job_number == None:
            self.job_number = jobid
        else:
            assert jobid == self.job_number
        self.prior = t[1]
        self.name = t[2]
        self.user = t[3]
        self.state = t[4]
        self.submit = t[5] + ' ' + t[6]
        if 'r' in self.state:
            self.queue, self.host = t[7].split('@')

        self.slots = t[-1]



def SGEGetRunningList(qstat = _qstat_bin):
    """Returns a list of running jobs from SGE using qstat"""

    P = Popen(qstat, stdout = PIPE)
    jobs = []
    for l in P.stdout:
        try:
            job_number = int(l.split()[0])
            thisjob = SGEJobData(job_number, qstat)
            thisjob.ParseQstatLine(l)
            if 'r' in thisjob.state: thisjob.CallQstatJ()
            jobs.append(thisjob)
        except (ValueError, IndexError):
            pass

    return jobs



def main(path = '.', extension = '.out'):
    jobs = SGEGetRunningList()

    Aborted = []
    Data = {}
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
    heapify(pq)
    print '\n\nSummary of excitation energies computed:'
    numthings = 0
    while True:
        try:
            jobroot, jobdata = heappop(pq)
            numthings += 1
        except IndexError:
            break
        print jobroot, 
        try: 
            print 'non-pol:',
            print ' %.3f eV' % ((jobdata['dscf-nonpol'][0] - jobdata['gs-nonpol'][0])/kcal_mol_to_ev),
            if jobdata['gs-nonpol'][2] or jobdata['dscf-nonpol'][2]: print SmallAbortFlag,
            elif jobdata['gs-nonpol'][1] or jobdata['dscf-nonpol'][1]: print SmallRunningFlag,
            else: print ' ',
        except KeyError:
            print NotAvailableFlag,
        try: 
            print 'pol:',
            print '%.3f eV' % ((jobdata['dscf'][0] - jobdata['gs'][0])/kcal_mol_to_ev),
            if jobdata['gs'][2] or jobdata['dscf'][2]: print SmallAbortFlag,
            elif jobdata['gs'][1] or jobdata['dscf'][1]: print SmallRunningFlag,
            else: print ' ',
        except KeyError:
            print NotAvailableFlag,
        print#if numthings % 3 == 0: print

if __name__ == '__main__':
    if not args.summary_only: print __doc__
    print 'Running in directory',args.workingdir
    main(args.workingdir)
        
