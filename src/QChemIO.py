#!/usr/bin/env python
"""
Low-level routines to deal with Q-Chem input decks and
output files

Jiahao Chen, 2009-12
non-update 2011-03
"""

_qchemcmd = 'qchem'

####################
# Code starts here #
####################
import logging, os, scipy

logging.basicConfig()
logger = logging.getLogger(__name__)

class QChemInput:
    def __init__(self, InputDeckFilename = None):
        logger.debug('Initializing new QChemInput() with filename: %s',
                     InputDeckFilename)
        
        self.filename = InputDeckFilename
        self.jobs = []
        if InputDeckFilename is not None and \
                os.path.exists(InputDeckFilename):
            self.ImportQChemOutput(InputDeckFilename)


    def isValid(self):
        "Checks that data contained can be used to make a valid input deck"
        #Delete possibly empty first job
        if repr(self.jobs[0]) == '':
            self.jobs.pop(0)
                
        return all([job.isValid() for job in self.jobs])



    def GetCurrentJob(self):
        if len(self.jobs) == 0: #Create empty job
            self.jobs.append(QChemInputJob())

        return self.jobs[-1]



    def write(self, InputDeckFilename = None):
        """
        Writes input file

        InputDeckFilename: optional. Name of Q-Chem input file to write.
        If None, will return the entire output file as a string,
        i.e. same as calling __repr__()
        """

        if InputDeckFilename != None:
            self.filename = InputDeckFilename

        assert self.isValid(), "Insufficient data for valid Q-Chem input deck."

        buf = []
        for job in self.jobs:
            for input_section_keyword, input_section in job.sections:
                buf.append('$' + input_section_keyword + '\n')
                buf.append(input_section)
                buf.append('\n$end\n\n')
            buf.append('@@@\n')

        #Delete last end-of-job marker
        buf.pop(-1)

        outputtext = ''.join(buf)
        if self.filename == None:
            return outputtext
        else:
            with open(self.filename, 'w') as f:
                f.write(outputtext)


    def ImportQChemOutput(self, filename):
        job = QChemInputJob()
        job.filename = filename
        input_section_buffer = []
        input_section_keyword = None
        state = 'none'
        jobid = 0

        filenameroot = '.'.join(filename.split('.')[:-1])
        filenameext = filename.split('.')[-1]

        for line in open(filename):
            if '$end' in line:
                state = 'none'

                if input_section_keyword == 'molecule':
                    oldline1 = input_section_buffer[0]
                    #buf = SanitizeMolecule(filename)
                    #Q-Chem geometry optimization output is broken in 3.2:
                    #1. Final geometry is always output in some broken Z-matrix
                    #2. Number of electrons output is incorrect with ECP -__-
                    #manual override
                    input_section_buffer[0] = oldline1

                elif input_section_keyword == 'rem':
                    input_section_buffer = SanitizeRem(input_section_buffer)

                job.write_input_section(input_section_keyword, \
                     ''.join(input_section_buffer))

                input_section_buffer = []

            elif '$' in line:
                input_section_keyword = line.split()[0][1:]
                state = 'readblock'

            elif state == 'readblock':
                input_section_buffer.append(line)

            elif '@@@' in line:
                job.filename = '.'.join([filenameroot + '-job' + str(jobid),
                                         filenameext])
                jobid += 1
                self.jobs.append(job)
                job = QChemInputJob()

        job.filename = '.'.join([filenameroot + '-job' + str(jobid),
                                 filenameext])

        self.jobs.append(job)



    def execute(self, outfilename = None, qchemcmd = _qchemcmd,
                parse = False, savedir = None, Overwrite = False):
        """Runs Q-Chem input file. Returns a QChemOutput class."""

        #Commits current data to disk
        self.write()
        if outfilename is None:
            outfilename = self.filename[:-2] + 'out'

        if savedir is not None:
            #if 'QCLOCALSCR' in os.environ and 'QCSCRATCH' in os.environ \
            #        and os.environ['QCLOCALSCR'] == os.environ['QCSCRATCH']:
            #    del os.environ['QCLOCALSCR']

            dosave = '-save'
        else:
            dosave = ''
            savedir = ''


        OutputExists = os.path.exists(outfilename)
        if Overwrite or not OutputExists:
            from subprocess import Popen, STDOUT, PIPE
            try:
                cmdstring = ' '.join([qchemcmd, dosave, self.filename,
                    outfilename, savedir])
            except TypeError, e:
                logger.error("""Could not generate valid shell command line:
qchemcmd    = %s
filename    = %s
outfilename = %s
""", qchemcmd, self.filename, outfilename)
                raise e

            p = Popen(cmdstring, stdout = PIPE, stderr = STDOUT, shell = True)
            p.wait()
            print 'Console output is: %s' % p.stdout.read()
            logger.info('Console output is: %s', p.stdout.read())
        else:
            if OutputExists:
                logger.warning('Output file exists: %s and overwrite was not \
specified. Refusing to run Q-Chem.', outfilename)
        if parse:
            return QChemOutput(outfilename)



    def __repr__(self):
        return '\n@@@\n'.join([x.__repr__() for x in self.jobs])


    def __str__(self):
        return self.__repr__()


class QChemInputJob:

    #Here is a subclass for the Rem block
    """
    class QChemInputRemBlock(dict):
        def __init__(self, *args, **kwargs):

            if 'inputstring' in kwargs:
                inputstring = kwargs['inputstring']
                del kwargs['inputstring']
            else:
                inputstring = None

            dict.__init__(self, args, kwargs)
            
            if inputstring is None:
                return

            for line in inputstring.split('\n'):
                if line[0] == '!':
                    #Throw away comment
                    continue
                
                t = line.split()
                if len(t) >= 2:
                    keyword = t[0]
                    value = t[1]
                    dict.__setitem__(self, keyword, value)
                    #Throw away comments
    """
    def __init__(self, filename = None):
        self.sections = []
        self.filename = filename
        #self.rem = self.QChemInputRemBlock() #TODO finish migration to RemBlock

    def isValid(self):
        "Checks that data contained can be used to make a valid input deck"
        return self.has_input_section('rem') and \
               self.has_input_section('molecule')



    def has_input_section(self, input_section_keyword):
        return (input_section_keyword in [x for x, _ in self.sections])



    def read_input_section(self, input_section_keyword):
        for keyword, section in self.sections:
            if input_section_keyword == keyword:
                return section

        return None



    def delete_input_section(self, input_section_keyword):
        newsections = []
        for keyword, section in self.sections:
            if input_section_keyword != keyword:
                newsections.append([keyword, section])
        self.sections = newsections



    def write_input_section(self, input_section_keyword, input_section,
                            Overwrite = True):
        if Overwrite and self.has_input_section(input_section_keyword):
            self.delete_input_section(input_section_keyword)
        elif not Overwrite and self.has_input_section(input_section_keyword):
            return #Do nothing
        self.sections.append([input_section_keyword, input_section])

        #Special handling for REM block
        if input_section_keyword == 'rem':
            pass

    def append_input_section(self, input_section_keyword, input_section):

        #find requested input section keyword
        input_section_keywords = [x for x, _ in self.sections]
        n = -1
        for i, this_keyword in enumerate(input_section_keywords):
            if this_keyword == input_section_keyword:
                n = i
                break

        if n == -1: #If we got here, could not find input section
            logger.debug('Creating new input section: %s',
                         input_section_keyword)
            self.write_input_section(input_section_keyword, input_section)
        else:
            self.sections[n][1] += input_section

    def get_input_section(self, input_section_keyword):
        for a, b in self.sections:
            if a.lower().strip() == input_section_keyword.lower().strip():
                return b


    #Specific to REM section
    
    def rem_normalize(self):
        "Cleans up rem block"
        rem_section = self.get_input_section('rem')
        rem_dict = {}
        for l in rem_section.split('\n'):
            t = l.split()
            try:
                keyword, value, comments = t[0], t[1], ' '.join(t[2:]) 
                rem_dict[keyword] = value, comments
            except IndexError:
                pass

        maxlen = max([len(w) for w in rem_dict.keys()])
        maxlen2 = max([len(w) for w, _ in rem_dict.items()])

        new_rem_section = []

        for keyword, (value, comments) in rem_dict.items():
            padlen = maxlen - len(keyword) + 1
            padlen2 = maxlen2 - len(value) + 1
            new_rem_section.append(keyword + ' '*padlen + value + ' '*padlen2\
                                       + comments)

        self.write_input_section('rem', '\n'.join(new_rem_section), \
                              Overwrite = True)


    def rem_get(self, keyword):
        rem_section = self.get_input_section('rem')
        for l in rem_section.split('\n'):
            t = l.split()
            if len(t) >= 2:
                if keyword.upper() == t[0].upper():
                    return t[1]
        
        return None

    def rem_set(self, keyword, value):
        rem_section = self.get_input_section('rem')
        new_rem_section = []
        found = False
        for l in rem_section.split('\n'):
            t = l.split()
            if len(t) >= 2:
                thiskeyword, thisvalue = t[0], t[1]
                if len(t) > 2:
                    thiscomment = ' '.join(t[2:])
                else:
                    thiscomment = ''

                if keyword.upper() == thiskeyword.upper():
                    found = True
                    new_rem_section.append(' '.join((thiskeyword, value,
                                                    thiscomment)))
                else:
                    new_rem_section.append(' '.join((thiskeyword, thisvalue,
                                                    thiscomment)))
        
        if not found: #add it in
            new_rem_section.append(' '.join((keyword, value)))

        self.write_input_section('rem', '\n'.join(new_rem_section))


    def rem_delete(self, keyword):
        "Comments out keyword from rem block"
        rem_section = self.get_input_section('rem')
        new_rem_section = rem_section.replace(keyword, '!' + keyword)
        self.write_input_section('rem', new_rem_section, Overwrite = True)
        


    def execute(self, outfilename = None, qchemcmd = _qchemcmd, parse = True,
                savedir = '', Overwrite = False):
        """Runs Q-Chem input file. Returns a QChemOutput class."""

        if outfilename == None:
            outfilename = self.filename[:-2] + 'out'

        #Commits current data to disk
        tmpjob = QChemInput(self.filename)
        tmpjob.jobs.append(self)
        tmpjob.write()
        return tmpjob.execute(outfilename = outfilename, qchemcmd = qchemcmd,
                              Overwrite = Overwrite, savedir = savedir,
                              parse = parse)



    def write(self, filename = None, doCheck = True):
        """
        Writes input file
        
        filename name of filename to write.
        If None, returns as string. Same as __repr__()
        """

        if doCheck:
            assert self.isValid(), "Insufficient data for valid Q-Chem input deck."

        outputtext = self.__repr__()
        if filename == None:
            return outputtext
        else:
            with open(filename, 'w') as f:
                f.write(outputtext)

            self.filename = filename



    def __repr__(self):
        buf = []
        for input_section_keyword, input_section in self.sections:
            buf.append('$' + input_section_keyword + '\n')
            buf.append(input_section)
            buf.append('\n$end\n\n')

        return ''.join(buf)

    def __str__(self):
        return self.__repr__()



class QChemOutput:
    def __init__(self, filename = ''):
        #Check that file exists
        self.filename = filename

        if os.path.isfile(filename):
            #Count basis functions
            self.BasisCount = None
            self.CountBasisFunctions()
        else:
            error_msg = "Filename " + filename + " does not appear to be valid"
            logger.error(error_msg)
            raise IOError, error_msg


    def ReadMatrix(self, title, readrange = None):
        """Parses a Q-Chem output file for a matrix named _title_."""
        state = 'FindTitle'

        Matrix = []
        Buffer = []
        thisinstance = -1
        Matrices = []
        if readrange == None: readrange = [0]
        for line in open(self.filename):
            if state == 'FindTitle':
                if title in line:
                    thisinstance += 1
                    if thisinstance in readrange:
                        state = 'ReadCols'

            elif state == 'ReadCols':
                w = line.split()

                #Are we done? Check if first and last tokens in line are numbers
                try:
                    int(w[0])
                    float(w[-1])
                except (IndexError, ValueError):
                    #Done reading in stuff, this is other output

                    #Transfer columns that have already been read in
                    for MatrixCol in Buffer:
                        Matrix.append(MatrixCol)

                    #Done reading in Q-Chem output, now do sanity checking
                    Matrix = scipy.matrix(Matrix)
                    if Matrix.shape[1] == 0:
                        #Matrix not found
                        assert False, 'Error parsing instance '+\
                            str(thisinstance)+' of matrix '+title
                        Matrices.append(None)
                    else:
                        Matrices.append(Matrix)

                    if thisinstance == max(readrange): #Done
                        break
                    else: #Continue reading
                        state = 'FindTitle'

                #Check if it's a header or not
                #See if last entry contains a decimal point
                if '.' not in w[-1]: #This is a header

                    #Transfer columns that have already been read in
                    for MatrixCol in Buffer:
                        Matrix.append(MatrixCol)

                    #Read in new header
                    FirstColId = int(w[0])
                    LastColId = int(w[-1])

                    Buffer = []
                    for i in range(FirstColId, LastColId + 1):
                        Buffer.append([])

                else: #This is not a header

                    #Q-Chem's format is rownum __ __ ... __            

                    #Skip first entry
                    w.pop(0)

                    for i, entry in enumerate(w):
                        Buffer[i].append(float(entry))

        if len(Matrices) == 0:
            return
        elif len(Matrices) == 1:
            Matrices = Matrices[0]
        
        return Matrices


    def CheckFatalError(self):
        """Reads out part of Q-Chem output file with fatal error.
        Returns an exception that contains the error message and
        everything after it.
        If no error, returns empty string."""

        output = []

        accumulate = False
        for line in open(self.filename):
            if 'fatal error' in line:
                accumulate = True

            if accumulate:
                output.append(line)

        error = ''.join(output)
        if error != '':
            raise ValueError, "Error detected in Q-Chem output file " +\
                self.filename + '\n' + error



    def GetMolden(self):
        """Reads a Molden input file from the Q-Chem output file"""

        state = 0
        buf = []
        for line in open(self.filename):
            if '======= MOLDEN-FORMATTED INPUT FILE FOLLOWS =======' in line:
                state = 1
            elif '======= END OF MOLDEN-FORMATTED INPUT FILE =======' in line:
                break
            elif state == 1:
                buf.append(line)
        return ''.join(buf)



    def CountBasisFunctions(self, AtomList = None):
        """Counts the number of basis functions associated with atoms in
        AtomList as printed in the Q-Chem output file.

        If called with no AtomList, counts everything

        Populates self.BasisCount if not already populated."""

        if self.BasisCount == None:
            self.BasisCount = []

            dspher = True
            #True = use Spherical d Gaussians, else use Cartesian d Gaussians

            state = 'Find header'
            for line in open(self.filename):

                #Check for hard-coded use of Spherical vs
                #Cartesian Gaussians per basis set
                if 'Requested basis set' in line:
                    if '6-31' in line and '6-311' not in line:
                        logger.debug("Detected Pople double zeta basis set. \
Setting Cartesian d Gaussians")
                        dspher = False

                if state == 'Find header':
                    if "Atom   I     L     Exponents" in line:
                        state = 'Skip line'

                elif state == 'Skip line':
                    state = 'Count'
                    NumBasis = 0

                elif state == 'Count':
                    if '----------------' in line: #Read to the end
                        self.BasisCount.append(NumBasis)
                        self.BasisCount.pop(0)
                        break

                    atomid = line[:5].strip()
                    if atomid != '':
                        self.BasisCount.append(NumBasis)
                        NumBasis = 0

                    #The only other information we need is the angular momentum
                    #quantum number
                    L = line[12:19].strip()

                    if L == '0-1':
                        NumBasis += 4
                    elif L == '':
                        pass
                    else:
                        ll = int(L)
                        if ll == 2: #d Gaussians
                            if dspher: #Using spherical d Gaussians?
                                NumBasis += 5
                            else:
                                NumBasis += 6
                        elif ll < 5:
                            #Q-Chem uses spherical Gaussians for d,f,g
                            NumBasis += 2 * ll + 1
                        elif ll >= 5:
                            #Q-Chem uses Cartesian Gaussians for h
                            NumBasis += (ll + 1) * (ll + 2) / 2

                if 'PURECART' in line.upper():
                    error_msg = """\
Hey! You specified something for PURECART!"
I don't know what to do!"""
                    logger.error(error_msg)
                    raise ValueError, error_msg

        #Done parsing, now calculate
        if AtomList == None:
            Total = sum(self.BasisCount)
        else:
            Total = sum([self.BasisCount[x - 1] for x in AtomList])
        return Total



    def ReadDensityMatrix(self):
        "Reads a density matrix out of the output"

        ## For testing purposes, get initial guess matrices
        ## Remember to put SCF_GUESS_PRINT 2 in Q-Chem $REM block
        
        PA = self.ReadMatrix('Final Alpha density matrix')
        if PA is None:
            error_msg = "I don't know where the alpha density matrix is."
            logger.error(error_msg)
            raise ValueError, error_msg

        PB = self.ReadMatrix('Final Beta density matrix')
        if PB is None: #Can't find beta density matrix
            logger.debug("I don't know where the beta density matrix is; \
assuming it's the same as alpha.")
            PB = PA

        return PA + PB


###

def SanitizeRem(buf):
    newbuf = []
    keyword = ''
    for l in buf:
        if l[0] == '!':
            continue
        t = l.split()
        if len(t) >= 2:
            keyword = t[0]
            if keyword != 'jobtype':
                newbuf.append(l)

    return newbuf


def SanitizeMolecule(filename):
    """Sometimes Q-Chem fails to write correct optimized geometry block
    For now the sticky electron code only works with Cartesian formatted
    geometry
    Super cheat - use babel to generate correct geometry block"""

    os.system('babel -iqcout ' + filename + ' -oqcin TMP')
    buf = []
    state = None
    for l in open('TMP'):
        if '$end' in l:
            state = None
        elif '$molecule' in l:
            state = 'add'
        elif state == 'add':
            buf.append(l)
    return buf
