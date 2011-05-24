#!/usr/bin/env python
"""
Archives data into a :term:`HDF5` database.

This is a standalone executable requiring :term:`PyTables`
and the homebrew :py:mod:`QChemIO`.

.. versionadded:: 0.1

"""

#TODO SHOULD WE ARCHIVE ELECTRONIC STRUCTURE CONVERGENCE INFORMATION?

from glob import glob
import argparse, numpy, os, os.path, shutil, sys, warnings


try:
    import tables
except ImportError:
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()

#Turn off NaturalNameWarning from PyTables
warnings.filterwarnings('ignore', category=tables.NaturalNameWarning)

#My custom module
try:
    import QChemIO
except ImportError:
    sys.stderr.write("""Could not find QChemIO module.
Please check that PYTHONPATH is set correctly.
""")
    exit()

#Internal modules
from ParseOutput import ParseOutput
from Units import kcal_mol


class CHARMM_CARD(tables.IsDescription):
    """
    :term:`PyTables` data structure for storing a :term:`CHARMM card file` in a :term:`HDF5` database.

    .. versionadded:: 0.1
    """
    #: Atomic number
    AtomNo = tables.UInt64Col()
    #: Residue number
    ResidNo = tables.UInt64Col()
    #: Residue label of up to 4 characters
    Res = tables.StringCol(4)
    #: Type
    Type = tables.StringCol(4)
    #: Coordinates
    Coord = tables.Float32Col(shape = 3,)
    #: SegID
    SegID = tables.StringCol(4)
    #: ResID
    ResID = tables.StringCol(4)
    #: Weighting
    Weighting = tables.Float32Col()

    def __init__(self):
        tables.IsDescription.__init__(self)



def LoadCHARMM_CARD(h5table, CARDfile):
    """
    Reads a :term:`CHARMM card file` and saves it into a :term:`HDF5` database.
    
    :argument h5table: Name of :term:`HDF5` table
    :type h5table: :term:`PyTables` :py:class:`tables.File` object
    
    :argument string CARDfile: Name of :term:`CHARMM card file`

    :returns: None

    .. versionadded:: 0.1
    """
    state = 'title'
    for l in open(CARDfile):
        if state == 'title':
            if not (len(l) > 0 and l[0] == '*'):
                state = 'numatoms'
        if state == 'numatoms':
            numatoms = int(l.split()[0])
            thisnumatoms = 0
            state = 'read'
        elif state == 'read':
            t = l.split()
            data = h5table.row
            data['AtomNo'] = int(t[0])
            data['ResidNo'] = int(t[1])
            data['Res'] = t[2]
            data['Type'] = t[3]
            data['Coord'] = float(t[4]), float(t[5]), float(t[6])
            data['SegID'] = t[7]
            data['ResID'] = t[8]
            data['Weighting'] = float(t[9])
            data.append()
            thisnumatoms += 1

    assert thisnumatoms == numatoms, 'Did not read in expected number of atoms.'
    print 'Loaded CHARMM file', CARDfile, 'into', h5table
    h5table.flush()



def CheckIfRunning(filename):
    """
    Is the process responsible for producing filename still running?

    Need to infer this from qstat

    Get list of running SGE processes
    then check working directories?
    """
    #TODO
    return False


def RunQChemAgain(filename):
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
    """


def LoadEmUp(path = '.', h5filename = 'h2pc-data.h5', Nuke = False, NukeInc = False, DoCheckPoint = False, DoOverwrite = False):
    """
    Recursively reads data from all files in a given path and loads it into a
    :term:`HDF5` database.
    
    :argument string path: Path to read data from recursively. (Default: '.')
    
    :argument string h5filename: Name of :term:`HDF5` database to create. \
        or update (Default: 'h2pc-data.h5')

    :argument Boolean Nuke: Delete files after successful archival (Default: False)
    :argument Boolean NukeInc: Delete incomplete output files detected (Default: False)
    :argument Boolean DoCheckPoint: Create new check point after successful
         archival (Default: False)

    :argument Boolean DoOverwrite: Force overwrite of existing data with data in existing files (Default: False)

    :returns: None

    .. TODO :: replace with a real way to test if node exists in Pytables
    .. TODO :: Failed to detect symlink

    .. versionadded:: 0.1
    """
    compress = tables.Filters(complevel = 9, complib = 'zlib')
    h5data = tables.openFile(h5filename, mode = 'a', filters = compress,
        title = 'H2Pc')
    if not h5data.isUndoEnabled(): h5data.enableUndo(filters = compress)
    for root, _, files in os.walk(path):
        for cwdfile in files:

            filename = os.path.join(root, cwdfile) #Relative path-qualified filename

            #Skip symlinks
            if os.path.islink(filename): continue

            ##############################
            # Parse CHARMM coordinate file
            ##############################

            if len(cwdfile) > 4 \
                and cwdfile[-4:] == '.cor' \
                and cwdfile[-6:-4] != '.d':

                # CARDs with Drude particles are skipped since they
                # will need to be re-equilibrated as part of the SCF process

                cardtitle = 'CHARMM CARD Coordinates'
                cardname = 'CHARMM_CARD'
                
                cwdh5name = cwdfile[:-4]#.replace('-', '_')

                try:
                    thischarmmcard = h5data.createTable(
                        where = os.path.join('/' + root, cwdh5name), \
                        name = cardname, \
                        description = CHARMM_CARD, \
                        title = cardtitle, createparents = True)
                    LoadCHARMM_CARD(thischarmmcard, cwdfile)
                except tables.exceptions.NodeError:
                    #print cwdfile, 'already exists; skipping.'
                    continue #Force continue
                    #TODO replace continue with resetting of table
                    if DoOverwrite:
                        thischarmmcard = h5data.getNode(
                        where = os.path.join('/' + root, cwdh5name), \
                        name = cardname)
                        LoadCHARMM_CARD(thischarmmcard, cwdfile)
                    else:
                        pass #print cwdfile, 'already exists; skipping.'
                


            ####################################
            # Parse CHARMM or Q-Chem output file
            ####################################

            if len(cwdfile) > 4 \
                and cwdfile[-4:] == '.out':

                #Assume path has the structure [.../]GeomName/Site/jobtype
                Site = root.split(os.sep)[-2]
                GeometryName = root.split(os.sep)[-3]
                calctype = os.path.split(root)[1]

                State = Energy = Energies = Dipole = Environment = None
                data = ParseOutput(filename)
                isDone = data[-1]
                EnergyList = data[1]
                Dipole = data[2]

                if not isDone:
                    #TODO Check if job is still running on SGE
                    isRunning = CheckIfRunning(filename)
                    if isRunning:
                        #If it is running, do nothing
                        print 'Skipping incomplete calculation in', filename
                        continue
                    
                    #If it isn't, check for errors
                    
                    #If there are convergence errors:
                    RunQChemAgain(filename)

                    if NukeInc:
                        print 'Deleting incomplete calculation in', root
                        shutil.rmtree(root)
                        continue

                if '-nonpol' in calctype:
                    Environment = 'Fixed'
                else:
                    Environment = 'WithDrude'
                    #TODO Record Drude particle positions

                if 'gs' in calctype:
                    State = 0
                    if len(EnergyList) > 0:
                        Energy = EnergyList[-1]
                    else: continue
                
                elif 'dscf' in calctype:
                    State = 1
                    if len(EnergyList) > 0:
                        Energy = EnergyList[-1]
                    else: continue
                    # DEPRECATED
                    # This fixed a bug in Q-Chem wrapper script that
                    # forgot to parse nonpolarizable
                    # delta-SCF correctly
                    # The script itself has been corrected and this
                    # should no longer be necessary.
                    """
                    if calctype == 'dscf-nonpol':
                        Environment = 'Fixed'
                        State = 1
                        Site = root.split(os.sep)[-2]
                        GeometryName = root.split(os.sep)[-3]#.replace('-', '_')
                        if isDone: Energy = data[1][-1]
                        for filename in glob(os.path.join(root, '*/qchem.out')):
                            #try: 
                            qcdata = ParseOutput(filename)
                            if qcdata[-1]:
                                Energies = qcdata[1]
                                #print filename, Energies
                                if len(Energies) >= 2:
                                    Correction = Energies[-1] - Energies[-2]
                                    Energy += Correction
                                    #print Energy, Correction;exit()
                                else:
                                    isDone = False
                            else:
                                isDone = False
                            #except TypeError: Energy = None
                    """
                elif cwdfile == 'td.out' and 'td' in calctype:
                    Dipole = data[2]

                elif 'tddft' in calctype:
                    from QChemTDDFTParser import QChemTDDFTParser
                    TDAEnergies, TDADipole, Energies, Dipole = QChemTDDFTParser(filename).GetEnergiesAndDipole()
                    if Energies != None: Energies /= kcal_mol

                elif cwdfile == 'qchem.out':
                    """Assume all energy calculations go through the
                    CHARMM/Q-Chem interface.

                    This avoids double processing of the polarizable simulations"""
                    continue

                #Update only if job is done
                if not isDone: continue

                #Update a single state's energy
                if Energy != None:
                    try:
                        energyarray = h5data.createEArray(\
                         where = os.path.join('/' + GeometryName, Environment, 'Sites', Site), \
                         name = 'Energy', atom = tables.atom.FloatAtom(), shape = (0,), title = 'Energy',
                         expectedrows = 2, createparents = True)
                    except tables.exceptions.NodeError:
                        #Already in there
                        if DoOverwrite:
                            energyarray = h5data.getNode(os.path.join('/' + GeometryName, Environment, 'Sites', Site), 'Energy')
                        else: continue

                    while State > len(energyarray) - 1:
                        energyarray.append([0.0])

                    energyarray[State] = Energy * kcal_mol
                    print 'Archiving into ', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Energy')+'['+str(State)+']'
                    print 'Archived', energyarray, '[' + str(State) + ']'


                #Update all energies
                if Energies != None:
                    try:
                        energyarray = h5data.createEArray(\
                         where = os.path.join('/' + GeometryName, Environment, 'Sites', Site), \
                         name = 'Energy', atom = tables.atom.FloatAtom(), shape = (0,), title = 'Energy',
                         expectedrows = 2, createparents = True)
                    except tables.exceptions.NodeError:
                        #Already in there
                        if DoOverwrite:
                            energyarray = h5data.getNode(os.path.join('/' + GeometryName, Environment, 'Sites', Site), 'Energy')
                        else: continue

                    while len(Energies) > len(energyarray):
                        energyarray.append([0.0])

                    print 'Archiving into ', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Energy')
                    for State, Energy in enumerate(Energies):
                        energyarray[State] = Energy * kcal_mol
                        print 'Archived', Energy, '[' + str(State) + ']'



                #Update a transition dipole
                if Dipole != None:
                    try:
                        location = os.path.join(os.sep, GeometryName, Environment, 'Sites', Site)
                        dipolearray = h5data.createCArray(where = location,
                         name = 'Dipole', atom = tables.atom.FloatAtom(), shape = Dipole.shape, createparents = True)
                    except tables.exceptions.NodeError:
                        #Already in there
                        if DoOverwrite:
                            dipolearray = h5data.getNode(location, 'Dipole')
                        else:
                            print 'Skipping data that is already present', location
                            continue

                    dipolearray[:,:,:] = Dipole
                    print 'Archiving', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Dipole')
                    print 'Archived transition dipole matrix'
                    print dipolearray[:,:,:]


                if 'tddft' in calctype: #Store also TDA energies and dipoles

                    try:
                        energyarray = h5data.createEArray(\
                         where = os.path.join('/' + GeometryName, Environment, 'Sites', Site), \
                         name = 'TDAEnergy', atom = tables.atom.FloatAtom(), shape = (0,), title = 'Energy',
                         expectedrows = 2, createparents = True)
                    except tables.exceptions.NodeError:
                        #Already in there
                        if DoOverwrite:
                            energyarray = h5data.getNode(os.path.join('/' + GeometryName, Environment, 'Sites', Site), 'Energy')
                        else: continue

                    while len(TDAEnergies) > len(energyarray):
                        energyarray.append([0.0])

                    energyarray[:] = TDAEnergies
                    print 'Archiving into ', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Energy')
                    for State, Energy in enumerate(TDAEnergies):
                        energyarray[State] = Energy
                        print 'Archived', Energy, '[' + str(State) + ']'


                    try:
                        location = os.path.join(os.sep, GeometryName, Environment, 'Sites', Site)
                        dipolearray = h5data.createCArray(where = location,
                         name = 'TDADipole', atom = tables.atom.FloatAtom(), shape = TDADipole.shape, createparents = True)
                    except tables.exceptions.NodeError:
                        #Already in there
                        if DoOverwrite:
                            dipolearray = h5data.getNode(location, 'Dipole')
                        else: continue

                    dipolearray[:,:,:] = TDADipole
                    print 'Archiving', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Dipole')
                    print 'Archived transition dipole matrix'
                    print dipolearray[:,:,:]


                #If specified, clean up
                if Nuke:
                    print 'Deleting', root
                    shutil.rmtree(root)



    #Update history of the database 
    if h5data.isUndoEnabled() and DoCheckPoint:
        print 'Checkpointing database at mark', h5data.mark()
    h5data.close()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Archives data into HDF5 database')
    parser.add_argument('workingdir', nargs = '?', default = '.', help = 'Root to recurse from')
    parser.add_argument('--h5file', action = 'store', default = 'h2pc-data.h5', help = 'Name of HDF5 database file')
    parser.add_argument('--nuke', action = 'store_true', default = False, help = 'Deletes files after parsing')
    parser.add_argument('--nukeinc', action = 'store_true', default = False, help = 'Deletes incomplete output files')
    parser.add_argument('--checkpoint', action = 'store_true', default = False, help = 'Create new Pytables checkpoint in HDF5 file')
    parser.add_argument('--force', action = 'store_true', default = False, help = 'Forces reloading of existing data')
    args = parser.parse_args()


    LoadEmUp(args.workingdir, args.h5file, args.nuke, args.nukeinc, args.checkpoint, args.force)
