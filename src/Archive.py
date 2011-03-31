#!/usr/bin/env python
"""
Archiver for adding data into a HDF5 database
"""

from glob import glob
import numpy
import sys
import os.path
import shutil

Nuke = True #False #True#False#True #Delete files after successful archival
DoCheckPoint = False

try:
    import tables
except ImportError:
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()

#My custom module
try:
    import QChemIO
except ImportError:
    import sys
    sys.stderr.write("""Could not find QChemIO module.
Please check that PYTHONPATH is set correctly.
""")
    exit()

#Internal modules
from ParseOutput import ParseOutput, Hartree_to_kcal_mol



class CHARMM_CARD(tables.IsDescription):
    #Title = tables.StringAtom(1024)
    #NumAtoms = tables.UIntAtom()
    AtomNo = tables.UInt64Col()
    ResidNo = tables.UInt64Col()
    Res = tables.StringCol(4)
    Type = tables.StringCol(4)
    Coord = tables.Float32Col(shape = 3,)
    SegID = tables.StringCol(4)
    ResID = tables.StringCol(4)
    Weighting = tables.Float32Col()



def LoadCHARMM_CARD(h5table, CARDfile):
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

    assert thisnumatoms == numatoms, 'Warning, did not read in expected number of atoms'
    print 'Loaded CHARMM file', CARDfile, 'into', h5table
    h5table.flush()



def LoadEmUp(path):
    filename = 'h2pc-data.h5'
    compress = tables.Filters(complevel = 9, complib = 'zlib')
    h5data = tables.openFile(filename, mode = 'a', filters = compress,
        title = 'H2Pc')
    if not h5data.isUndoEnabled(): h5data.enableUndo(filters = compress)
    for root, _, files in os.walk(path):
        #Create recursive structure
        #h5data.createGroup(h5data.root, root)
        #If root not in datafile, create groups to reproduce tree
        for cwdfile in files:

            ##############################
            # Parse CHARMM coordinate file
            ##############################

            if False and not os.path.islink(os.path.join(root, cwdfile)) \
                and len(cwdfile) > 4 \
                and cwdfile[-4:] == '.cor' \
                and cwdfile[-6:-4] != '.d': #Is a CHARMM coordinate file
                #Add coordinates to H5file
                #TODO replace with a real way to test if node exists
                try: #TODO THIS FAILED TO DETECT SYMLINK.
                    cwdh5name = cwdfile[:-4].replace('-', '_')
                    thischarmmcard = h5data.createTable(
                        where = os.path.join('/' + root, cwdh5name), \
                        name = 'CHARMM_CARD', \
                        description = CHARMM_CARD, \
                        title = 'CHARMM CARD Coordinates', createparents = True)
                    LoadCHARMM_CARD(thischarmmcard, cwdfile)
                except tables.exceptions.NodeError:
                    pass #print cwdfile, 'already exists; skipping.'

            # Save the one with Drudes
            if False and not os.path.islink(os.path.join(root, cwdfile)) \
                and len(cwdfile) > 6 \
                and cwdfile[-6:] == '.d.cor':
                #Add coordinates to H5file
                #TODO replace with a real way to test if node exists
                try: #TODO THIS FAILED TO DETECT SYMLINK.
                    cwdh5name = cwdfile[:-6].replace('-', '_')
                    thischarmmcard = h5data.createTable(
                        where = os.path.join('/' + root, cwdh5name), \
                        name = 'CHARMM_CARD_WITHDRUDE', \
                        description = CHARMM_CARD, \
                        title = 'CHARMM CARD Coordinates with Drude particles', createparents = True)
                    LoadCHARMM_CARD(thischarmmcard, cwdfile)
                except tables.exceptions.NodeError:
                    pass #print cwdfile, 'already exists; skipping.'



            ####################################
            # Parse CHARMM or Q-Chem output file
            ####################################
            if not os.path.islink(os.path.join(root, cwdfile)) \
                and len(cwdfile) > 4 \
                and cwdfile[-4:] == '.out':
                filename = os.path.join(root, cwdfile)
                calctype = os.path.split(root)[1]
                GeometryName = Site = State = Energy = Dipole = Environment = None
                data = ParseOutput(filename)
                isDone = data[-1]
                if not isDone:
                    print 'Skipping incomplete calculation in', filename
                    #if Nuke:
                    #    print 'Deleting', root
                    #    shutil.rmtree(root)
                elif calctype == 'gs':
                    Environment = 'WithDrude'
                    State = 0
                    Site = root.split(os.sep)[-2]
                    GeometryName = root.split(os.sep)[-3]
                    if isDone: Energy = data[1][-1]
                    #TODO Record Drude particle positions
                elif calctype == 'dscf':
                    Environment = 'WithDrude'
                    State = 1
                    Site = root.split(os.sep)[-2]
                    GeometryName = root.split(os.sep)[-3]#.replace('-', '_')
                    if isDone: Energy = data[1][-1]
                    #TODO Record Drude particle positions
                elif calctype == 'gs-nonpol':
                    Environment = 'Fixed'
                    State = 0
                    Site = root.split(os.sep)[-2]
                    GeometryName = root.split(os.sep)[-3]#.replace('-', '_')
                    if isDone: Energy = data[1][-1]
                elif calctype == 'dscf-nonpol':
                    Environment = 'Fixed'
                    State = 1
                    Site = root.split(os.sep)[-2]
                    GeometryName = root.split(os.sep)[-3]#.replace('-', '_')
                    if isDone: Energy = data[1][-1]
                    #This fixes a bug in Q-Chem wrapper script that forgot to parse nonpolarizable
                    #delta-SCF correctly
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
                elif cwdfile == 'qchem.out':
                    pass
                elif cwdfile == 'td.out':
                    calctype = os.path.split(root)[1]
                    Site = root.split(os.sep)[-2]
                    GeometryName = root.split(os.sep)[-3]#.replace('-', '_')
                    #calctype = os.path.split(os.path.normpath(os.path.join(root, '..')))[1]
                    if calctype == 'td-nonpol':
                        Environment = 'Fixed'
                        q = QChemIO.QChemOutput(filename)
                        #Dipole is now a 3-tensor with coordinates (direction, state, state)
                        try:
                            Dipole = numpy.array([q.ReadMatrix('dipole x component in diabatic basis')[0], \
                            q.ReadMatrix('dipole y component in diabatic basis')[0], \
                            q.ReadMatrix('dipole z component in diabatic basis')[0]])
                        except TypeError: pass
                    elif calctype == 'td':
                        Environment = 'WithDrude'
                        q = QChemIO.QChemOutput(filename)
                        #Dipole is now a 3-tensor with coordinates (direction, state, state)
                        try:
                            Dipole = numpy.array([q.ReadMatrix('dipole x component in diabatic basis')[0], \
                            q.ReadMatrix('dipole y component in diabatic basis')[0], \
                            q.ReadMatrix('dipole z component in diabatic basis')[0]])
                        except TypeError: pass
                    else:
                        print 'Unknown job type', calctype, 'and file', cwdfile
                else:
                    print 'Unknown job type', calctype, 'and file', cwdfile

                #print 'Data', filename, GeometryName, Site, Environment, State, Energy, Dipole != None
                if isDone and Energy != None:
                    try:
                        energyarray = h5data.createEArray(\
                         where = os.path.join('/' + GeometryName, Environment, 'Sites', Site), \
                         name = 'Energy', atom = tables.atom.FloatAtom(), shape = (0,), title = 'Energy',
                         expectedrows = 2, createparents = True)
                    except tables.exceptions.NaturalNameWarning:
                        pass
                    except tables.exceptions.NodeError:
                        #Already in there
                        energyarray = h5data.getNode(os.path.join('/' + GeometryName, Environment, 'Sites', Site), 'Energy')
                    while State > len(energyarray) - 1:
                        energyarray.append([0.0])
                    energyarray[State] = Energy / Hartree_to_kcal_mol
                    #print 'Archived', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Energy')+'['+str(State)+']'
                    print 'Archived', energyarray, '[' + str(State) + ']'
                    if Nuke:
                        print 'Deleting', root
                        shutil.rmtree(root)
                if isDone and Dipole != None:
                    try:
                        dipolearray = h5data.createCArray(\
                         where = os.path.join(os.sep, GeometryName, Environment, 'Sites', Site), \
                         name = 'Dipole', atom = tables.atom.FloatAtom(), shape = Dipole.shape)
                        dipolearray = Dipole
                        #print 'Archived', os.path.join(os.sep, GeometryName, Environment, 'Sites', Site, 'Dipole')
                        print 'Archived', dipolearray
                    except tables.exceptions.NaturalNameWarning:
                        pass
                    except tables.exceptions.NodeError:
                        pass #Already in there
                    if Nuke:
                        print 'Deleting', root
                        shutil.rmtree(root)

    #Update history of the database 
    if h5data.isUndoEnabled() and DoCheckPoint: print 'Checkpointing database at mark', h5data.mark()
    h5data.close()

if __name__ == '__main__':
    LoadEmUp('.')