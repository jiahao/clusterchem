#!/usr/bin/env python

import numpy
from numpy import arccos, dot, pi
from numpy.linalg import norm
from Units import Angstrom

try:
    import tables
except ImportError:
    import sys
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()

def ExtractEnergyAndTdip(h5filename = 'h2pc-data.h5'):
    
    data = dict()
    h5data = tables.openFile(h5filename, mode = 'r')
    for node in h5data.walkNodes():
        if node._v_name == 'Dipole':
            dipole = node[:, 0, 1]
            
            #Look up geometry
            path = node._v_pathname.split('/')
            resid  = path[-2]
            snapshot = path[1]
            snapshotidx = int(snapshot[5:])
            resididx = int(resid)

            #Extract transition dipole from matrix
            try:
                geometry = h5data.getNode('/'+snapshot+'/CHARMM_CARD')
            except tables.exceptions.NoSuchNodeError: continue
            
            #Extract basic hydrogens
            atoms = [ x['Coord'] for x in geometry.iterrows() \
                      if x['ResID'] == resid and 'NQ' in x['Type'] ]

            if len(atoms) == 2: #XXX I have arbitrarily assumed that this is H2PC
                continue #silently fail
                assert len(atoms) == 2, """\
Wrong number of free base nitrogens in %r:%r
Expected 2 but found %d
Coordinates:
""" % (snapshot, resid, len(atoms)) + '\n'.join([str(x) for x in atoms])
                
                atomvec = atoms[1] - atoms[0]
    
                #Extract energies
                try:
                    energies = h5data.getNode('/'.join(path[:-1]+['Energy']))
                    excite = (energies[1] - energies[0])
                except (IndexError, tables.exceptions.NoSuchNodeError):
                    excite = 0.0
    
                if norm(dipole) != 0.0:
                    angle = arccos(dot(atomvec, dipole)/(norm(atomvec)*\
                                                             norm(dipole)))
                    angle = min(angle, pi - angle)
    
                    if len(energies) > 2: #Try next higher state
                        newexcite = (energies[2] - energies[0])
                        newdipole = node[:, 0, 2]
                        newangle = arccos(dot(atomvec, newdipole)/\
                                              (norm(atomvec)*norm(newdipole)))
                        newangle = min(newangle, pi - newangle)
                        if newangle < angle:
                            #print 'Root flip at', snapshot, resid
                            angle = newangle
                            excite = newexcite
                            dipole = newdipole
            else: #XXX No special check for the transition dipole
                try:
                    energies = h5data.getNode('/'.join(path[:-1]+['Energy']))
                    excite = (energies[1] - energies[0])
                except (IndexError, tables.exceptions.NoSuchNodeError):
                    excite = 0.0
    

            if snapshotidx not in data:
                data[snapshotidx] = dict()
            print snapshotidx, resididx, excite, dipole[0], dipole[1], \
                dipole[2]
            data[snapshotidx][resididx] = (excite, dipole)

    h5data.close()

    return data



def ExtractMonomerDistances(h5file = 'h2pc-data.h5', data = None):
    distances = {}
    h5data = tables.openFile(h5file, mode = 'r')
    for node in h5data.walkNodes():
        if node._v_name == 'CHARMM_CARD':
            #Look up geometry
            geometry = node
            path = node._v_pathname.split('/')
            snapshot = path[1]
            snapshot_idx = int(snapshot[5:])
            if snapshot_idx not in distances:
                distances[snapshot_idx] = {}

            if snapshot_idx not in data:
                print 'No data found for ', snapshot_idx
                continue

            for resid in data[snapshot_idx]:
                #Find centroid of molecule
                atoms = [ x['Coord'] for x in geometry.iterrows() \
                    if x['ResID'] == str(resid)]
                centroid = numpy.average(atoms, 0) * Angstrom
                distances[snapshot_idx][resid] = centroid
                print snapshot_idx, resid, centroid[0], centroid[1], \
                    centroid[2]
    return distances


def PrintForsterDipoleCoupledHamiltonian(data, distances,
                                         MatlabFilename = 'Hamiltonians.mat'):

    #Create sorted list to map indices
    snapshot_map = [snapshot for snapshot in data]
    snapshot_map.sort()
    snapshot_mapper = {}
    for idx, entry in enumerate(snapshot_map):
        snapshot_mapper[entry] = idx

    for snapshot, snapshot_data in data.items():
        resid_map = [resid for resid in snapshot_data]
        break

    resid_map.sort()
    resid_mapper = {}

    for idx, entry in enumerate(resid_map):
        resid_mapper[entry] = idx

    system_size = len(resid_mapper)
    num_samples = len(snapshot_mapper)
    
    hamiltonians = numpy.empty((system_size, system_size, num_samples))
    coordinates = numpy.empty((system_size, 3, num_samples))
    dipoles = numpy.empty((system_size, 3, num_samples))
    hamiltonians[:] = numpy.NAN
    coordinates[:] = numpy.NAN
    dipoles[:] = numpy.NAN

    for snapshot, snapshot_data in data.items():
        snapshot_idx = snapshot_mapper[snapshot]

        for residA, (energyA, dipoleA) in snapshot_data.items():
            if residA not in resid_mapper:
                continue
            resA_idx = resid_mapper[residA]
            for residB, (energyB, dipoleB) in snapshot_data.items():
                if residB not in resid_mapper:
                    continue
                resB_idx = resid_mapper[residB]
                if residA == residB:
                    print snapshot, residA, residB, energyA
                    hamiltonians[resA_idx, resB_idx, snapshot_idx] = energyA
                    dipoles[resA_idx, :, snapshot_idx]= dipoleA
                    coordinates[resA_idx, :, snapshot_idx]=distances[snapshot][residA]
                else:
                    R = (distances[snapshot][residA] - \
                         distances[snapshot][residB])
                    distance = norm(R)
                    nR = R / distance
                    #Here we compute the Forster coupling matrix element
                    #which assumes dipole---dipole coupling
                    V = (3 * dot(dipoleA, nR) * dot(dipoleB, nR) - \
                            dot(dipoleA, dipoleB)) / distance**3
                    hamiltonians[resA_idx, resB_idx, snapshot_idx] = V
                    print snapshot, residA, residB, V

    #Write Matlab file
    from scipy.io import savemat
    savemat(MatlabFilename, {'Hamiltonians':hamiltonians, 'Coordinates':coordinates, 'Dipoles':dipoles})

    return hamiltonians


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Generates Hamiltonians')
    parser.add_argument('--h5file', action = 'store', default = 'h2pc-data.h5',
                        help = 'Name of HDF5 database file')
    args = parser.parse_args()

    print '\n\nEnergy and Transition dipole moment\n'
    data = ExtractEnergyAndTdip(args.h5file)
    print '\n\nCentroid positions\n'
    distances = ExtractMonomerDistances(args.h5file, data)
    print '\n\nHamiltonians\n'
    hamiltonians = PrintForsterDipoleCoupledHamiltonian(data, distances)
    
