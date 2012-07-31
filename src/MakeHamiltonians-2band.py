#!/usr/bin/env python
# pylint: disable=R0903,W0212
"""
Constructs tight-binding one-exciton Hamiltonian from QM/MM data in HDF5 file.
"""

from numpy import arccos, array, average, dot, empty, NAN, pi, zeros
from numpy.linalg import norm
from Units import Angstrom, eV

try:
    import tables
except ImportError:
    import sys
    sys.stderr.write('This file requires PyTables to be installed.\n')
    exit()

class StateData:
    """
    Data container for excited state information
    """
    def __init__(self, Energy = None, Dipole = None):
        if Energy is not None:
            self.Energy = Energy
        if Dipole is not None:
            self.Dipole = Dipole

def ExtractEnergyAndTdip(h5filename = 'h2pc-data.h5', numStates = 2):
    """
    Iterates over HDF5 data file to extract excitation energies and transition
    dipoles.

    For phthalocyanines or other molecules with two NQ and two NR atoms, will
    attempt to disentangle Qx and Qy states and return them as index 0 and 1
    respectively.
    """
    data = dict()
    h5data = tables.openFile(h5filename, mode = 'r')
    for node in h5data.walkNodes():
        if node._v_name != 'Dipole':
            continue

        #Look up geometry
        path = node._v_pathname.split('/')
        resid  = path[-2]
        snapshot = path[1]
        snapshotidx = int(snapshot[5:])
        resididx = int(resid)

        geometry = h5data.root._f_getChild(snapshot).CHARMM_CARD
        
        #Extract energies
        energies = h5data.getNode('/'.join(path[:-1]+['Energy']))
        excites = (energies - energies[0])[1:(1+numStates)]
        #Extract transition dipoles
        dipoles = node[:, 0, 1:(1+numStates)].transpose()
   
        # Now, label the states as Qx (0) or Qy (1).
        StateAssignments = range(numStates)

        # Generally, Qx is lower in energy, but this is not enough.
        # Consider the transition dipoles to them: by definition,
        #      Qx has a transition dipole parallel to the protonated nitrogens
        #  and Qy has a transition dipole parallel to the unprotonated nitrogens

        #Extract basic, protonated hydrogens
        atoms1 = [ x['Coord'] for x in geometry.iterrows() \
                  if x['ResID'] == resid and 'NQ' in x['Type'] ]
        
        #Extract basic, unprotonated hydrogens
        atoms2 = [ x['Coord'] for x in geometry.iterrows() \
                  if x['ResID'] == resid and 'NR' in x['Type'] ]

        #Do analysis only for phthalocyanines
        if len(atoms1) == len(atoms2) == 2:
            atomvec1 = atoms1[1]-atoms1[0]
            atomvec1 /= norm(atomvec1)
            atomvec2 = atoms2[1]-atoms2[0]
            atomvec2 /= norm(atomvec2)
   
            angles = []
            for i in range(numStates):
                angles.append(dot(atomvec1, dipoles[i])/norm(dipoles[i]))
                angles.append(dot(atomvec2, dipoles[i])/norm(dipoles[i]))
            angles = arccos(array(angles))
            for i in range(numStates*2):
                angles[i] = min(angles[i], pi-angles[i])

            #Here, reassign lowest two excited states so that Qx = 0, Qy = 1
            if numStates >= 2:
                #Is transition dipole is more parallel to NH-NH axis than N-N?
                if angles[0] < angles[1]:
                    StateAssignments[0] = 0
                else:
                    StateAssignments[0] = 1 
                if angles[2] < angles[3]:
                    StateAssignments[1] = 0
                else:
                    StateAssignments[1] = 1
                #If need be, check which transition dipole is more parallel
                if StateAssignments[0] == StateAssignments[1] == 1:
                    if angles[1] < angles[3]:
                        StateAssignments[1] = 0
                    else:
                        StateAssignments[0] = 0
                elif StateAssignments[0] == StateAssignments[1] == 0:
                    if angles[0] < angles[2]:
                        StateAssignments[1] = 1
                    else:
                        StateAssignments[0] = 1

        print node._v_pathname, excites/eV, angles/pi, StateAssignments
        assert StateAssignments[0] != StateAssignments[1]

        data[snapshotidx, resididx] = [StateData(excites[idx], dipoles[idx]) 
                                         for idx in StateAssignments]
    
    return data



def ExtractMonomerDistances(h5file = 'h2pc-data.h5', data = None):
    """
    Iterates over HDF5 data file to extract centroids of molecules.
    """
    distances = {}
    h5data = tables.openFile(h5file, mode = 'r')

    snapshots = set([x[0] for x in data.keys()])
    for snapshot in snapshots:
        geometry = h5data.root._f_getChild('h2zn-%s' % snapshot).CHARMM_CARD
        resids = [x[1] for x in data.keys() if x[0] == snapshot]
        for resid in resids:
            #Find centroid of molecule
            atoms = [ x['Coord'] for x in geometry.iterrows() \
                if x['ResID'] == str(resid)]
            centroid = average(atoms, 0) * Angstrom
            distances[snapshot, resid] = centroid
            print snapshot, resid, centroid[0], centroid[1], \
                centroid[2]

    return distances

def ForsterDipoleCoupling(dipole1, dipole2, nR, distance):
    """
    Calculates Forster dipole coupling matrix element

    nR: unit vector pointing along R
    distance: magnitude of dipole-dipole separation
    """
    return (3 * dot(dipole1, nR) * dot(dipole2, nR) - \
               dot(dipole1, dipole2)) / distance**3

def WriteForsterDipoleCoupledHamiltonian(data, distances,
                                         MatlabFilename = 'Hamiltonians.mat'):
    """
    Writes exciton Hamiltonians to Matlab file
    """

    #How many states are we reading? Count length of data
    for states in data.values():
        num_states = len(states)
        break
    #Complication: must map the arbitrary labels for snapshots and resids
    #               onto 0...N-1 respectively.
    snapshots = sorted(set([x[0] for x in data.keys()]))
    resids    = sorted(set([x[1] for x in data.keys()]))
    system_size = len(resids)
    num_samples = len(snapshots)

    #Initialize to NAN
    hamiltonians = empty((system_size*num_states, system_size*num_states,
        num_samples))
    hamiltonians.fill(NAN)
    coordinates = empty((system_size, 3, num_samples))
    coordinates.fill(NAN)
    dipoles = empty((system_size, 3, num_samples, num_states))
    dipoles.fill(NAN)

    for (snapshot, resid), states in data.items():
        snapshot_id = snapshots.index(snapshot)
        resid_id = resids.index(resid)
        for state_id, state in enumerate(states):
            dipoles[resid_id, :, snapshot_id, state_id] = state.Dipole
        coordinates[resid_id, :, snapshot_id] = distances[snapshot, resid]

    for (snapshot1, resid1), states1 in data.items():
        resid_id1 = resids.index(resid1)
        snapshot_id = snapshots.index(snapshot1)
        for (snapshot2, resid2), states2 in data.items():
            if snapshot1 != snapshot2:
                continue

            submatrix = zeros((num_states, num_states))
            if resid1 == resid2:
                for state_id, state in enumerate(states1):
                    submatrix[state_id, state_id] = state.Energy
            else:
                R = distances[snapshot1, resid1] - distances[snapshot1, resid2]
                distance = norm(R)
                nR = R / distance
                for state_id1, state1 in enumerate(states1):
                    dipole1 = state1.Dipole
                    for state_id2, state2 in enumerate(states2):
                        #Here we compute the Forster coupling matrix element
                        #which assumes dipole---dipole coupling
                        submatrix[state_id1, state_id2] = \
                            ForsterDipoleCoupling(dipole1, state2.Dipole,
                                    nR, distance)

            resid_id2 = resids.index(resid2)
            hamiltonians[resid_id1*num_states:(resid_id1+1)*num_states,
                         resid_id2*num_states:(resid_id2+1)*num_states,
                         snapshot_id] = submatrix
            print snapshot1, resid1, resid2, submatrix.flatten()

    #Write Matlab file
    from scipy.io import savemat
    savemat(MatlabFilename, {'Hamiltonians':hamiltonians,
                'Coordinates':coordinates, 'Dipoles':dipoles,
                'Snapshots':snapshots, 'Resids':resids})


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = 'Generates Hamiltonians')
    parser.add_argument('--h5file', action = 'store', default = 'znpc-data.h5',
                        help = 'Name of HDF5 database file')
    args = parser.parse_args()

    nStates = 2
    print '\n\nEnergy and Transition dipole moment\n'
    statedata = ExtractEnergyAndTdip(args.h5file, nStates)
    print '\n\nCentroid positions\n'
    thedistances = ExtractMonomerDistances(args.h5file, statedata)
    print '\n\nHamiltonian subblocks\n'
    WriteForsterDipoleCoupledHamiltonian(statedata, thedistances,
            MatlabFilename='Hamiltonians-%dband.mat' % nStates)
    
