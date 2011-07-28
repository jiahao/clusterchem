#!/usr/bin/env python

#TODO Integrate into QChemIO!
import numpy
class DataRecorder:
    def __init__(self):
        self.ground = []
        self.tddft_tda = []
        self.tddft = []

    def Record(self, mode, StateIdx, Energy, Multiplicity, Dipole = None,
               OscillatorStrength = None, Amplitudes = None):
        """
        Record parsed TDDFT data:
        :param string moder: Type of data ('dft', 'tddft-tda', or 'tddft')
        :param integer StateIdx: Index of state (0 = ground, 1 = first excited
           ...)
        :param float Energy: Energy in atomic units
        :param integer Multiplicity: Spin multiplicity (1 = singlet)
        :param float(3) Dipole: Transition dipole from ground state (Optional)
        :param float OscillatorStrength: Oscillator strength of transition
        (Optional)
        :param dict Amplitudes: Transition amplitudes (Optional)
        :type Amplitudes: the key of the dictionary is a integer 2-tuple
        representing the originating and destination orbitals respectively.
        The value is the amplitude.
        """

        assert mode in ['dft', 'tddft-tda', 'tddft']
        if mode == 'dft':
            self.ground.append([Energy, Multiplicity])
        elif mode == 'tddft-tda':
            self.tddft_tda.append([StateIdx, Energy, Multiplicity, Dipole,
                                   OscillatorStrength, Amplitudes])
        elif mode == 'tddft':
            self.tddft.append([StateIdx, Energy, Multiplicity, Dipole,
                               OscillatorStrength, Amplitudes])


    def GetEnergiesAndDipole(self):
        try:
            maxstates = max([datum[0] for datum in self.tddft_tda]) + 1
        except ValueError:
            #Not enough data
            return (None, )*4

        Energies = numpy.ndarray((maxstates,))
        Dipole = numpy.ndarray((3, maxstates, maxstates))
        TDAEnergies = numpy.ndarray((maxstates,))
        TDADipole = numpy.ndarray((3, maxstates, maxstates))

        #Initialize to NaN
        #TODO Update all code to use NaN convention
        TDAEnergies[:] = Energies[:] = numpy.NaN
        TDADipole[:, :, :] = Dipole[:, :, :] = numpy.NaN

        try:
            Energies[0] = self.ground[-1][0]
            TDAEnergies[0] = self.ground[-1][0]
        except IndexError:
            print 'Could not find ground state energy'
            return (None, )*4

        for stateid, energy, _, dipole, _, _ in self.tddft_tda:
            TDAEnergies[stateid] = energy
            TDADipole[:, 0, stateid] = TDADipole[:, stateid, 0] = dipole

        for stateid, energy, _, dipole, _, _ in self.tddft:
            Energies[stateid] = energy
            Dipole[:, 0, stateid] = Dipole[:, stateid, 0] = dipole

        return TDAEnergies, TDADipole, Energies, Dipole


    def GetOtherData(self):
        "If TDDFT available, return it else return TDA energies"
        try:
            maxstates = max([datum[0] for datum in self.tddft_tda]) + 1
        except ValueError:
            #Not enough data
            return (None, )*4

        TDAMultiplicities = numpy.ndarray((maxstates,), numpy.int32)
        TDAOscillatorStrengths = numpy.ndarray((maxstates,))
        Multiplicities = numpy.ndarray((maxstates,), numpy.int32)
        OscillatorStrengths = numpy.ndarray((maxstates,))

        #XXX Initialize ground state to be singlet
        Multiplicities[0] = 0
        TDAMultiplicities[0] = 0
        #TDATransitionAmplitudes = numpy.ndarray((maxstates, maxamplitudestda))
        #TransitionAmplitudes = numpy.ndarray((maxstates, maxamplitudestddft))

        for stateid, _, multiplicity, _, strength, amplitudes in \
                self.tddft_tda:
            TDAMultiplicities[stateid] = multiplicity
            TDAOscillatorStrengths[stateid] = strength
        
        for stateid, _, multiplicity, _, strength, amplitudes in self.tddft:
            Multiplicities[stateid] = multiplicity
            OscillatorStrengths[stateid] = strength

        return TDAMultiplicities, TDAOscillatorStrengths, Multiplicities, \
            OscillatorStrengths



def QChemTDDFTParser(filename):
    mode = 'scan'
    skiplines = 0

    Data = DataRecorder()

    for line in open(filename):
        #print mode, line,
        if skiplines > 0:
            skiplines -= 1
            continue

        t = line.split()

        if mode == 'scan':
            if 'Convergence criterion met' in line:
                Data.Record('dft', 0, Energy = float(t[1]), Multiplicity = 1)
                #XXX Hard-coded singlet ground state!!

            elif 'TDDFT/TDA Excitation Energies' in line:
                skiplines = 2
                mode = 'tddft-tda'
            elif 'TDDFT Excitation Energies' in line:
                skiplines = 2
                mode = 'tddft'

        elif 'tddft' in mode:
            if 'Total energy for state' in line:
                Energy = float(t[-1])
                StateIdx = int(t[4][:1])

            elif 'Multiplicity' in line:
                if t[1] == 'Singlet':
                    Multiplicity = 0
                elif t[1] == 'Triplet':
                    Multiplicity = 2
                else:
                    raise ValueError, 'Unknown multiplicity:\n'+line
                
            elif 'Trans. Mom.' in line:
                Dipole = tuple([float(t[n]) for n in (2, 4, 6)])
            elif 'Strength' in line:
                OscillatorStrength = float(t[2])
                Amplitudes = {}
            #elif '-->' in line: #transition amplitudes
            #    #    D(130) --> V(  1) amplitude =  0.2599
            #    #012345678901234567890
            #    print line,
            #    orbital_from = int(line[6:9])
            #    orbital_to   = int(line[17:20])
            #    amplitude = float(t[-1])
            #    Amplitudes[(orbital_from, orbital_to)] = amplitude
            elif len(t) == 0: #Blank line
                Data.Record(mode, StateIdx, Energy, Multiplicity, Dipole, 
                            OscillatorStrength, Amplitudes)
                del StateIdx, Energy, Multiplicity, Dipole, \
                    OscillatorStrength, Amplitudes
            elif '---------------------------------------------------' in line:
                mode = 'scan'

    return Data



if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = 'qchem.out'

    QChemTDDFTParser(filename)
