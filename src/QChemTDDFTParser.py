#!/usr/bin/env python

#TODO Integrate into QChemIO!
class DataRecorder:
    def __init__(self):
        self.ground = []
        self.tddft_tda = []
        self.tddft = []

    def Record(self, mode, StateIdx, Energy, Multiplicity, Dipole = None, OscillatorStrength = None):
        """
        Record parsed TDDFT data:
        :param string moder: Type of data ('dft', 'tddft-tda', or 'tddft')
        :param integer StateIdx: Index of state (0 = ground, 1 = first excited...)
        :param float Energy: Energy in atomic units
        :param integer Multiplicity: Spin multiplicity (1 = singlet)
        :param float(3) Dipole: Transition dipole from ground state (Optional)
        :param float OscillatorStrength: Oscillator strength of transition (Optional)
        """

        assert mode in ['dft', 'tddft-tda', 'tddft']
        if mode == 'dft':
            self.ground.append([Energy, Multiplicity])
        elif mode == 'tddft-tda':
            self.tddft_tda.append([StateIdx, Energy, Multiplicity, Dipole, OscillatorStrength])
        elif mode == 'tddft':
            self.tddft.append([StateIdx, Energy, Multiplicity, Dipole, OscillatorStrength])


    def GetEnergiesAndDipole(self):
        "If TDDFT available, return it else return TDA energies"
        try:
            maxstates = max([datum[0] for datum in self.tddft_tda]) + 1
        except ValueError:
            #Not enough data
            return (None, )*4

        import numpy
        Energies = numpy.ndarray((maxstates,))
        Dipole = numpy.ndarray((3, maxstates, maxstates))
        TDAEnergies = numpy.ndarray((maxstates,))
        TDADipole = numpy.ndarray((3, maxstates, maxstates))

        #Initialize to NaN
        #TODO Update all code to use NaN convention
        TDAEnergies[:] = Energies[:] = numpy.NaN
        TDADipole[:,:,:] = Dipole[:,:,:] = numpy.NaN

        try:
            Energies[0] = self.ground[-1][0]
            TDAEnergies[0] = self.ground[-1][0]
        except IndexError:
            print 'Could not find ground state energy'
            return (None, )*4

        for stateid, energy, _, dipole, _ in self.tddft_tda:
            TDAEnergies[stateid] = energy
            TDADipole[:, 0, stateid] = TDADipole[:, stateid, 0] = dipole

        for stateid, energy, _, dipole, _ in self.tddft:
            Energies[stateid] = energy
            Dipole[:, 0, stateid] = Dipole[:, stateid, 0] = dipole

        return TDAEnergies, TDADipole, Energies, Dipole



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
                    Multiplicity = 1
                elif t[1] == 'Triplet':
                    Multiplicity = 3
                else:
                    raise ValueError, 'Unknown multiplicity:\n'+line
                
            elif 'Trans. Mom.' in line:
                Dipole = tuple([float(t[n]) for n in (2, 4, 6)])
            elif 'Strength' in line:
                OscillatorStrength = float(t[2])
                #Skip parsing of transition amplitudes
            elif len(t) == 0: #Blank line
                Data.Record(mode, StateIdx, Energy, Multiplicity, Dipole, OscillatorStrength)
                del StateIdx, Energy, Multiplicity, Dipole, OscillatorStrength
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
