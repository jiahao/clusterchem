#!/usr/bin/env python
"""
:term:`Q-Chem` interface module
Creates :term:`Q-Chem` input decks for term:`QM/MM` calculations.

If run as a standalone script, will execute QChemInputForElectrostaticEmbedding
passing command-line arguments to it.

Requires the homebrew :py:mod:`QChemIO` module.
"""

from copy import deepcopy
import QChemIO

def QChemInputForElectrostaticEmbedding(ResID, CHARMM_CARD_file,
    CHARMM_RTFile, InputDeck = None, QChemInputFileName = None):
    '''
    Create a Q-Chem input deck to calculate the electrostatic field at each MM
    atom due to the quantum-mechanical region.

    :parameter ResID: CHARMM ResID which defines the embedded QM region 
    :type ResID: string

    :parameter CHARMM_CARD_file: Name of CHARMM CARD file containing molecular
    coordinates
    :type CHARMM_CARD_file: string

    :parameter CHARMM_RTFile: Name of CHARMM Residue Topology File containing
    force field parameters (notably, atomic charges)
    :type CHARMM_RTFile: string or iterable

    :parameter InputDeck: Q-Chem input deck. (Optional.) \
    If value is a string, is treated as a filename for QChemInput object. \
    If it is a QChemInput object, then it will be modified. \
    If not specified, a new QChemInput object will be created.

    :type InputDeck: string or QChemInputObject

    :parameter NewFileName: A new filename for InputDeck. (Optional.)
    :type NewFileName: string

    :returns: A Q-Chem input deck that is ready to be written and executed.
    :rtype: :py:class:`QChemInput`

    NOTE: The file will NOT be written to disk. issue a .write() manually as
    necessary!               
    '''

    #Load fixed charge specification in RTF file
    ChargeParameter = {}
    try:
        for l in open(CHARMM_RTFile):
            if l[:4].upper() == 'ATOM':
                t = l.split()
                AtomType, Charge = t[1].upper(), float(t[3])
                ChargeParameter[AtomType] = Charge

    except TypeError:
        #Not a string, assume it's an iterable of strings
        for filename in CHARMM_RTFile:
            for l in open(filename):
                if l[:4].upper() == 'ATOM':
                    t = l.split()
                    AtomType, Charge = t[1].upper(), float(t[3])
                    ChargeParameter[AtomType] = Charge

    #Generate $external_charges and $molecule blocks for Q-Chem input
    QMbuf = ['0 1'] #XXX Hard-coded charge and spin!
    MMbuf = []
    for l in open(CHARMM_CARD_file):
        t = l.split()
        if len(t) > 8: #Need at least 9 columns
            AtomType, x, y, z, thisResID = t[3].upper(), float(t[4]), \
                float(t[5]), float(t[6]), t[8]
            if thisResID == ResID: #QM region
                QMbuf.append(AtomType[0]+('%15.8f'*3) % (x, y, z))
            else:
                assert AtomType in ChargeParameter, 'ERROR: Could not find \
charge parameter for '+AtomType+' in '+CHARMM_RTFile
                charge = ChargeParameter[AtomType]
                MMbuf.append(('%15.8f\t'*4) % (x, y, z, charge)) 

    #Make Q-Chem input deck
    if InputDeck == None:
        InputDeck = QChemIO.QChemInput(QChemInputFileName)
    elif type(InputDeck) == type('1'): #A string, treat as filename
        InputDeck = QChemIO.QChemInput(InputDeck)
    elif isinstance(InputDeck, QChemIO.QChemInput):
        pass
    else:
        assert False, "I don't know what to do with InputDeck ="+str(InputDeck)
    InputDeck.GetCurrentJob().append_input_section('rem', 'igdesp\t%d\n' %
                                                   len(MMbuf))
    InputDeck.GetCurrentJob().write_input_section('molecule', '\n'.join(QMbuf),
                                                  Overwrite = True)
    InputDeck.GetCurrentJob().write_input_section('external_charges',
                                    '\n'.join(MMbuf), Overwrite = True)
    if QChemInputFileName != None: InputDeck.filename = QChemInputFileName
    return InputDeck #Warning, not written to disk yet!!



def QChemInputForTransitionDipole(InputDeck1, InputDeck2,
                                  NewFileName = 'tmp.in'):
    '''
    Create a :term:`Q-Chem` input deck to calculate the transition dipole
    moment between two electronic states.

    This makes use of CDFT-CI's ability to calculate transition dipole moments.

    :param InputDeck1: :term:`Q-Chem` input deck containing the first
    electronic state
     
    If value is a string, is treated as a filename for the QChemInput() output.
    If it is a :py:class:`QChemInput` object, then it will be modified.
    If not specified, a new :py:class:`QChemInput` object will be created.
    :type InputDeck1: string or :py:class:`QChemInput`

    :param InputDeck2: :term:`Q-Chem` input deck containing the second
    electronic state
    :type InputDeck2: string or :py:class:`QChemInput`

    :param string NewFileName: A new filename for the resultant QChemInput().
    (Default: tmp.in)

    :returns: A Q-Chem input deck that is ready to be written and executed.
    :rtype: :py:class:`QChemInput`
    
    NOTE: The file will NOT be written to disk. issue a .write() manually as
    necessary!   
    '''


    Q = QChemIO.QChemInput(NewFileName)
    del Q.jobs[0] #Delete empty job that is initialized

    if type(InputDeck1) == type('1'): #A string, treat as filename
        Q.ImportQChemOutput(InputDeck1)
    elif isinstance(InputDeck1, QChemIO.QChemInput):
        Q.jobs.append(deepcopy(InputDeck1.GetCurrentJob()))
    else:
        assert False, "I don't know what to do with InputDeck1 ="+\
            str(InputDeck1)

    Q.GetCurrentJob().append_input_section('rem', """\
symmetry                  off
sym_ignore                true
cdftci                    true
cdftci_stop               1
cdftci_skip_promolecules  true
print_input               false
qm_mm                     true
qmmm_print                true
skip_charge_self_interact true
print_orbitals            true
""")

    Q.GetCurrentJob().write_input_section('cdft', """\
0.
1. 0 0
0.
1. 0 0 s
---
0.
1. 0 0
0.
1. 0 0 s
""")

    if type(InputDeck2) == type('1'): #A string, treat as filename
        Q.ImportQChemOutput(InputDeck2)
    elif isinstance(InputDeck2, QChemIO.QChemInput):
        Q.jobs.append(deepcopy(InputDeck2.GetCurrentJob()))
    else:
        assert False, "I don't know what to do with InputDeck2 ="+\
str(InputDeck2)


    Q.GetCurrentJob().append_input_section('rem', """\
symmetry                  off
sym_ignore                true
cdftci                    true
cdftci_restart            1
cdftci_print              2
cdftci_skip_promolecules  true
print_input               false
qm_mm                     true
qmmm_print                true
skip_charge_self_interact true
""")

    Q.GetCurrentJob().write_input_section('cdft', """\
0.
1. 0 0
0.
1. 0 0 s
---
0.
1. 0 0
0.
1. 0 0 s
""")


    return Q



if __name__ == '__main__':
    import os, sys
    QCInput = QChemInputForElectrostaticEmbedding(*sys.argv[1:])
    if os.path.exists(QCInput.filename):
        print 'Output file', QCInput.filename, 'exists! Not overwriting.'
        print 'Writing to console instead.'
        print QCInput
    else:
        QCInput.write()
