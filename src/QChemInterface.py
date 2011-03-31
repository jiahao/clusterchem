#!/usr/bin/env python
"""Creates Q-Chem input scripts for computation
"""

from copy import deepcopy
import QChemIO

def QChemInputForElectrostaticEmbedding(ResID, CHARMM_CARD_file, CHARMM_RTFile, InputDeck = None, QChemInputFileName = None):
    '''
    Create a Q-Chem input deck to calculate the electrostatic field at each MM atom
    due to the quantum-mechanical region.

    Inputs
    ------
    ResID            - CHARMM ResID which defines the embedded QM region 
    CHARMM_CARD_file - CHARMM CARD file containing molecular coordinates
    CHARMM_RTFile    - CHARMM Residue Topology File
    InputDeck        - Q-Chem input deck. (Optional.)
                       If value is a string, is treated as a filename for QChemInput object.
                       If it is a QChemInput object, then it will be modified.
                       If not specified, a new QChemInput object will be created.
    NewFileName      - A new filename for InputDeck. (Optional.)

    Returns
    -------
    QChemInput       - A QChemInput object.
                       NOTE: The file will NOT be written to disk. issue a .write() manually as necessary!               
    '''
    #Load fixed charge specification in RTF file
    ChargeParameter = {}
    for l in open(CHARMM_RTFile):
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
            AtomType, x, y, z, thisResID = t[3].upper(), float(t[4]), float(t[5]), float(t[6]), t[8]
            if thisResID == ResID: #QM region
                QMbuf.append(AtomType[0]+('%15.8f'*3) % (x, y, z))
            else:
                assert AtomType in ChargeParameter, 'ERROR: Could not find charge parameter for '+AtomType+' in '+CHARMM_RTFile
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
        assert False, "I don't know what to do with InputDeck ="+repr(InputDeck)
    InputDeck.GetCurrentJob().append_input_section('rem', 'igdesp\t%d\n' % len(MMbuf))
    InputDeck.GetCurrentJob().write_input_section('molecule', '\n'.join(QMbuf), Overwrite = True)
    InputDeck.GetCurrentJob().write_input_section('external_charges', '\n'.join(MMbuf), Overwrite = True)
    if QChemInputFileName != None: InputDeck.filename = QChemInputFileName
    return InputDeck #Warning, not written to disk yet!!



def QChemInputForTransitionDipole(InputDeck1, InputDeck2, NewFileName = 'tmp.in'):
    '''
    Create a Q-Chem input deck to calculate the transition dipole moment between
    two electronic states.

    This makes use of CDFT-CI's ability to calculate transition dipole moments.

    Inputs
    ------
    InputDeck1        - Q-Chem input deck containing the first electronic state
                       If value is a string, is treated as a filename for QChemInput object.
                       If it is a QChemInput object, then it will be modified.
                       If not specified, a new QChemInput object will be created.
    InputDeck2        - Q-Chem input deck containing the second electronic state
    NewFileName       - A new filename for the resultant QChemInput. (Optional. Default: tmp.in)

    Returns
    -------
    QChemInput       - A QChemInput object.
                       NOTE: The file will NOT be written to disk. issue a .write() manually as necessary!   
            
    '''


    Q = QChemIO.QChemInput(NewFileName)

    if type(InputDeck1) == type('1'): #A string, treat as filename
        Q.ImportQChemOutput(InputDeck1)
    elif isinstance(InputDeck1, QChemIO.QChemInput):
        Q.jobs.append(deepcopy(InputDeck1.GetCurrentJob()))
    else:
        assert False, "I don't know what to do with InputDeck1 ="+repr(InputDeck1)

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
        assert False, "I don't know what to do with InputDeck2 ="+repr(InputDeck2)


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
