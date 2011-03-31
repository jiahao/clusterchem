#!/usr/env/bin python
"""
Routines for parsing CHARMM and Q-Chem output files
"""

import numpy, sys
#My custom module
try:
    import QChemIO
except ImportError:
    sys.stderr.write("""Could not find QChemIO module.
Please check that PYTHONPATH is set correctly.
""")
    exit()

Hartree_to_kcal_mol = 627.5095

def ParseOutput(filename):
    """Parses a CHARMM or Q-Chem output file"""

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
                q = QChemIO.QChemOutput(filename)
                x = q.ReadMatrix('dipole x component in diabatic basis')[0]
                y = q.ReadMatrix('dipole y component in diabatic basis')[0]
                z = q.ReadMatrix('dipole z component in diabatic basis')[0]
                tdip = numpy.asarray([x, y, z])
                break

    #Check for large fluctuations in RMS gradient
    #The heuristic here is to detect spikes of at least twice the
    #height of the average of the GRMS immediately before and after
    isgrmsweird = False
    for i in range(1, charmmiter - 2):
        try:
            if grms_history[i] > grms_history[i + 1] + grms_history[i - 1]:
                isgrmsweird = True
                break
        except IndexError:
            pass

    #Check for non-monotonic convergence in energy
    isenergyweird = False
    tol = 1e-2
    for i in range(len(lastmm_e_history) - 2):
        if e_history[i] + tol < e_history[i + 1]:
            isenergyweird = True
            break

    return charmmiter, e_history, tdip, qm_iterthresh, isenergyweird, isgrmsweird, isdone


if __name__ == '__main__':
    for fname in sys.argv[1:]:
        print ParseOutput(fname)
