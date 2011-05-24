#!/usr/env/bin python
"""
Routines for parsing CHARMM and Q-Chem output files

If executed as a standalone script, will run :py:func:`ParseOutput`
on each file that is specified on the command line. (Wildcards
supported.)

.. versionadded:: 0.1
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

from Units import kcal_mol

def ParseOutput(filename):
    """
    Parses a :term:`CHARMM` or :term:`Q-Chem` output file containing
    :term:`QM/MM` information

    :param string filename: Name of output file to parse.
    :returns: A tuple containing:

         #. (integer) Number of :term:`CHARMM` iterations, 
         #. (list of floats) History of QM energies in kcal/mol
            (NOTE: *not* atomic units!),
         #. (:py:class:`numpy.ndarray`) transition dipole matrix,
         #. (integer) Q-Chem SCF convergence threshold n in 10^-n,
         #. (Boolean) Flag indicating unusual energy convergence,
         #. (Boolean) Flag indicating unusual RMS gradient convergence,
         #. (Boolean) Flag indicating if :term:`QM/MM` calculation is complete.

    Heuristics for unusual behavior
         #. Unusual energy convergence is defined as energy that is not 
            monotonically decreasing.
         #. Unusual RMS gradient convergence is defined as having at least
            a spike of at least twice the height of the average of the GRMS
            immediately before and after.

    .. versionadded:: 0.1
    """

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
                e = float(l.split()[1]) / kcal_mol
                e_history.append(e)

            if 'PURIFY final SCF energy' in l:
                e = float(l.split()[-1]) / kcal_mol
                e_history.append(e)

            if '*** MISSION COMPLETED -- STARFLEET OUT ***' in l:
                isdone = True

            #For batch inputs, lines of the form 'Job m of n' will be emitted
            if 'Job' in l and 'of' in l:
                isdone = False

            #Read transition dipole moment
            if 'dipole x component in diabatic basis' in l:
                tdip = None
                q = QChemIO.QChemOutput(filename)
                try:
                    x = q.ReadMatrix('dipole x component in diabatic basis')[0]
                    y = q.ReadMatrix('dipole y component in diabatic basis')[0]
                    z = q.ReadMatrix('dipole z component in diabatic basis')[0]
                    tdip = numpy.asarray([x, y, z])
                except (ValueError, TypeError):
                    #For some reason, failed to parse
                    print 'Could not parse transition dipole in', filename
                    tdip = None
                #break



    #Check for large fluctuations in RMS gradient
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
    from glob import glob
    for blob in glob(sys.argv[1:]):
        for fname in blob:
            print ParseOutput(fname)
