#!/usr/bin/env python

"""
Parses Q-Chem output files with FED calculations

For FED calculations, we grep TDDFT calculations
for the energies and transition dipoles of bright
excited states
"""
from numpy import array, asscalar, average, dot, ones, zeros, arccos, pi, log, where, vstack, lexsort
from numpy.linalg import norm

try:
    from openqube import *
except ImportError:
    print """Unable to load openqube Q-Chem parser. 

If you are running this script on caspian.mit.edu, run

export PYTHONPATH=$PYTHONPATH:/home/cjh/local/lib/python2.6/site-packages

and try again. Alternatively, install the python portion of openqube from

https://github.com/jiahao/openqube
"""
    exit()

#
# Parameters
#

BrightStateStrengthThreshold = 0.00000002 
ChargeChangeThreshold = 0.03 #At least this much charge must change

#
# Helper functions
#

def MomentOfInertiaTensor(R, Weights = None, Center = True):
    "Moment of inertia tensor for a bunch of points, optionally weighted"
    from numpy import zeros
    I = zeros((3,3))
    if Center == True:
        C = Centroid(R, Weights)
    elif Center is None:
        C = zeros(3)
    else:
        C = Center

    for idx, p in enumerate(R):
        weight = 1 if Weights is None else Weights[idx]
        for i in range(3):
            for j in range(3):
                I[i,j] += weight * (p[i]-C[i]) * (p[j]-C[j])
    return I

def Centroid(R, Weights = None):
    "Centroid for a bunch of points, optionally weighted"
    from numpy import zeros
    C = zeros((3,))
    
    for idx, p in enumerate(R):
        weight = 1 if Weights is None else Weights[idx]
        C += weight * p
 
    normalization = len(R) if Weights is None else sum(Weights)
    if abs(normalization) < 0.001:
        normalization = 1.0

    return C/normalization

#
# The main function
#
def ParseFED(filename, doPlot = False):
    """
    See Q-Chem 3.2 manual,
    10.17 Electronic Couplings for Electron Transfer and Energy Transfer
    """

    #Pickle
    try:
        import pickle
        f = open(filename+'.pkl', 'rb')
        Geometry = pickle.load(f)
        States = pickle.load(f)
        Charges = pickle.load(f)
        f.close()
        print 'Loaded from pickle', filename
    except (EOFError, IOError, AttributeError):
        mode = 'seek'
        skiplines = 0
        fed_trigger = False

        X = QChemOutput(filename)
        
        #TODO pickle X
        #Pickle
        #f = open(filename+'.pkl', 'wb')
        #pickle.dump(Geometry, f)
        #pickle.dump(States, f)
        #pickle.dump(QCModel, f)
        #f.close()
        #print 'Cached parsed output in', filename+'.pkl'

        #Electronic excited states
        States = [data for label, data in X.Data if 'TDDFT' in label][-1]
        
        #Molecular geometry
        Geometry = [data for label, data in X.Data if 'Geometry' in label][-1]
        
        #FED coupling matrix; will be initialized later
        FEDData = [data for label, data in X.Data if 'FED' in label]
        if len(FEDData) == 0:
            FED, FEDCouplings = None, None
        else:
            FED, FEDCouplings = FEDData[-1]

    ###################
    # Post-processing #
    ###################
    if doPlot:
        print 'Perceived the following state information:'
        for state in States:
            print state
            state.isValid()

    Coords = array([l[1:] for l in Geometry])

    if len(States) == 0:
        print filename, 'No data found'
        return False

    if True or FED is None: #Assume this is a monomer calculation
        if doPlot:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D

            print 'Plotting geometry in Figure 1'
            fig = plt.figure(1)
            ax = Axes3D(fig)
            for x1, y1, z1 in Coords/Angstrom:
                for x2, y2, z2 in Coords/Angstrom:
                    if 0 < ((x1-x2)**2 +  (y1-y2)**2 + (z1-z2)**2)**0.5 < 1.5:
                        ax.plot((x1,x2), (y1,y2), (z1,z2), 'k-', marker='.')
            plt.title('Geometry in Angstroms')
            #ax = fig.add_subplot(111, projection = '3d')
            c = Centroid(Coords)/Angstrom
            for state in States:
                tdip = state.TransitionDipole/Angstrom

                if state.OscillatorStrength > BrightStateStrengthThreshold:
                    d = c + tdip
                    ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]))
                    ax.text(d[0], d[1], d[2], str(state.Index))
            
            print 'Plotting absorption spectrum in Figure 2'
            fig = plt.figure(2)
            x = [HartreeToNm/state.ExcitationEnergy for state in States]
            y = [state.OscillatorStrength for state in States]
            plt.stem(x, y, markerfmt='*', linefmt='k-')
            plt.plot([min(x), max(x)], [BrightStateStrengthThreshold]*2, 'r:')
            plt.annotate('Bright state threshold', (0, BrightStateStrengthThreshold),
                xycoords=('axes fraction', 'data'), color='grey')
            plt.annotate('', xy=(1, BrightStateStrengthThreshold),
                xycoords=('axes fraction', 'data'), xytext=(0, BrightStateStrengthThreshold),
                arrowprops={'arrowstyle':'-', 'linestyle':'dashdot', 'color':'grey'} )
            for idx, state in enumerate(States):
               plt.text( x[idx], y[idx], str(state.Index), horizontalalignment='left')
                
            plt.xlabel('Wavelength (nm)')
            plt.ylabel('Oscillator strength')
            plt.title('Absorption spectrum predicted')

            #Plot couplings
            plt.show()

        from numpy.linalg import eig
        Inertia = MomentOfInertiaTensor(Coords)

        #Save monomer data
        f = open('monomer.pkl', 'wb')
        pickle.dump(Inertia, f)
        #Save bright states only
        States = filter(lambda x: x.OscillatorStrength > BrightStateStrengthThreshold, States)
        pickle.dump(len(States), f)
        for state in States:
            pickle.dump(state.Index, f)
            pickle.dump(state.ExcitationEnergy, f)
            pickle.dump(state.OscillatorStrength, f)
            pickle.dump(state.TransitionDipole, f)
        f.close()

        print 'Saved monomer data in monomer.pkl:'
        print 'Moment of inertia tensor:\n', Inertia
        print 'Transition dipoles:'
        for state in States:
            print state.Index, state.TransitionDipole
        return False
    else:
        #####################
        # DIMER CALCULATION #
        #####################
        #Load data from monomer pickle
        try:
            f = open('monomer.pkl', 'rb')
            #print 'Loading monomer data from monomer.pkl:'
            MonomerInertia = pickle.load(f)
            #print 'Moment of inertia tensor:\n', MonomerInertia
            numstates = pickle.load(f)
            MonomerData = []
            #print 'Transition dipoles:'
            for i in range(numstates):
                idx = pickle.load(f)
                energy = pickle.load(f)
                strength = pickle.load(f)
                tdip = pickle.load(f)
                MonomerData.append((idx, energy, strength, tdip))
            f.close()

        except (EOFError, IOError):
            raise ValueError, 'Tried to analyze FED calculation without monomer data in monomer.pkl.'

    # Step 1. Calculate geometric parameters

    # Find distance between centroids
    numatoms = len(Geometry)/2
    r1 = Centroid(Coords[:numatoms,:])
    r2 = Centroid(Coords[numatoms:,:])
    r = r1 - r2
    distance = norm(r)

    # Calculate local inertia tensors
    I1 = MomentOfInertiaTensor(Coords[:numatoms,:])
    I2 = MomentOfInertiaTensor(Coords[numatoms:,:])
   
    #Calculate rotation matrices for monomer 1 (M1) and 2 (M2)
    from numpy.linalg import svd, eig
    Um, Sm, Vm = svd(MonomerInertia)
    U1, S1, V1 = svd(I1)
    U2, S2, V2 = svd(I2)
    M1 = dot(U1, Vm)
    M2 = dot(U2, Vm)
    
    #Calculate all possible Forster couplings
    if doPlot:
        print 'Possible Forster couplings of monomer states:'

    import numpy.ma as ma
    num_monomer_states = max([x[0] for x in MonomerData])
    ForsterCouplings = ma.array(zeros((num_monomer_states, num_monomer_states)),
            mask=ones((num_monomer_states, num_monomer_states)))

    for idx, _, _, td in MonomerData:
        for idx2, _, _, td2 in MonomerData:
            c = abs(ForsterCoupling(dot(M1, td), dot(M2, td2), r))
            ForsterCouplings[idx-1, idx2-1] = c
            if doPlot and idx <= idx2:
                print 'M%d . M%d :' % (idx, idx2), c/eV, 'eV'
    
    #Prune dimer states too far or too low from monomer states
    #Limit margin to maximum displacement by strongest possible coupling with 10% fudge factor
    margin = 1.1*max(FEDCouplings.max(), -FEDCouplings.min(), \
                ForsterCouplings.max(), -ForsterCouplings.min())
    
    energies = [state[1] for state in MonomerData]
    min_nrg = min(energies) - margin
    max_nrg = max(energies) + margin
    
    if doPlot:
        print 'Matching dimer states in the energy range %.3f - %.3f eV' \
            % (min_nrg/eV, max_nrg/eV)
    
    """
    DiabaticHamiltonian = zeros((numstates, numstates))
    for (idx1, idx2), (change, coupling) in FEDCoupling.items():
        DiabaticHamiltonian[idx1, idx2] = coupling
        DiabaticHamiltonian[idx2, idx1] = coupling
    for state in States:
        idx = state.Index - 1
        DiabaticHamiltonian[idx, idx] = state.ExcitationEnergy
    
    from numpy.linalg import eig
    evals, evecs = eig(DiabaticHamiltonian)
    MatchingStates=dict()
    for i in range(numstates):
        evec = evecs[i,:]
        #Match by the two greatest components of the eigenvector
        lex = lexsort((evec, abs(evec)))
        #The first inequality corrects for duplicates
        #the second checks that there actually some coupling
        #the third checks that some FED density change actually occurred
        #the next two filters out states that fall outside a sensible energy
        if lex[-2] > lex[-1]: continue
        j, k = lex[-2:]
        print lex[-2:], FEDCoupling[j, k][0]
        if abs(FEDCoupling[j, k][0]) > 1e-6:
            MatchingStates[j, k] = FEDCoupling[j, k][0]
            print j, k, FEDCoupling[j, k][0]
    """
    #Research note: FEDCouplings is NOT the off-diagonal part of v. E. v.T or
    #v.T . E . v; I think they literally did a pairwise thing for all the
    #states. I have no idea what their diabatization scheme is. JC 2012-01-03
    
    #Match states by diagonalizing the FED matrix and picking out the coupled
    #states.
    u, v = eig(FED) #FED should have diagonal zeroed out
    MatchStates = []
    MatchCouplings = []
    lex_eval = lexsort((u, abs(u)))
    for eigval in u[lex_eval][::-1]:
        if abs(eigval) < ChargeChangeThreshold: break #done
        i = asscalar(array(where(u==eigval))) #some type munging :/
        evec = v[:,i]
        #Match by the two greatest components of the eigenvector
        lex = lexsort((evec, abs(evec)))
        if abs(evec[lex[-2]]) < 0.01: continue
        pair = sorted(lex[-2:])
        if pair in MatchStates: continue
        MatchStates.append(pair)
        #print 'Matched states::',pair
        MatchCouplings.append(FEDCouplings[tuple(pair)])
    MatchStates=array(MatchStates)
    MatchCouplings=array(MatchCouplings)
    lex = lexsort((MatchCouplings, abs(MatchCouplings)))
    matched_states = []
    MatchingStates = []
    for i, j in MatchStates[lex[::-1]]:
        if i not in matched_states and j not in matched_states \
                and min_nrg < States[i].ExcitationEnergy < max_nrg \
                and min_nrg < States[j].ExcitationEnergy < max_nrg:
            MatchingStates.append((i, j))
            matched_states.append(i)
            matched_states.append(j)
            if doPlot: print 'Matched states:',i,j
        elif i in matched_states or j in matched_states:
            if doPlot: print 'Warning: multiple states are coupled',i,j

    #Match by FED strength
    """
        j, k = lex[-2:].min(), lex[-2:].max()
        if j < k and \
                abs(DiabaticHamiltonian[j, k]) > 1e-6 and \
                abs(FEDCoupling[j, k][0]) > 1e-6 and \
                min_nrg < States[j].ExcitationEnergy < max_nrg and \
                min_nrg < States[k].ExcitationEnergy < max_nrg:

            #Check if already matched, and if it is, take the pair with greater
            #coupling strength
            if j in matched_states:
                the_matches = [matches for matches in MatchingStates if (matches[0]==j or matches[1]==j)]
                print the_matches

            matched_states.append(j)
            matched_states.append(k)
            matched_states = list(set(matched_states)) #uniqfy
            MatchingStates.append((j, k))
        if j < k and doPlot:
                print 'States',j+1,'and',k+1,'are coupled', \
                    'with strength', DiabaticHamiltonian[j,k]/eV, \
                    'eV and FED charge change', FEDCoupling[j,k][0], \
                    'and mutual eigenvector population', evec[lex[-2]]**2
        """
    #Match states by amount of change in electron and hole densities in excitation
    #This has the benefit of automatically removing spurious dark CT states which
    #contaminate TDDFT excitation spectra
    """
    numstates = len(States)
    FED_ChargeChanges = zeros((numstates, numstates))
    FED_Couplings = zeros((numstates, numstates))
    for (idx1, idx2), (change, coupling) in FEDCoupling.items():
        FED_ChargeChanges[idx1, idx2] = change
        FED_Couplings[idx1, idx2] = coupling
        FED_ChargeChanges[idx2, idx1] = change
        FED_Couplings[idx2, idx1] = coupling

    
    allowed_state_idxs = [state.Index-1 for state in States]
    matched_states = []
    MatchingStates = []
    for idx in allowed_state_idxs:
        if idx in matched_states: #Already paired
            continue
        x =FED_ChargeChanges[idx,:]
        xx=lexsort((x, abs(x)))
        x=x[xx] #Sorted by largest absolute value
        x=x[-1:] #XXX TAKE LARGEST ONE FOR NOW
        xx=xx[-1:]
        #Retrieve indices of matches
        matches = xx
        for idx2, change in zip(matches, x):
            if idx2 not in matched_states and \
                    abs(change) > ChargeChangeThreshold:  
                pair = min((idx, idx2)), max((idx, idx2))
                if doPlot:
                    print 'State', pair[0]+1, 'is exciton-coupled to', \
                        pair[1]+1, 'with density change', change
                matched_states.append(idx)
                matched_states.append(idx2)
                MatchingStates.append(pair)
    """
    #Match Forster couplings by comparing energies
    #The heuristic is to match the average of the dimer state energies
    #to the average of the monomer state energies
    #post-mortem 2012-01-02: This seems to work well most of the time
    #but sometimes assigns couplings to multiple dimer states.
    monomer_energies = {}
    for idx1, energy1, _, _ in MonomerData:
        for idx2, energy2, _, _ in MonomerData:
            if idx1 == idx2:
                mono_average = (energy1+energy2)/2
                monomer_energies[mono_average] = (idx1, idx2)
    
    DimerStates = States[:]
    for idx1, idx2 in MatchingStates:
        #Match by closest energy
        dimer_average = (DimerStates[idx1].ExcitationEnergy+\
            DimerStates[idx2].ExcitationEnergy)/2 
        monomer_average = nearest(dimer_average, monomer_energies)
        monomers = monomer_energies[monomer_average]
        TheForsterCoupling = ForsterCouplings[monomers[0]-1,monomers[1]-1]
        TheFEDCoupling = abs(FEDCouplings[idx1, idx2])
        #Matching using closest value 
        #TheForsterCoupling = nearest(TheFEDCouplings, ForsterCoupling.values())

        print filename, 
        print distance/Angstrom, 
        print idx1+1, idx2+1, TheFEDCoupling/eV,
        print TheForsterCoupling/eV,
        print 'M%d M%d' % monomers,
        td1 = [x[3] for x in MonomerData if x[0] == monomers[0]][0]
        td2 = [x[3] for x in MonomerData if x[0] == monomers[1]][0]
        print ForsterOrientationFactor(td1, td2, r),
        print dimer_average/eV, monomer_average/eV
    
    if doPlot:
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D

        print 'Plotting geometry in Figure 1'
        fig = plt.figure(1)
        #ax = fig.add_subplot(111, projection = '3d')
        ax = Axes3D(fig)
        for x1, y1, z1 in Coords/Angstrom:
            for x2, y2, z2 in Coords/Angstrom:
                if 0 < ((x1-x2)**2 +  (y1-y2)**2 + (z1-z2)**2)**0.5 < 1.5:
                    ax.plot((x1,x2), (y1,y2), (z1,z2), 'k-', marker='.')
        plt.title('Geometry in Angstroms')

        #Transition dipoles on first monomer
        c = r1/Angstrom
        for idx, _, strength, td in MonomerData:
            tdip = dot(M1, td/Angstrom)
            if strength > BrightStateStrengthThreshold: #XXX Don't need this test anymore if we assume threshold is unchanged
                d = c + tdip
                ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]), 'b-')
                ax.text(d[0], d[1], d[2], 'M'+str(idx))
        
        #Transition dipoles on second monomer
        c = r2/Angstrom
        for idx, _, strength, td in MonomerData:
            tdip = dot(M2, td/Angstrom)
            if strength > BrightStateStrengthThreshold:
                d = c + tdip
                ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]), 'r-')
                ax.text(d[0], d[1], d[2], 'M'+str(idx))

        #Transition dipoles on dimer
        c = Centroid(Coords)/Angstrom
        for state in DimerStates:
            tdip = state.TransitionDipole/Angstrom
            if state.OscillatorStrength > BrightStateStrengthThreshold:
                d = c + tdip
                ax.plot((c[0], d[0]), (c[1], d[1]), (c[2], d[2]), 'g-')
                ax.text(d[0], d[1], d[2], state.Index)
        
        print 'Plotting absorption spectrum in Figure 2'
        fig = plt.figure(2)

        #Data for dimer
        x = [HartreeToNm/state.ExcitationEnergy for state in DimerStates]
        y = [state.OscillatorStrength for state in DimerStates]
        plt.stem(x, y, markerfmt='*', linefmt='k-')
        plt.annotate('Bright', (0.02, BrightStateStrengthThreshold),
                xycoords=('axes fraction', 'data'), xytext=(0.02, BrightStateStrengthThreshold*1.3),
                arrowprops={'arrowstyle':'<-', 'relpos':(0,0), 'color':'grey'}, color='grey')
        plt.annotate('', xy=(1, BrightStateStrengthThreshold),
                xycoords=('axes fraction', 'data'), xytext=(0, BrightStateStrengthThreshold),
                arrowprops={'arrowstyle':'-', 'linestyle':'dashdot', 'color':'grey'} )
        for idx, state in enumerate(DimerStates):
           plt.text( x[idx], y[idx], str(state.Index), horizontalalignment='left')
       
        #Draw couplings
        for idx1, idx2 in MatchingStates:
            xs = [x[idx1], x[idx1], x[idx2], x[idx2]]
            barheight = max(y[idx1], y[idx2])*1.05
            ys = [y[idx1], barheight, barheight, y[idx2]]
            plt.plot(xs, ys, 'g:')
        #Data for monomer
        x = [HartreeToNm/state[1] for state in MonomerData]
        y = [state[2] for state in MonomerData]
        plt.stem(x, y, markerfmt='.', linefmt='b--')
        for idx, state in enumerate(MonomerData):
           plt.text(x[idx], y[idx], 'M'+str(state[0]), horizontalalignment='left')

        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Oscillator strength')
        plt.title('Absorption spectrum predicted')

        plt.show()




def nearest(x, things):
    y = sorted([(abs(z-x), z) for z in things])
    return y[0][1]



def ForsterCoupling(d1, d2, r):
    """
    This function calculates the Forster coupling between two chromophores
    using the dipole-dipole coupling approximation.

    @param d1 Transition dipole of 1st chromophore in atomic units (3-vector)
    @param d2 Transition dipole of 2nd chromophore in atomic units (3-vector)
    @param r  Displacement vector between the two chromophores in atomic units
    @returns The coupling matrix element in atomic units
    @f[
    V = \frac {3 (d_1 \cdot \hat r) (d_2 \cdot \hat r) - (d_1 \cdot d_2)} {\vert r \vert^3 }
    @f]

    @note the formula doesn't care which direction the displacemnent vector is in
    """
    normr = norm(r)
    rn = r / normr ##Normalized distance vector
    Coupling = (3 * dot(d1, rn) * dot(d2, rn) - dot(d1, d2)) / normr**3
    return Coupling

def ForsterOrientationFactor(d1, d2, r):
    """
    This function calculates the Forster orientation factor between two chromophores
    using the dipole-dipole coupling approximation.

    @param d1 Transition dipole of 1st chromophore in atomic units (3-vector)
    @param d2 Transition dipole of 2nd chromophore in atomic units (3-vector)
    @param r  Displacement vector between the two chromophores in atomic units
    @returns The coupling matrix element in atomic units
    @f[
    \kappa
    @f]

    @note the formula doesn't care which direction the displacemnent vector is in
    """
    rn =  r / norm(r) ##Normalized distance vector
    d1n = d1/ norm(d1)
    d2n = d2/ norm(d2)
    Factor = 3 * dot(d1n, rn) * dot(d2n, rn) - dot(d1n, d2n)
    return Factor

if __name__ == '__main__':
    import glob, sys
    if len(sys.argv) == 1:
        print """I will now proceed to process *.out in this directory.
If you did not want this behavior, abort now and rerun
with specified output file(s) on command line"""
        args = ['*.out']
    else:
        args = sys.argv[1:]

    #Get glob of globs
    filenames = []
    for fileglob in args:
        filenames += glob.glob(fileglob)

    #Uniqfy
    filenames = dict().fromkeys(filenames).keys()
    filenames.sort()

    doPlot = True if len(filenames) == 1 else False
    for fname in filenames:
        data = ParseFED(fname, doPlot)

        
