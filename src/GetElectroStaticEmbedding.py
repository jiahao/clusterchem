

def LoadGromacsTopology(filename):
    """
    Parses Gromacs topology file for atomic charge definitions
    @param filename Name of Gromacs topology file (usually .top or .itp)
    to parse.
    @returns a dictionary of atomtypes with (residue, atomtype) as key
    and charge as value.
    """
    AtomTypes = {}
    mode = 'seek'
    for line in open(filename):
        if mode == 'seek' and '[ atoms ]' in line:
            mode = 'read'
        elif mode == 'read':
            #The GROMACS topology file format specifies lines to have the form
            #;   nr    type   resnr  residu    atom    cgnr  charge
            theline = line[:line.find(';')] #Ignore comments between semicolon and newline
            assert '\\' not in theline, 'Cannot handle line continuations at this time.'
            t = theline.split()
            try:
                residu = t[3]
                atom = t[4]
                charge = float(t[6])
            except (ValueError, IndexError):
                #Could not parse, assume we are done with this section
                mode = 'seek'
                continue
            AtomTypes[residu, atom] = charge
    return AtomTypes


def LoadGromacsGeometry(h5table, filename):
    """
    Parses Gromacs topology file for atomic charge definitions
    @param h5table HDF5 table to populate. Must be in CHARMM_CARD format.
    @param filename Name of Gromacs topology file (usually .top or .itp)
    to parse.
    """
    mode = 'title'
    for line in open(filename):
        if mode == 'title':
            mode = 'numatoms'
        elif mode == 'numatoms':
            numatoms = int(line)
            thisnumatoms = 0
            mode = 'read'
        elif mode == 'read':
            try:
                data = h5table.row
                #The GROMACS format does not contain ResID and SegId
                #We fill them in with ResidNo.
                data['ResID'] = data['SegID'] = data['ResidNo'] = int(line[:5])
                data['Res'] = line[5:10].strip()
                data['Type'] = line[10:15].strip()
                data['AtomNo'] = int(line[15:20])
                numbers = map(float, line[20:].split())
                data['Coord'] = numbers[:3] #Discard velocities if present
                #The GROMACS format does not contain Weighting
                #Set to dummy value
                data['Weighting'] = 0
                data.append()
                thisnumatoms += 1
            except (ValueError, IndexError): #Assume this is the last line with box vectors
                #boxvectors = map(float, line.split())
                break
    assert thisnumatoms == numatoms, 'Wrong number of atoms read: expected %d but read %d' % (numatoms, thisnumatoms)
    h5table.flush()

