#!/usr/bin/env python

def dedrude(infilename, outfilename):
    "Removes Drude particles from CHARMMM coordinate file"
    buf = []    
    state = 'header'
    for l in open(infilename):
        t = l.split()
        if l[0] == '*':
            buf.append(l)
        elif state == 'header':
            if t[1] == 'EXT':
                numatomslinenum = len(buf) - 1
                buf.append(l)
                atom = 0
                state = 'atom'
        elif state == 'atom':
            if t[3][0] != 'D': #If not Drude
                atom += 1
                buf.append('%10d' % atom)
                buf.append(l[10:])
    buf[numatomslinenum] = '%10d' % atom + '  EXT'

    f = open(outfilename, 'w')
    f.write(''.join(buf))
    f.close()

if __name__ == '__main__':
    import sys
    dedrude(sys.argv[1], sys.argv[2])
