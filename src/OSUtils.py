#!/usr/bin/env python

"""
Miscellaneous OS Utilities
"""

from fnmatch import fnmatch
from glob import glob
from heapq import heappop, heappush
import logging, os, shutil

def chdirn(thedir):
    """
    Changes into a directory. If directory does not exist, make it.

    :param string thedir: The name of the directory to go to.
    :returns: None 
    """
    try:
        os.mkdir(thedir)
    except OSError: #Directory exists
        pass
    os.chdir(thedir)



def find_files(directory, pattern):
    """
    Iterator that recursively traverses directory tree *in order*
    for files matching pattern
    
    .. NOTE:: does NOT recurse into symbolic links.

    :param string directory: Name of directory to traverse
    :param string pattern: File pattern to match

    :returns: Iterator
    :rtype: tuple(current_directory, filename)

    .. versionadded:: 0.1
    """

    stack = [directory] #A priority queue
    while stack:
        nowdir = heappop(stack)
        for base in os.listdir(nowdir):
            name = os.path.normpath(os.path.join(nowdir, base))
            if os.path.isdir(name):
                if not os.path.islink(name): #Avoid infinite loop with symlinks
                    heappush(stack, name)
                    #Use default lexicographical order for priority (sort)
            elif fnmatch(name, pattern):
                yield nowdir, name



def WildCardExpandedFileList(iterable):
    """Takes an iterable of filenames or wildcarded selections parseable by
    glob and returns a uniq-ed, sorted list of filenames"""

    filenames = set()
    for blob in iterable:
        filenames = filenames.union(set(glob(blob)))
    return sorted(filenames)


def copywild(src, dest):
    """Wildcard copy.
    NOT recursive."""

    logger = logging.getLogger('OSUtils.copywild')

    for filename in glob(src):
        if os.path.isfile(filename):
            try:
                shutil.copy(filename, os.path.join(src, dest))
            except IOError, e:
                logger.error('Encountered error while copying file %s: %s', 
                    filename, e)


def linkwild(src, dest, overwrite=False):
    """
    Symbolic link creation with wildcard support.
    
    @note Symlinks matched by src will *not* be copied to avoid possible
    circular references.

    @param src A specification of source files parseable by glob.glob()

    @param dest A destination (directory).

    @param overwrite: Whether to overwrite files that exist in :py:var:dest.
    """

    logger = logging.getLogger('linkwild')

    for fname in glob(src):
        destfname = os.path.join(dest, os.path.basename(fname))
        if overwrite:
            doit = True
        else:
            doit = not os.path.exists(destfname)
        if (os.path.isfile(fname) or os.path.islink(fname)) and doit:
            if overwrite and os.path.exists(destfname):
                os.unlink(destfname)
            try:
                os.symlink(fname, destfname)
                logger.info('Made symlink %s --> %s', os.path.abspath(fname),
                        os.path.abspath(destfname))
            except OSError:
                logger.warning('Warning: could not make symlink: %s -X-> %s',
                        os.path.abspath(fname), os.path.abspath(destfname))

