#!/usr/bin/env python

"""
Miscellaneous OS Utilities

"""
import os


def chdirn(thedir):
    """
    Changes into a directory. If directory does not exist, make it.

    :param string thedir: The name of the directory to go to.
    :returns: None 
    """
    if not os.path.exists(thedir):
        os.mkdir(thedir)
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
    from heapq import heappop, heappush
    from fnmatch import fnmatch

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

    from glob import glob
    filenames = set()
    for blob in iterable:
        filenames = filenames.union(set(glob(blob)))
    return sorted(filenames)
    

