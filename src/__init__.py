#!/usr/bin/env python

try:
    import tables
except ImportError:
    log.info('PyTables not found. HDF5 archival routines will be disabled.')

