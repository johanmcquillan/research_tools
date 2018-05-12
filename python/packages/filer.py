"""Methods for creating nested directories."""

import os
import errno


def make_directory(path):
    """Make a directory, creating any non-existent intermediate directories as necessary,"""
    
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def safe_open(filename, read_write):
    """Open a file and create subdirectories as necessary.
    
    Args:
        filename (str): Path to file
        read_write (str): Read/write option (same as used in built-in open())
    """
    
    if os.path.dirname(filename) != '':
        make_directory(os.path.dirname(filename))
    return open(filename, read_write)

