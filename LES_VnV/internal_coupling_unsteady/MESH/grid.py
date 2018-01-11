#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys
import string
from optparse import OptionParser

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------                       
# Application modules import
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """Processes the passed command line arguments."""
    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-c", "--case", dest="case", type="string",
                      help="Directory of the current case")

    parser.add_option("-p", "--param", dest="param", type="string",
                      help="Name f the file of parameters")

    parser.add_option("-m", "--mesh", dest="mesh", type="string",
                      help="Name of the new mesh")

    (options, args) = parser.parse_args(argv)

    return options

#-------------------------------------------------------------------------------

def main(options):
    fp = os.path.join(options.case, "DATA", options.param)
    f = open(fp, 'r');
    lines = f.readlines()
    f.close()
    os.remove(fp)
    f = open(fp, 'w');
    for l in lines:
        if l.rfind("#") < 0:
            if (l.rfind('domain.mesh_input') > 0):
                l = "        domain.mesh_input = '"+str(options.mesh)+"'\n"
        f.write(l)
    f.close()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    options = process_cmd_line(sys.argv[1:])
    main(options)

#-------------------------------------------------------------------------------
