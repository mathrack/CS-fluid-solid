#!/usr/bin/env python
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys
from optparse import OptionParser
import math
from numpy.fft import fft
#from scipy import arange
#from pylab import plot, show, xlabel, ylabel, subplot, ylim
#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """
    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-c", "--cases", dest="cases", type="string",
                      help="String which contains the list of the cases")

    parser.add_option("-d", "--directories", dest="directories",
                      help="String which contains the list of the directories of result")

    parser.add_option("-s", "--study", dest="study", type="string",
                      help="label of the current study")

    parser.add_option("-e", "--exp", dest="exp", type="string",
                      help="Name of the configuration")

    (options, args) = parser.parse_args(argv)

    return options
#-------------------------------------------------------------------------------

def main(options):
# data for shock tube problem

# File reading and writing for one configuration (SOD, DDS, DCS or DCI)
# we read L1 error for each mesh and write it in a single file

    ny = []
    e1 = []
    e2 = []
    ei = []
    e1grad = []
    e2grad = []
    eigrad = []

    runs = options.directories
    runs = runs.split()
    study = os.path.join(runs[0], '../../..')
    study = os.path.normpath(study)
    post = os.path.join(study, 'POST')

    n = 0
    for i in range(len(runs)):
        case = os.path.join(runs[i], '../..')
        case = os.path.normpath(case)
        head, case = os.path.split(case)

        if case == options.exp:
            n += 1
            ff = os.path.join(runs[i], 'error.dat')
            f = open(ff)
            flines = f.readlines()
            f.close()

            for l in flines:
                if l.rfind(" #") != 0:
                    d = l.split()
                    ny.append(float(d[0]))
                    e1.append(float(d[1]))
                    e2.append(float(d[2]))
                    ei.append(float(d[3]))

            ff = os.path.join(runs[i], 'error_grad.dat')
            f = open(ff)
            flines = f.readlines()
            f.close()

            for l in flines:
                if l.rfind(" #") != 0:
                    d = l.split()
                    e1grad.append(float(d[1]))
                    e2grad.append(float(d[2]))
                    eigrad.append(float(d[3]))

    out_name = options.exp + "_hexa.tex"
    out_f = os.path.join(post, out_name)
    ncells = open(out_f, 'w')
    for m in range(n/2):
        ncells.write( '$ ' + '%7.4e' %(ny[m]) + ' $ & $ ' + '%7.4e' %(e2[m]) + ' $ & $ ' + '%7.4e' %(e2grad[m]) + ' $ \\\\ \n' )
    ncells.close()
    out_name = options.exp + "_tetra.tex"
    out_f = os.path.join(post, out_name)
    ncells = open(out_f, 'w')
    for m in range(n/2,n):
        ncells.write( '$ ' + '%7.4e' %(ny[m]) + ' $ & $ ' + '%7.4e' %(e2[m]) + ' $ & $ ' + '%7.4e' %(e2grad[m]) + ' $ \\\\ \n' )
    ncells.close()

    out_name = options.exp + "_errors.dat"
    out_f = os.path.join(post, out_name)
    out = open(out_f, 'w')
    out.write('\n')
    out.write(' ======\n')
    out.write(' Meshtype '+'|'+' Erreur L1 '+'|'+' Erreur L2 '+'|'+' Erreur Linf \n')
    out.write(' ======\n')
    out.write('\n')
    for m in range(n):
        if m == 0:
            out.write( ' hexa.05 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
        if m == 1:
            orde1 = math.log(e1[m-1]/e1[m])/math.log(ny[m]/ny[m-1])
            orde2 = math.log(e2[m-1]/e2[m])/math.log(ny[m]/ny[m-1])
            ordei = math.log(ei[m-1]/ei[m])/math.log(ny[m]/ny[m-1])
            out.write( ' hexa.10 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1)+' | '+str(orde2)+' | '+str(ordei)+' |\n')
        if m == 2:
            orde1 = math.log(e1[m-1]/e1[m])/math.log(ny[m]/ny[m-1])
            orde2 = math.log(e2[m-1]/e2[m])/math.log(ny[m]/ny[m-1])
            ordei = math.log(ei[m-1]/ei[m])/math.log(ny[m]/ny[m-1])
            out.write( ' hexa.15 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1)+' | '+str(orde2)+' | '+str(ordei)+' |\n')
        if m == 3:
            orde1 = math.log(e1[m-1]/e1[m])/math.log(ny[m]/ny[m-1])
            orde2 = math.log(e2[m-1]/e2[m])/math.log(ny[m]/ny[m-1])
            ordei = math.log(ei[m-1]/ei[m])/math.log(ny[m]/ny[m-1])
            out.write( ' hexa.20 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1)+' | '+str(orde2)+' | '+str(ordei)+' |\n')
        if m == 4:
            out.write( ' tetr.05 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
        if m == 5:
            orde1 = math.log(e1[m-1]/e1[m])/math.log(ny[m]/ny[m-1])
            orde2 = math.log(e2[m-1]/e2[m])/math.log(ny[m]/ny[m-1])
            ordei = math.log(ei[m-1]/ei[m])/math.log(ny[m]/ny[m-1])
            out.write( ' tetr.10 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1)+' | '+str(orde2)+' | '+str(ordei)+' |\n')
        if m == 6:
            orde1 = math.log(e1[m-1]/e1[m])/math.log(ny[m]/ny[m-1])
            orde2 = math.log(e2[m-1]/e2[m])/math.log(ny[m]/ny[m-1])
            ordei = math.log(ei[m-1]/ei[m])/math.log(ny[m]/ny[m-1])
            out.write( ' tetr.15 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1)+' | '+str(orde2)+' | '+str(ordei)+' |\n')
        if m == 7:
            orde1 = math.log(e1[m-1]/e1[m])/math.log(ny[m]/ny[m-1])
            orde2 = math.log(e2[m-1]/e2[m])/math.log(ny[m]/ny[m-1])
            ordei = math.log(ei[m-1]/ei[m])/math.log(ny[m]/ny[m-1])
            out.write( ' tetr.20 '+' | '+repr(e1[m])+' | '+repr(e2[m])+' | '+repr(ei[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1)+' | '+str(orde2)+' | '+str(ordei)+' |\n')
    out.close()

    out_name = options.exp + "_errors_grad.dat"
    out_f = os.path.join(post, out_name)
    out = open(out_f, 'w')
    out.write('\n')
    out.write(' ======\n')
    out.write(' Meshtype '+'|'+' Erreur L1 '+'|'+' Erreur L2 '+'|'+' Erreur Linf \n')
    out.write(' ======\n')
    out.write('\n')
    for m in range(n):
        if m == 0:
            out.write( ' hexa.05 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
        if m == 1:
            orde1grad = math.log(e1grad[m-1]/e1grad[m])/math.log(ny[m]/ny[m-1])
            orde2grad = math.log(e2grad[m-1]/e2grad[m])/math.log(ny[m]/ny[m-1])
            ordeigrad = math.log(eigrad[m-1]/eigrad[m])/math.log(ny[m]/ny[m-1])
            out.write( ' hexa.10 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1grad)+' | '+str(orde2grad)+' | '+str(ordeigrad)+' |\n')
        if m == 2:
            orde1grad = math.log(e1grad[m-1]/e1grad[m])/math.log(ny[m]/ny[m-1])
            orde2grad = math.log(e2grad[m-1]/e2grad[m])/math.log(ny[m]/ny[m-1])
            ordeigrad = math.log(eigrad[m-1]/eigrad[m])/math.log(ny[m]/ny[m-1])
            out.write( ' hexa.15 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1grad)+' | '+str(orde2grad)+' | '+str(ordeigrad)+' |\n')
        if m == 3:
            orde1grad = math.log(e1grad[m-1]/e1grad[m])/math.log(ny[m]/ny[m-1])
            orde2grad = math.log(e2grad[m-1]/e2grad[m])/math.log(ny[m]/ny[m-1])
            ordeigrad = math.log(eigrad[m-1]/eigrad[m])/math.log(ny[m]/ny[m-1])
            out.write( ' hexa.20 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1grad)+' | '+str(orde2grad)+' | '+str(ordeigrad)+' |\n')
        if m == 4:
            out.write( ' tetr.05 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
        if m == 5:
            orde1grad = math.log(e1grad[m-1]/e1grad[m])/math.log(ny[m]/ny[m-1])
            orde2grad = math.log(e2grad[m-1]/e2grad[m])/math.log(ny[m]/ny[m-1])
            ordeigrad = math.log(eigrad[m-1]/eigrad[m])/math.log(ny[m]/ny[m-1])
            out.write( ' tetr.10 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1grad)+' | '+str(orde2grad)+' | '+str(ordeigrad)+' |\n')
        if m == 6:
            orde1grad = math.log(e1grad[m-1]/e1grad[m])/math.log(ny[m]/ny[m-1])
            orde2grad = math.log(e2grad[m-1]/e2grad[m])/math.log(ny[m]/ny[m-1])
            ordeigrad = math.log(eigrad[m-1]/eigrad[m])/math.log(ny[m]/ny[m-1])
            out.write( ' tetr.15 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1grad)+' | '+str(orde2grad)+' | '+str(ordeigrad)+' |\n')
        if m == 7:
            orde1grad = math.log(e1grad[m-1]/e1grad[m])/math.log(ny[m]/ny[m-1])
            orde2grad = math.log(e2grad[m-1]/e2grad[m])/math.log(ny[m]/ny[m-1])
            ordeigrad = math.log(eigrad[m-1]/eigrad[m])/math.log(ny[m]/ny[m-1])
            out.write( ' tetr.20 '+' | '+repr(e1grad[m])+' | '+repr(e2grad[m])+' | '+repr(eigrad[m])+' |\n')
            out.write( '   ordre '+' | '+str(orde1grad)+' | '+str(orde2grad)+' | '+str(ordeigrad)+' |\n')
    out.close()

    out_name = options.exp + ".dat"
    out_f = os.path.join(post, out_name)
    out = open(out_f, 'w')
    for m in range(n):
        out.write( repr(ny[m]) + ' ' + repr(e2[m]) + ' ' + repr(e2grad[m]) + ' \n')
    out.close()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    options = process_cmd_line(sys.argv[1:])
    main(options)

#-------------------------------------------------------------------------------
