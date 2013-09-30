#! /usr/bin/env python
import subprocess, shlex, time, pickle, math, sys, os
from array import array
import ROOT as r


if __name__ == '__main__':

    ### Get command line arguements
    if len(sys.argv) > 1:
        batch = sys.argv[1]
    else:
        print 'Must provide information about input file!'
        exit()

    # Logon not automatically loaded through PyROOT (logon loads TMVA library) load also GUI
    r.gROOT.SetMacroPath("${ROOTSYS}/tmva/test/.") 
    r.gROOT.Macro       ("${ROOTSYS}/tmva/test/TMVAlogon.C")
    r.gROOT.LoadMacro   ("${ROOTSYS}/tmva/test/TMVAGui.C")
    r.gROOT.ProcessLine('TMVAGui(\"mvaOutput/{0}.root\")'.format(batch))
    r.gApplication.Run() 
