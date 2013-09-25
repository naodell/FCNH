#! /usr/bin/env python

import os
import sys
import fileinput
import string

folders = [
   "DYJets",
]

for line in folders:
    #os.system("mkdir /eos/uscms/store/user/akub19/TREES_SEP_2013/"+line)
    os.system("srmls  -2  srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/naodell/nuTuples_v7_4/"+line+" > files_"+line+".log") 
    files = open("files_"+line+".log")
    for file in files:
      if file.find("root") > 0:
        file = file.rstrip("\n")
        things = file.split(" ")
        #print "./"+line,things[-1]
        os.system("srmcp -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN="+things[-1]+" file:////tthome/naodell/storage/data/nuTuples_v7_4/"+line+"/")
    files.close()
    os.system("rm files_"+line+".log")
