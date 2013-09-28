#! /usr/bin/env python

import os
import sys
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> recover
import subprocess
import fileinput
import string

<<<<<<< HEAD
filePath = '/tthome/naodell/storage/data/nuTuples_v7_4'
sampleList = [
    'DYJets',
    'DYJets_M-10To50',
    'WJetsToLNu',
    'ZbbToLL',
    'WbbToLL',

    'DoubleMu_Run2012A',
    'DoubleMu_Run2012B',
    'DoubleMu_Run2012C',
    'DoubleMu_Run2012D',
    'DoubleElectron_Run2012A',
    'DoubleElectron_Run2012B',
    'DoubleElectron_Run2012C',
    'DoubleElectron_Run2012D',
    'MuEG_Run2012A',
    'MuEG_Run2012B',
    'MuEG_Run2012C',
    'MuEG_Run2012D',

    #'TTJets',
    #'T_tW',
    #'Tbar_tW',
    #'T_t',
    #'Tbar_t',

    #'WWGJets',
    #'WWWJets',
    #'ZZZNoGstarJets',
    #'WZZNoGstarJets',
    #'WWZNoGstarJets',

    #'ZZTo4e',
    #'ZZTo4mu',
    #'ZZTo4tau',
    #'ZZTo2e2mu',
    #'ZZTo2e2tau',
    #'ZZTo2mu2tau',
    #'WGstarToLNu2E',
    #'GluGluToZZTo2L2L',
    #'GluGluToZZTo4L',

    #'TTZJets',
    #'TTWJets',
    #'TTGJets',

    #'ZZJetsTo2L2Nu',
    #'ZZJetsTo2L2Q',
    #'WWJetsTo2L2Nu',
    #'WZJetsTo3LNu',
    #'WZJetsTo2L2Q',

    #'QCD_Pt_20_MuEnrichedPt_15',
    #'QCD_Pt_20_30_EMEnriched',
    #'QCD_Pt_30_80_EMEnriched',
    #'QCD_Pt_80_170_EMEnriched',
    #'QCD_Pt_170_250_EMEnriched',
    #'QCD_Pt_250_350_EMEnriched',
    #'QCD_Pt_350_EMEnriched',

    #'GluGluToHToWWTo2LAndTau2Nu_M-125',
    #'GluGluToHToZZTo4L_M-125',
    #'WH_HToWW_3l_M-125',
    #'TTH_Inclusive_M-125'
]

for sample in sampleList:

    path = '{0}/{1}'.format(filePath, sample)
    if not os.path.exists(path):
        os.system('mkdir -p {0}'.format(path))
    elif len(os.listdir(filePath)) is not 0 and False:
        os.system('rm -r {0}'.format(path))

    files = [line.split(' ')[-1] for line in os.popen('srmls -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/naodell/nuTuples_v7_4/{0}'.format(sample)).read().split('\n')]

    print files
    print 'Transferring {0} files from {1} dataset from lpc...'.format(len(files), sample)

    for i,dataset in enumerate(files):
        if dataset == '' : continue
        print 'Transferring file {0}...'.format(i+1, len(files))
        os.system('srmcp -num_streams=10 -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN={0} file:////tthome/naodell/storage/data/nuTuples_v7_4/{1}/'.format(dataset, sample))

    #os.system('srmcp -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/naodell/nuTuples_v7_4/{0}/ file:////tthome/naodell/storage/data/nuTuples_v7_4/{0}/. -recursive'.format(sample))
=======
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
>>>>>>> 18ffb5e611afdfafac090643e0c149608c892737
=======
inPath  = '/11/store/user/naodell/nuTuples_v7_4'
outPath = '/tthome/naodell/storage/data/nuTuples_v7_4'

sampleList = [line.split(' ')[-1] for line in os.popen('srmls -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN={0}'.format(inPath)).read().split('\n')]

for sample in sampleList[1:]:
    if sample == '': continue

    path = '{0}/{1}'.format(outPath, sample)
    if not os.path.exists(path):
        os.system('mkdir -p {0}'.format(path))
    elif len(os.listdir(outPath)) is not 0 and False:
        os.system('rm -r {0}'.format(path))

    files = [line.split(' ')[-1] for line in os.popen('srmls -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN={0}'.format(sample)).read().split('\n')]

    print 'Transferring {0} files from {1}...'.format(len(files), sample)

    for i,dataset in enumerate(files[1:]):
        if dataset == '' : continue

        print 'Transferring file {0}...'.format(i+1, len(files))
        os.system('srmcp -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN={0} file:///{1}/{2}/'.format(dataset, outPath, sample))

    #os.system('srmcp -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN=/11/store/user/naodell/nuTuples_v7_4/{0}/ file:////tthome/naodell/storage/data/nuTuples_v7_4/{0}/. -recursive'.format(sample))
>>>>>>> recover
