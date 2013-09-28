#! /usr/bin/env python

import os
import sys
import subprocess
import fileinput
import string

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
