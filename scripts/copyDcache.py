#! /usr/bin/env python

import os
import sys
import string

inPath  = '/11/store/user/naodell/nuTuples_v7_4'
outPath = '/tthome/naodell/storage/data/nuTuples_v7_4'

sampleList = [line.split(' ')[-1] for line in os.popen('srmls -2 srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN={0}'.format(inPath)).read().split('\n')]

for sample in sampleList[1:]:
    os.system('srm-copy srm://cmssrm.fnal.gov:8443/srm/managerv2?SFN={0} -td file:///{1}/. -recursive'.format(sample, outPath))
